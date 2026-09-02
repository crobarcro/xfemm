#include "AnalysisSession.h"
#include "FSolverAnalysisBackend.h"
#include "FemmReader.h"
#include "TangleMesherBackend.h"
#include "TriangleMesherBackend.h"
#include "fpproc.h"
#include "mesh/SolverMeshValidator.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using Periodicity = femm::mesh::SolverMesh::Periodicity;

struct Point { double x; double y; };
struct Tolerance { double absolute; double relative; };
struct Fixture {
    std::string name;
    std::string path;
    std::vector<Point> samples;
    Periodicity periodicity;
    bool hasPeriodicConstraints;
    Tolerance potentialTolerance;
    Tolerance fieldTolerance;
    Tolerance integralTolerance;
};
struct Run {
    std::unique_ptr<femm::AnalysisSession> session;
    std::shared_ptr<femm::FSolverAnalysisBackend> solver;
    std::unique_ptr<FPProc> postprocessor;
};

std::unique_ptr<femm::FemmProblem> readProblem(const std::string &path)
{
    auto problem = std::make_unique<femm::FemmProblem>(femm::FileType::MagneticsFile);
    auto parserView = std::shared_ptr<femm::FemmProblem>(problem.get(), [](femm::FemmProblem *) {});
    femm::MagneticsReader reader(parserView, std::cerr);
    if (reader.parse(path) != femm::F_FILE_OK)
        throw std::runtime_error("could not parse " + path);
    return problem;
}

Run run(const Fixture &fixture, const std::string &backendName)
{
    std::shared_ptr<fmesher::MesherBackend> mesher;
    if (backendName == "Triangle")
        mesher = std::make_shared<fmesher::TriangleMesherBackend>();
    else
        mesher = std::make_shared<fmesher::TangleMesherBackend>();
    auto solver = std::make_shared<femm::FSolverAnalysisBackend>();
    auto session = std::make_unique<femm::AnalysisSession>(
        femm::ModelDefinition(readProblem(fixture.path)), mesher, solver);
    const auto mesh = session->ensureMesh();
    const auto validation = femm::mesh::validateSolverMesh(*mesh);
    if (!validation.valid())
        throw std::runtime_error(fixture.name + "/" + backendName + " produced invalid mesh");
    if (mesh->nodes.empty() || mesh->elements.empty() || mesh->edges.empty())
        throw std::runtime_error(fixture.name + "/" + backendName + " produced empty topology");
    if (!mesh->airGaps.empty())
        throw std::runtime_error(fixture.name + "/" + backendName + " unexpectedly produced AGE topology");
    if (mesh->periodicConstraints.empty() != !fixture.hasPeriodicConstraints)
        throw std::runtime_error(fixture.name + "/" + backendName + " has wrong PBC presence");
    std::vector<unsigned> pairedNodeUses(mesh->nodes.size(), 0);
    constexpr double boundaryTolerance = 1e-11;
    for (const auto &constraint : mesh->periodicConstraints) {
        if (constraint.periodicity != fixture.periodicity)
            throw std::runtime_error(fixture.name + "/" + backendName + " has wrong PBC type");
        const auto &first = mesh->nodes[constraint.first];
        const auto &second = mesh->nodes[constraint.second];
        const bool firstLowerSecondUpper =
            std::abs(first.y + 1.0) <= boundaryTolerance &&
            std::abs(second.y - 1.0) <= boundaryTolerance;
        const bool firstUpperSecondLower =
            std::abs(first.y - 1.0) <= boundaryTolerance &&
            std::abs(second.y + 1.0) <= boundaryTolerance;
        if ((!firstLowerSecondUpper && !firstUpperSecondLower) ||
            std::abs(first.x - second.x) > boundaryTolerance)
            throw std::runtime_error(fixture.name + "/" + backendName +
                                     " has incorrect paired-boundary correspondence");
        ++pairedNodeUses[constraint.first];
        ++pairedNodeUses[constraint.second];
    }
    if (fixture.hasPeriodicConstraints) {
        for (std::size_t i = 0; i < mesh->nodes.size(); ++i) {
            const auto &node = mesh->nodes[i];
            const bool onPairedBoundary =
                (std::abs(node.y + 1.0) <= boundaryTolerance ||
                 std::abs(node.y - 1.0) <= boundaryTolerance) &&
                std::abs(node.x) < 1.0 - boundaryTolerance;
            if (onPairedBoundary && pairedNodeUses[i] != 1)
                throw std::runtime_error(fixture.name + "/" + backendName +
                                         " does not pair every non-corner seam node exactly once");
        }
    }

    session->solve();
    auto postprocessor = std::make_unique<FPProc>();
    if (!postprocessor->OpenDocument(session->model().problem(), solver->solvedSolver(),
                                     solver->solvedSystem()))
        throw std::runtime_error(fixture.name + "/" + backendName + " post-processing failed");
    for (auto &label : postprocessor->blocklist)
        label.IsSelected = true;
    return {std::move(session), std::move(solver), std::move(postprocessor)};
}

bool finiteComparison(double triangle, double tangle, Tolerance tolerance,
                      double &difference, double &permitted)
{
    if (!std::isfinite(triangle) || !std::isfinite(tangle))
        return false;
    difference = std::abs(triangle - tangle);
    permitted = tolerance.absolute + tolerance.relative *
        std::max(std::abs(triangle), std::abs(tangle));
    return std::isfinite(difference) && std::isfinite(permitted);
}

void compare(const Fixture &fixture, const Point &point, const char *quantity,
             double triangle, double tangle, Tolerance tolerance)
{
    double difference = 0.0;
    double permitted = 0.0;
    if (!finiteComparison(triangle, tangle, tolerance, difference, permitted)) {
        std::cerr << fixture.name << " (" << point.x << ", " << point.y << ") "
                  << quantity << ": non-finite comparison: Triangle=" << triangle
                  << " Tangle=" << tangle << '\n';
        throw std::runtime_error(fixture.name + " " + quantity + " is non-finite");
    }
    std::cout << fixture.name << " (" << point.x << ", " << point.y << ") "
              << quantity << ": Triangle=" << triangle << " Tangle=" << tangle
              << " abs-difference=" << difference << " tolerance=" << permitted << '\n';
    if (difference > permitted)
        throw std::runtime_error(fixture.name + " " + quantity + " differential exceeded tolerance");
}

void compareFixture(const Fixture &fixture)
{
    auto triangle = run(fixture, "Triangle");
    auto tangle = run(fixture, "Tangle");
    double triangleMaxB = 0.0;
    double tangleMaxB = 0.0;
    for (const auto &point : fixture.samples) {
        CMPointVals tv, gv;
        if (!triangle.postprocessor->GetPointValues(point.x, point.y, tv) ||
            !tangle.postprocessor->GetPointValues(point.x, point.y, gv))
            throw std::runtime_error(fixture.name + " sample is outside a generated mesh");
        compare(fixture, point, "A", tv.A.re, gv.A.re, fixture.potentialTolerance);
        compare(fixture, point, "B1", tv.B1.re, gv.B1.re, fixture.fieldTolerance);
        compare(fixture, point, "B2", tv.B2.re, gv.B2.re, fixture.fieldTolerance);
        triangleMaxB = std::max(triangleMaxB, std::hypot(tv.B1.re, tv.B2.re));
        tangleMaxB = std::max(tangleMaxB, std::hypot(gv.B1.re, gv.B2.re));
    }
    const Point global{0, 0};
    const double triangleEnergy = triangle.postprocessor->BlockIntegral(2).re;
    const double tangleEnergy = tangle.postprocessor->BlockIntegral(2).re;
    if (!std::isfinite(triangleEnergy) || !std::isfinite(tangleEnergy) ||
        triangleEnergy <= 1e-12 || tangleEnergy <= 1e-12 ||
        !std::isfinite(triangleMaxB) || !std::isfinite(tangleMaxB) ||
        triangleMaxB <= 1e-9 || tangleMaxB <= 1e-9)
        throw std::runtime_error(fixture.name + " produced a non-finite or trivial solution");
    compare(fixture, global, "energy", triangleEnergy, tangleEnergy,
            fixture.integralTolerance);
    compare(fixture, global, "coenergy", triangle.postprocessor->BlockIntegral(17).re,
            tangle.postprocessor->BlockIntegral(17).re, fixture.integralTolerance);
}

} // namespace

int main(int argc, char **argv)
{
    if (argc != 4) {
        std::cerr << "expected non-periodic, periodic, and antiperiodic fixture paths\n";
        return 2;
    }
    try {
        double difference = 0.0;
        double permitted = 0.0;
        const Tolerance finiteProbe{1.0, 1.0};
        if (finiteComparison(std::numeric_limits<double>::quiet_NaN(), 0.0,
                             finiteProbe, difference, permitted) ||
            finiteComparison(std::numeric_limits<double>::infinity(), 0.0,
                             finiteProbe, difference, permitted) ||
            finiteComparison(std::numeric_limits<double>::max(),
                             -std::numeric_limits<double>::max(), finiteProbe,
                             difference, permitted) ||
            finiteComparison(1.0, 1.0,
                             {std::numeric_limits<double>::max(),
                              std::numeric_limits<double>::max()},
                             difference, permitted))
            throw std::runtime_error("finite comparison guard accepted a non-finite case");

        // The off-centre current rectangle breaks both x and y symmetry. Samples
        // stay at least 0.15 m from its interface and 0.25 m from the outer seams.
        const std::vector<Point> samples{{-0.60, -0.55}, {-0.55, 0.40},
                                         {-0.20, 0.20}, {0.70, 0.45},
                                         {0.275, -0.25}};
        const std::vector<Fixture> fixtures{
            {"non-periodic", argv[1], samples, Periodicity::Periodic, false,
             {5e-5, 0.002}, {1.5e-3, 0.02}, {0.1, 0.0025}},
            {"periodic", argv[2], samples, Periodicity::Periodic, true,
             {5e-5, 0.0035}, {1.3e-3, 0.015}, {0.1, 0.0055}},
            {"antiperiodic", argv[3], samples, Periodicity::Antiperiodic, true,
             {5e-5, 0.004}, {1.3e-3, 0.02}, {0.1, 0.0065}}
        };
        for (const auto &fixture : fixtures)
            compareFixture(fixture);
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
