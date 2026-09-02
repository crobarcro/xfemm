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
    for (const auto &constraint : mesh->periodicConstraints)
        if (constraint.periodicity != fixture.periodicity)
            throw std::runtime_error(fixture.name + "/" + backendName + " has wrong PBC type");

    session->solve();
    auto postprocessor = std::make_unique<FPProc>();
    if (!postprocessor->OpenDocument(session->model().problem(), solver->solvedSolver(),
                                     solver->solvedSystem()))
        throw std::runtime_error(fixture.name + "/" + backendName + " post-processing failed");
    for (auto &label : postprocessor->blocklist)
        label.IsSelected = true;
    return {std::move(session), std::move(solver), std::move(postprocessor)};
}

void compare(const Fixture &fixture, const Point &point, const char *quantity,
             double triangle, double tangle, Tolerance tolerance)
{
    const double difference = std::abs(triangle - tangle);
    const double permitted = tolerance.absolute + tolerance.relative *
        std::max(std::abs(triangle), std::abs(tangle));
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
    for (const auto &point : fixture.samples) {
        CMPointVals tv, gv;
        if (!triangle.postprocessor->GetPointValues(point.x, point.y, tv) ||
            !tangle.postprocessor->GetPointValues(point.x, point.y, gv))
            throw std::runtime_error(fixture.name + " sample is outside a generated mesh");
        compare(fixture, point, "B1", tv.B1.re, gv.B1.re, fixture.fieldTolerance);
        compare(fixture, point, "B2", tv.B2.re, gv.B2.re, fixture.fieldTolerance);
    }
    const Point global{0, 0};
    compare(fixture, global, "energy", triangle.postprocessor->BlockIntegral(2).re,
            tangle.postprocessor->BlockIntegral(2).re, fixture.integralTolerance);
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
        // The 0.08 m target element size gives several hundred linear elements.
        // Point fields tolerate local element-gradient differences; domain integrals
        // are appreciably less mesh-sensitive and therefore use tighter limits.
        const std::vector<Point> samples{{-0.55, -0.35}, {-0.25, 0.30},
                                         {0.20, -0.20}, {0.55, 0.40}};
        const std::vector<Fixture> fixtures{
            {"non-periodic", argv[1], samples, Periodicity::Periodic, false,
             {1e-3, 0.03}, {1.0, 0.002}},
            {"periodic", argv[2], samples, Periodicity::Periodic, true,
             {2e-3, 0.02}, {1.0, 0.003}},
            {"antiperiodic", argv[3], samples, Periodicity::Antiperiodic, true,
             {2e-3, 0.025}, {1.0, 0.005}}
        };
        for (const auto &fixture : fixtures)
            compareFixture(fixture);
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
