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
using Mesh = femm::mesh::SolverMesh;

struct Tolerance { double absolute; double relative; };
struct Observation { double angle; double torque; double energy; };
struct Run {
    std::string backend;
    std::unique_ptr<femm::AnalysisSession> session;
    std::shared_ptr<femm::FSolverAnalysisBackend> solver;
    std::shared_ptr<const Mesh> mesh;
    std::vector<Observation> observations;
};

[[noreturn]] void fail(const std::string &backend, const std::string &gap,
                       const std::string &invariant)
{
    throw std::runtime_error(backend + "/AGE '" + gap + "': " + invariant);
}

std::unique_ptr<femm::FemmProblem> readProblem(const std::string &path)
{
    auto problem = std::make_unique<femm::FemmProblem>(femm::FileType::MagneticsFile);
    auto view = std::shared_ptr<femm::FemmProblem>(problem.get(), [](femm::FemmProblem *) {});
    femm::MagneticsReader reader(view, std::cerr);
    if (reader.parse(path) != femm::F_FILE_OK)
        throw std::runtime_error("could not parse " + path);
    return problem;
}

std::shared_ptr<fmesher::MesherBackend> mesher(const std::string &name)
{
    if (name == "Triangle") return std::make_shared<fmesher::TriangleMesherBackend>();
    if (name == "Tangle") return std::make_shared<fmesher::TangleMesherBackend>();
    throw std::logic_error("unknown backend");
}

const Mesh::AirGap &onlyGap(const Mesh &mesh, const std::string &backend)
{
    if (mesh.airGaps.size() != 1)
        fail(backend, "AGE", "expected exactly one AGE, got " +
             std::to_string(mesh.airGaps.size()));
    if (mesh.airGaps.front().boundaryName != "AGE")
        fail(backend, mesh.airGaps.front().boundaryName, "unexpected boundary name");
    return mesh.airGaps.front();
}

void validateRing(const Mesh &mesh, const Mesh::AirGap &gap,
                  const std::vector<Mesh::AirGapRingPoint> &ring,
                  double radius, const std::string &backend, const char *side)
{
    if (ring.empty() || ring.size() < gap.totalArcElements)
        fail(backend, gap.boundaryName, std::string(side) + " ring has incompatible cardinality");
    bool positive = false, negative = false;
    double previous = -1e100;
    for (std::size_t i = 0; i < ring.size(); ++i) {
        const auto &point = ring[i];
        if (point.node >= mesh.nodes.size())
            fail(backend, gap.boundaryName, std::string(side) + " ring node " +
                 std::to_string(i) + " is invalid");
        if (!std::isfinite(point.elementPosition) || !std::isfinite(point.weight))
            fail(backend, gap.boundaryName, std::string(side) + " ring value is non-finite");
        if (point.elementPosition + 1e-10 < previous)
            fail(backend, gap.boundaryName, std::string(side) + " ring positions are unordered");
        previous = point.elementPosition;
        if (std::abs(std::abs(point.weight) - 1.0) > 1e-12)
            fail(backend, gap.boundaryName, std::string(side) + " ring weight is not a unit sign");
        positive = positive || point.weight > 0; negative = negative || point.weight < 0;
        const auto &node = mesh.nodes[point.node];
        const double actual = std::hypot(node.x - gap.centerX, node.y - gap.centerY);
        if (std::abs(actual - radius) > 2e-10)
            fail(backend, gap.boundaryName, std::string(side) + " ring node is off its circle");
    }
    if (gap.periodicity == Mesh::Periodicity::Antiperiodic && (!positive || !negative))
        fail(backend, gap.boundaryName, std::string(side) +
             " antiperiodic ring does not contain both signs");
    if (gap.periodicity == Mesh::Periodicity::Periodic && negative)
        fail(backend, gap.boundaryName, std::string(side) +
             " periodic ring contains a negative sign");
}

void validateGap(const Mesh &mesh, const std::string &backend)
{
    const auto validation = femm::mesh::validateSolverMesh(mesh);
    if (!validation.valid()) fail(backend, "AGE", "validateSolverMesh rejected mesh");
    if (mesh.nodes.empty() || mesh.elements.empty() || mesh.edges.empty())
        fail(backend, "AGE", "empty base mesh");
    if (mesh.periodicConstraints.empty())
        fail(backend, "AGE", "expected ordinary periodic constraints");
    for (const auto &constraint : mesh.periodicConstraints)
        if (constraint.periodicity != Mesh::Periodicity::Periodic)
            fail(backend, "AGE", "ordinary constraint has wrong periodicity");
    const auto &gap = onlyGap(mesh, backend);
    if (gap.periodicity != Mesh::Periodicity::Periodic)
        fail(backend, gap.boundaryName, "expected periodic AGE metadata");
    validateRing(mesh, gap, gap.innerRing, gap.innerRadius, backend, "inner");
    validateRing(mesh, gap, gap.outerRing, gap.outerRadius, backend, "outer");

    // The AGE replaces the ordinary finite-element annulus. Check centroids
    // comfortably away from either ring so boundary roundoff is irrelevant.
    const double annulusTolerance = 2e-10;
    for (std::size_t i = 0; i < mesh.elements.size(); ++i) {
        const auto &element = mesh.elements[i];
        double x = 0.0, y = 0.0;
        for (const auto node : element.nodes) {
            x += mesh.nodes[node].x;
            y += mesh.nodes[node].y;
        }
        const double radius = std::hypot(x / 3.0 - gap.centerX, y / 3.0 - gap.centerY);
        if (radius > gap.innerRadius + annulusTolerance &&
            radius < gap.outerRadius - annulusTolerance)
            fail(backend, gap.boundaryName, "ordinary element " + std::to_string(i) +
                 " occupies the AGE annulus");
    }
    if (gap.innerRing.size() != gap.outerRing.size())
        fail(backend, gap.boundaryName, "inner and outer full-ring cardinalities differ");
    if (gap.quadraturePoints.size() != gap.totalArcElements + 1)
        fail(backend, gap.boundaryName, "quadrature cardinality is not arc elements + 1");
    const auto ringContains = [](const std::vector<Mesh::AirGapRingPoint> &ring,
                                 femm::mesh::MeshIndex node, double weight) {
        return std::any_of(ring.begin(), ring.end(), [node, weight](const auto &point) {
            return point.node == node && std::abs(point.weight - weight) <= 1e-12;
        });
    };
    for (std::size_t q = 0; q < gap.quadraturePoints.size(); ++q) {
        const auto &point = gap.quadraturePoints[q];
        for (std::size_t i = 0; i < 4; ++i) {
            if (point.nodes[i] >= mesh.nodes.size() || !std::isfinite(point.weights[i]))
                fail(backend, gap.boundaryName, "invalid quadrature entry " + std::to_string(q));
            // FEMM AGE quadrature stores the periodic image sign, rather than
            // barycentric coefficients. Each side's adjacent values are unit signs.
            if (std::abs(std::abs(point.weights[i]) - 1.0) > 1e-12)
                fail(backend, gap.boundaryName, "quadrature weight is not a unit sign");
        }
        if (point.weights[0] != point.weights[2] || point.weights[1] != point.weights[3])
            fail(backend, gap.boundaryName, "inner/outer quadrature signs disagree");
        for (std::size_t i = 0; i < 4; ++i) {
            const bool inner = i < 2;
            const auto &ring = inner ? gap.innerRing : gap.outerRing;
            const double expectedRadius = inner ? gap.innerRadius : gap.outerRadius;
            if (!ringContains(ring, point.nodes[i], point.weights[i]))
                fail(backend, gap.boundaryName, "quadrature node/sign is absent from its " +
                     std::string(inner ? "inner" : "outer") + " ring");
            const auto &node = mesh.nodes[point.nodes[i]];
            const double radius = std::hypot(node.x-gap.centerX, node.y-gap.centerY);
            if (std::abs(radius-expectedRadius) > 2e-10)
                fail(backend, gap.boundaryName, "quadrature node is off its physical ring");
        }
    }
}

bool sameTopology(const Mesh::AirGap &a, const Mesh::AirGap &b)
{
    const auto sameRing = [](const std::vector<Mesh::AirGapRingPoint> &x,
                             const std::vector<Mesh::AirGapRingPoint> &y) {
        if (x.size() != y.size()) return false;
        for (std::size_t i = 0; i < x.size(); ++i)
            if (x[i].node != y[i].node || x[i].elementPosition != y[i].elementPosition ||
                x[i].weight != y[i].weight) return false;
        return true;
    };
    return sameRing(a.innerRing, b.innerRing) && sameRing(a.outerRing, b.outerRing) &&
           a.nodeIndices == b.nodeIndices;
}

Run sweep(const std::string &path, const std::string &backend)
{
    auto solver = std::make_shared<femm::FSolverAnalysisBackend>();
    auto session = std::make_unique<femm::AnalysisSession>(
        femm::ModelDefinition(readProblem(path)), mesher(backend), solver);
    const auto mesh = session->ensureMesh();
    validateGap(*mesh, backend);
    const auto canonicalGap = onlyGap(*mesh, backend);
    const auto topologyIdentity = session->meshTopologyIdentity();
    // 3.75 degrees is exactly one element for this 96-element full ring; it
    // directly exercises the zero/full-ring endpoint convention.
    const std::vector<double> angles{0.0, 3.75, 10.0, 30.0, 70.0};
    std::vector<Observation> observations;
    for (double angle : angles) {
        if (angle != 0.0)
            session->setAirGapAngle(session->model().airGap("AGE"), angle, 0.0);
        session->solve();
        FPProc postprocessor;
        if (!postprocessor.OpenDocument(session->model().problem(), solver->solvedSolver(),
                                        solver->solvedSystem()))
            fail(backend, "AGE", "in-memory post-processing failed");
        for (auto &label : postprocessor.blocklist) label.IsSelected = true;
        double torque = 0.0;
        if (postprocessor.gapDCTorqueIntegral("AGE", torque) != FPProcError::NoError)
            fail(backend, "AGE", "gapDCTorqueIntegral failed");
        const double energy = postprocessor.BlockIntegral(2).re;
        if (!std::isfinite(torque) || !std::isfinite(energy) || energy <= 1e-12)
            fail(backend, "AGE", "non-finite or trivial physical result");
        observations.push_back({angle, torque, energy});
    }
    if (session->meshGenerationCount() != 1 || solver->topologyImportCount() != 1 ||
        solver->orderingCount() != 1 || solver->solveCount() != angles.size() ||
        solver->couplingRegenerationCount() != angles.size() - 1 ||
        session->meshTopologyIdentity() != topologyIdentity)
        fail(backend, "AGE", "angle sweep rebuilt mesh/topology or has wrong evaluation counts");
    if (!sameTopology(canonicalGap, onlyGap(*session->mesh(), backend)))
        fail(backend, "AGE", "canonical reusable ring topology was mutated");
    const auto torqueRange = std::minmax_element(observations.begin(), observations.end(),
        [](const Observation &a, const Observation &b) { return a.torque < b.torque; });
    if (torqueRange.second->torque - torqueRange.first->torque < 0.25)
        fail(backend, "AGE", "torque did not materially respond to AGE movement");
    std::cout << backend << " AGE: ring-size=" << canonicalGap.innerRing.size()
              << " arc-elements=" << canonicalGap.totalArcElements
              << " radii=(" << canonicalGap.innerRadius << ',' << canonicalGap.outerRadius
              << ") center=(" << canonicalGap.centerX << ',' << canonicalGap.centerY
              << ") coupling-regenerations=" << solver->couplingRegenerationCount() << '\n';
    return {backend, std::move(session), std::move(solver), mesh, std::move(observations)};
}

void compareAnalyticalTorque(const Run &run)
{
    // This benchmark is normalized so its analytical torque is sin(angle).
    // Retain the established Lua benchmark's 0.02 absolute allowance while
    // requiring finite values at every session-positioned angle.
    constexpr double tolerance = 0.02;
    const double pi = std::acos(-1.0);
    for (const auto &observation : run.observations) {
        const double expected = std::sin(observation.angle * pi / 180.0);
        const double difference = std::abs(observation.torque - expected);
        std::cout << run.backend << " AGE angle=" << observation.angle
                  << " analytical torque: value=" << observation.torque
                  << " expected=" << expected << " abs-difference=" << difference
                  << " tolerance=" << tolerance << '\n';
        if (!std::isfinite(observation.torque) || !std::isfinite(difference) ||
            difference > tolerance)
            throw std::runtime_error(run.backend + " analytical AGE torque comparison failed");
    }
}

void compareValue(const char *quantity, double angle, double triangle, double tangle,
                  Tolerance tolerance)
{
    const double difference = std::abs(triangle - tangle);
    const double permitted = tolerance.absolute + tolerance.relative *
        std::max(std::abs(triangle), std::abs(tangle));
    std::cout << "AGE angle=" << angle << ' ' << quantity << ": Triangle=" << triangle
              << " Tangle=" << tangle << " abs-difference=" << difference
              << " tolerance=" << permitted << '\n';
    if (!std::isfinite(triangle) || !std::isfinite(tangle) || !std::isfinite(difference) ||
        !std::isfinite(permitted) || difference > permitted)
        throw std::runtime_error(std::string("AGE angle ") + std::to_string(angle) +
                                 " " + quantity + " differential failed");
}

void compareStructure(const Run &triangle, const Run &tangle)
{
    const auto &a = onlyGap(*triangle.mesh, triangle.backend);
    const auto &b = onlyGap(*tangle.mesh, tangle.backend);
    if (a.periodicity != b.periodicity || a.totalArcElements != b.totalArcElements ||
        a.innerRing.size() != b.innerRing.size() || a.outerRing.size() != b.outerRing.size() ||
        a.quadraturePoints.size() != b.quadraturePoints.size())
        throw std::runtime_error("AGE discrete Triangle/Tangle topology differs");
    const auto close = [](double x, double y) { return std::abs(x-y) <= 2e-12; };
    if (!close(a.innerRadius,b.innerRadius) || !close(a.outerRadius,b.outerRadius) ||
        !close(a.centerX,b.centerX) || !close(a.centerY,b.centerY) ||
        !close(a.totalArcLengthDegrees,b.totalArcLengthDegrees) ||
        !close(a.innerAngleDegrees,b.innerAngleDegrees) ||
        !close(a.outerAngleDegrees,b.outerAngleDegrees))
        throw std::runtime_error("AGE geometric Triangle/Tangle topology differs");
}

void periodicWithoutAge(const std::string &path, const std::string &backend)
{
    auto solver = std::make_shared<femm::FSolverAnalysisBackend>();
    femm::AnalysisSession session(femm::ModelDefinition(readProblem(path)), mesher(backend), solver);
    const auto mesh = session.ensureMesh();
    if (!femm::mesh::validateSolverMesh(*mesh).valid() || mesh->periodicConstraints.empty() ||
        !mesh->airGaps.empty())
        throw std::runtime_error(backend + " periodic-without-AGE semantics failed");
    session.solve();
    if (solver->topologyImportCount() != 1 || solver->solveCount() != 1 ||
        solver->couplingRegenerationCount() != 0)
        throw std::runtime_error(backend + " periodic-without-AGE solver path failed");
}
} // namespace

int main(int argc, char **argv)
{
    if (argc != 3) {
        std::cerr << "expected AGE and periodic-without-AGE FEM paths\n";
        return 2;
    }
    try {
        periodicWithoutAge(argv[2], "Triangle");
        periodicWithoutAge(argv[2], "Tangle");
        auto triangle = sweep(argv[1], "Triangle");
        auto tangle = sweep(argv[1], "Tangle");
        compareStructure(triangle, tangle);
        compareAnalyticalTorque(triangle);
        compareAnalyticalTorque(tangle);
        for (std::size_t i = 0; i < triangle.observations.size(); ++i) {
            const auto &a = triangle.observations[i];
            const auto &b = tangle.observations[i];
            compareValue("torque", a.angle, a.torque, b.torque, {7e-6, 1e-4});
            compareValue("energy", a.angle, a.energy, b.energy, {5e-5, 3e-5});
        }
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
