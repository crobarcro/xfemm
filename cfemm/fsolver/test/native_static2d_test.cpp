#include "AnalysisSession.h"
#include "FSolverAnalysisBackend.h"

#include "CBlockLabel.h"
#include "CCircuit.h"
#include "CMaterialProp.h"
#include "CPointProp.h"
#include "mesh/LogicalMeshView.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using namespace femm;
using namespace femm::mesh;

constexpr double Pi = 3.141592653589793238462643383279502884;

#define REQUIRE(condition)                                                          \
    do {                                                                            \
        if (!(condition)) {                                                         \
            std::cerr << __FILE__ << ':' << __LINE__ << ": " << #condition << '\n'; \
            return false;                                                           \
        }                                                                           \
    } while (false)

std::unique_ptr<FemmProblem> makeProblem()
{
    auto problem = std::make_unique<FemmProblem>(FileType::MagneticsFile);

    auto point = std::make_unique<CMPointProp>();
    point->PointName = "fixed";
    point->A = 0;
    point->J = 0;
    problem->nodeproplist.push_back(std::move(point));

    auto material = std::make_unique<CMMaterialProp>();
    material->BlockName = "steel";
    material->mu_x = material->mu_y = 1000;
    material->H_c = 0;
    material->Cduct = 0;
    material->LamType = 0;
    material->LamFill = 1;
    problem->blockproplist.push_back(std::move(material));

    auto circuit = std::make_unique<CMCircuit>();
    circuit->CircName = "phase";
    circuit->CircType = 1;
    circuit->Amps = 1;
    problem->circproplist.push_back(std::move(circuit));

    auto label = std::make_unique<CMBlockLabel>();
    label->BlockType = 0;
    label->BlockTypeName = "steel";
    label->InCircuit = 0;
    label->Turns = 1;
    problem->labellist.push_back(std::move(label));
    return problem;
}

MeshTemplate makeWedgeTemplate(double innerRadius, double outerRadius, double angleDegrees)
{
    const double radians = angleDegrees * Pi / 180.0;
    const double cosine = std::cos(radians);
    const double sine = std::sin(radians);
    SolverMesh mesh;
    // Node 0 carries point marker 2, which pins A = 0 and fixes the gauge.
    mesh.nodes = {
        {innerRadius, 0.0, 2},
        {outerRadius, 0.0, 0},
        {innerRadius * cosine, innerRadius * sine, 0},
        {outerRadius * cosine, outerRadius * sine, 0},
    };
    mesh.elements.push_back({{{0, 1, 3}}, 1});
    mesh.elements.push_back({{{0, 3, 2}}, 1});
    mesh.edges = {{0, 1, 0}, {1, 3, 0}, {2, 3, 0}, {0, 2, 0}, {0, 3, 0}};
    MeshTemplate meshTemplate;
    meshTemplate.localMesh = mesh;
    meshTemplate.seams = {{"left", {0, 1}}, {"right", {2, 3}}};
    return meshTemplate;
}

InstancedMesh makeRing(std::size_t count)
{
    const double angleStep = 360.0 / static_cast<double>(count);
    InstancedMesh instanced;
    instanced.templates.push_back(makeWedgeTemplate(1.0, 2.0, angleStep));
    for (std::size_t k = 0; k < count; ++k) {
        MeshInstance instance;
        instance.templateIndex = 0;
        instance.transform.rotationDegrees = static_cast<double>(k) * angleStep;
        SeamConnection connection;
        connection.seam = 1;
        connection.otherInstance = (k + 1) % count;
        connection.otherSeam = 0;
        instance.seamConnections.push_back(connection);
        instanced.instances.push_back(instance);
    }
    return instanced;
}

std::map<std::pair<long long, long long>, std::size_t> indexByCoordinate(const FSolver &solver)
{
    std::map<std::pair<long long, long long>, std::size_t> lookup;
    for (std::size_t i = 0; i < solver.meshnode.size(); ++i) {
        const long long x = std::llround(solver.meshnode[i].x * 1e6);
        const long long y = std::llround(solver.meshnode[i].y * 1e6);
        lookup[{x, y}] = i;
    }
    return lookup;
}

/** Compare the native field with the materialised field node by node. */
bool compareNativeToMaterialized(const std::shared_ptr<FSolverAnalysisBackend> &solver,
                                 AnalysisSession &session, double tolerance)
{
    const std::size_t materializedNodes = solver->solvedSolver().NumNodes;
    const auto materializedLookup = indexByCoordinate(solver->solvedSolver());
    const auto &materializedSystem = solver->solvedSystem();
    const std::vector<double> materializedA(materializedSystem.rhs().begin(),
                                            materializedSystem.rhs().end());

    std::vector<MaterializationDiagnostic> diagnostics;
    const LogicalMeshView view = LogicalMeshView::build(*session.instancedMesh(), diagnostics);
    REQUIRE(diagnostics.empty());
    REQUIRE(view.valid());
    REQUIRE(view.nodeCount() == materializedNodes);

    REQUIRE(solver->solveNative(view, session.instanceLabelBases(),
                                session.airGapPositioning()));

    const auto &nativeSystem = solver->nativeSystem();
    REQUIRE(nativeSystem.dimension() == static_cast<int>(view.nodeCount()));

    double largest = 0.0;
    for (std::size_t node = 0; node < view.nodeCount(); ++node) {
        double x = 0.0, y = 0.0;
        view.nodeCoordinates(node, x, y);
        const auto found = materializedLookup.find(
            {std::llround(x * 100.0 * 1e6), std::llround(y * 100.0 * 1e6)});
        REQUIRE(found != materializedLookup.end());
        const double nativeValue = nativeSystem.rhs()[node];
        const double reference = materializedA[found->second];
        largest = std::max(largest, std::abs(reference));
        REQUIRE(std::abs(nativeValue - reference) < tolerance * (1.0 + std::abs(reference)));
    }
    REQUIRE(largest > 1e-9);
    return true;
}

// G3/G5: the native assembly reproduces the materialised linear solver field.
bool testRingEquivalence()
{
    auto solver = std::make_shared<FSolverAnalysisBackend>();
    AnalysisSession session(ModelDefinition(makeProblem()), solver);
    session.setInstancedMesh(makeRing(6));
    session.solve();
    REQUIRE(solver->solveCount() == 1);
    REQUIRE(compareNativeToMaterialized(solver, session, 1e-7));
    return true;
}

// G3/G5: the native Newton iteration reproduces a nonlinear materialised solve.
bool testNonlinearEquivalence()
{
    auto problem = makeProblem();
    auto *material = dynamic_cast<CMMaterialProp *>(problem->blockproplist[0].get());
    material->mu_x = material->mu_y = 1000;
    material->Bdata = {0.0, 0.5, 1.0, 1.5, 2.0};
    material->Hdata = {0.0, 100.0, 500.0, 2000.0, 10000.0};
    material->BHpoints = 5;
    auto *circuit = dynamic_cast<CMCircuit *>(problem->circproplist[0].get());
    circuit->Amps = 50;

    auto solver = std::make_shared<FSolverAnalysisBackend>();
    AnalysisSession session(ModelDefinition(std::move(problem)), solver);
    session.setInstancedMesh(makeRing(6));
    session.solve();
    REQUIRE(solver->solveCount() == 1);
    REQUIRE(compareNativeToMaterialized(solver, session, 1e-4));
    return true;
}

// G3/E: the session-native mode skips materialisation and reuses the view.
bool testSessionNativeMode()
{
    auto materializedSolver = std::make_shared<FSolverAnalysisBackend>();
    AnalysisSession materialized(ModelDefinition(makeProblem()), materializedSolver);
    materialized.setInstancedMesh(makeRing(6));
    const TrialSolution reference = materialized.solve();
    REQUIRE(materializedSolver->topologyImportCount() == 1);
    REQUIRE(reference.real);

    std::map<std::pair<long long, long long>, std::size_t> lookup;
    for (std::size_t i = 0; i < reference.real->nodal.x.size(); ++i) {
        lookup[{std::llround(reference.real->nodal.x[i] * 1e6),
                std::llround(reference.real->nodal.y[i] * 1e6)}] = i;
    }

    auto nativeSolver = std::make_shared<FSolverAnalysisBackend>();
    AnalysisSession native(ModelDefinition(makeProblem()), nativeSolver);
    native.setInstancedMesh(makeRing(6));
    native.setNativeInstanced(true);
    const TrialSolution trial = native.solve();
    REQUIRE(nativeSolver->nativeSolveCount() == 1);
    REQUIRE(native.materializationCount() == 0);
    REQUIRE(native.viewGenerationCount() == 1);
    REQUIRE(trial.real);
    REQUIRE(trial.real->nodal.x.size() == reference.real->nodal.x.size());

    double largest = 0.0;
    for (std::size_t i = 0; i < trial.real->nodal.x.size(); ++i) {
        const auto found = lookup.find({std::llround(trial.real->nodal.x[i] * 1e6),
                                        std::llround(trial.real->nodal.y[i] * 1e6)});
        REQUIRE(found != lookup.end());
        const double value = trial.real->nodal.magneticVectorPotential[i];
        const double expected = reference.real->nodal.magneticVectorPotential[found->second];
        largest = std::max(largest, std::abs(expected));
        REQUIRE(std::abs(value - expected) < 1e-7 * (1.0 + std::abs(expected)));
    }
    REQUIRE(largest > 1e-9);
    REQUIRE(nativeSolver->nativeOrderingBuildCount() == 1);

    // Physics changes reuse both the view and its ordering; only a transform
    // rebuilds them. Translate the whole ring rigidly so welded seams stay
    // coincident.
    native.setCircuitCurrent(native.model().circuit("phase"), CComplex(2, 0));
    native.solve();
    REQUIRE(native.viewGenerationCount() == 1);
    REQUIRE(nativeSolver->nativeOrderingBuildCount() == 1);
    for (std::size_t i = 0; i < native.instancedMesh()->instances.size(); ++i) {
        RigidTransform2D moved = native.instancedMesh()->instances[i].transform;
        moved.translationXMetres += 0.01;
        native.setInstanceTransform(i, moved);
    }
    native.solve();
    REQUIRE(native.viewGenerationCount() == 2);
    REQUIRE(nativeSolver->nativeOrderingBuildCount() == 2);
    return true;
}

} // namespace

int main()
{
    const std::vector<std::pair<std::string, bool (*)()>> tests = {
        {"ring-equivalence", testRingEquivalence},
        {"nonlinear-equivalence", testNonlinearEquivalence},
        {"session-native-mode", testSessionNativeMode},
    };
    bool success = true;
    for (const auto &test : tests) {
        if (!test.second()) {
            std::cerr << "FAILED: " << test.first << '\n';
            success = false;
        }
    }
    if (success)
        std::cout << "native static2d tests passed\n";
    return success ? 0 : 1;
}
