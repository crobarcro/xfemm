#include "AnalysisSession.h"
#include "FSolverAnalysisBackend.h"
#include "NativeEquivalence.h"
#include "NativeMatrixFieldEquivalence.h"

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

// G3/G5: the independent native and materialized session pipelines must
// reproduce the same physical field, regardless of their solver node ordering.
bool testRingEquivalence()
{
    auto referenceSolver = std::make_shared<FSolverAnalysisBackend>();
    auto nativeSolver = std::make_shared<FSolverAnalysisBackend>();
    AnalysisSession reference(ModelDefinition(makeProblem()), referenceSolver);
    AnalysisSession native(ModelDefinition(makeProblem()), nativeSolver);
    auto instanced = makeRing(6);
    reference.setInstancedMesh(instanced);
    native.setInstancedMesh(std::move(instanced));
    native.setNativeInstanced(true);

    xfemm_test::EquivalenceOptions options;
    options.absoluteTolerance = options.relativeTolerance = 1e-7;
    xfemm_test::MatrixFieldOptions diagnostics;
    diagnostics.nodal = options;
    diagnostics.maximumRelativeResidual = 1e-4;
    const auto comparison = xfemm_test::compareMatrixFieldSessions(
        reference, *referenceSolver, native, *nativeSolver, diagnostics);
    if (!comparison.passed)
        std::cerr << "linear ring G5.2: " << comparison.message << '\n';
    REQUIRE(comparison.passed);
    REQUIRE(referenceSolver->solveCount() == 1);
    REQUIRE(nativeSolver->nativeSolveCount() == 1);
    REQUIRE(referenceSolver->topologyImportCount() == 1);
    return true;
}

// G3/G5: the native Newton iteration reproduces the materialized solution.
bool testNonlinearEquivalence()
{
    auto nonlinearProblem = [] {
        auto problem = makeProblem();
        auto *material = dynamic_cast<CMMaterialProp *>(problem->blockproplist[0].get());
        material->mu_x = material->mu_y = 1000;
        material->Bdata = {0.0, 0.5, 1.0, 1.5, 2.0};
        material->Hdata = {0.0, 100.0, 500.0, 2000.0, 10000.0};
        material->BHpoints = 5;
        auto *circuit = dynamic_cast<CMCircuit *>(problem->circproplist[0].get());
        circuit->Amps = 50;
        return problem;
    };

    auto referenceSolver = std::make_shared<FSolverAnalysisBackend>();
    auto nativeSolver = std::make_shared<FSolverAnalysisBackend>();
    AnalysisSession reference(ModelDefinition(nonlinearProblem()), referenceSolver);
    AnalysisSession native(ModelDefinition(nonlinearProblem()), nativeSolver);
    auto instanced = makeRing(6);
    reference.setInstancedMesh(instanced);
    native.setInstancedMesh(std::move(instanced));
    native.setNativeInstanced(true);

    xfemm_test::EquivalenceOptions options;
    options.absoluteTolerance = options.relativeTolerance = 1e-4;
    xfemm_test::MatrixFieldOptions diagnostics;
    diagnostics.nodal = options;
    diagnostics.matrixRelativeTolerance = 1e-3;
    diagnostics.rhsRelativeTolerance = 1e-3;
    diagnostics.maximumRelativeResidual = 1e-2;
    diagnostics.fieldRelativeTolerance = 1e-4;
    diagnostics.hRelativeTolerance = 1e-3;
    const auto comparison = xfemm_test::compareMatrixFieldSessions(
        reference, *referenceSolver, native, *nativeSolver, diagnostics);
    if (!comparison.passed)
        std::cerr << "nonlinear ring G5.2: " << comparison.message << '\n';
    REQUIRE(comparison.passed);
    REQUIRE(referenceSolver->solveCount() == 1);
    REQUIRE(nativeSolver->nativeSolveCount() == 1);
    return true;
}

/** Two triangles with ordinary PBC, independent of an air-gap coupling.
 * This checks the constrained matrix and RHS for both sign conventions.
 */
InstancedMesh makePeriodicSquare(bool antiperiodic)
{
    SolverMesh mesh;
    mesh.nodes = {{0., 0., 2}, {1., 0., 0}, {1., 1., 0}, {0., 1., 0}};
    mesh.elements = {{{{0, 1, 2}}, 1}, {{{0, 2, 3}}, 1}};
    mesh.edges = {{0, 1, 0}, {1, 2, 0}, {2, 3, 0}, {3, 0, 0}, {0, 2, 0}};
    mesh.periodicConstraints.push_back({
        1, 3, antiperiodic ? SolverMesh::Periodicity::Antiperiodic
                          : SolverMesh::Periodicity::Periodic
    });
    MeshTemplate square;
    square.localMesh = std::move(mesh);
    InstancedMesh instanced;
    instanced.templates.push_back(std::move(square));
    instanced.instances.push_back({});
    return instanced;
}

bool testOrdinaryPeriodicEquivalence(bool antiperiodic)
{
    auto referenceSolver = std::make_shared<FSolverAnalysisBackend>();
    auto nativeSolver = std::make_shared<FSolverAnalysisBackend>();
    AnalysisSession reference(ModelDefinition(makeProblem()), referenceSolver);
    AnalysisSession native(ModelDefinition(makeProblem()), nativeSolver);
    const auto instanced = makePeriodicSquare(antiperiodic);
    reference.setInstancedMesh(instanced);
    native.setInstancedMesh(instanced);
    native.setNativeInstanced(true);
    xfemm_test::MatrixFieldOptions options;
    options.nodal.absoluteTolerance = options.nodal.relativeTolerance = 1e-7;
    options.maximumRelativeResidual = 1e-4;
    const auto comparison = xfemm_test::compareMatrixFieldSessions(
        reference, *referenceSolver, native, *nativeSolver, options);
    if (!comparison.passed)
        std::cerr << (antiperiodic ? "antiperiodic" : "periodic")
                  << " PBC G5.2: " << comparison.message << '\n';
    REQUIRE(comparison.passed);
    REQUIRE(reference.mesh()->periodicConstraints.size() == 1);
    REQUIRE(native.logicalView()->periodicConstraints().size() == 1);
    REQUIRE(comparison.matrix.numericNonzeros > 0);
    return true;
}

bool testPeriodicEquivalence()
{
    return testOrdinaryPeriodicEquivalence(false);
}

bool testAntiperiodicEquivalence()
{
    return testOrdinaryPeriodicEquivalence(true);
}

// G3/E: a session-native solve must not materialize, and physics changes must
// reuse topology-derived ordering until the instance layout changes.
bool testSessionNativeMode()
{
    auto referenceSolver = std::make_shared<FSolverAnalysisBackend>();
    auto nativeSolver = std::make_shared<FSolverAnalysisBackend>();
    AnalysisSession reference(ModelDefinition(makeProblem()), referenceSolver);
    AnalysisSession native(ModelDefinition(makeProblem()), nativeSolver);
    auto instanced = makeRing(6);
    reference.setInstancedMesh(instanced);
    native.setInstancedMesh(std::move(instanced));
    native.setNativeInstanced(true);

    const auto comparison = xfemm_test::compareSessionSolutions(reference, native);
    if (!comparison.passed)
        std::cerr << "native-session: " << comparison.message << '\n';
    REQUIRE(comparison.passed);
    REQUIRE(referenceSolver->topologyImportCount() == 1);
    REQUIRE(nativeSolver->nativeSolveCount() == 1);
    REQUIRE(native.materializationCount() == 0);
    REQUIRE(native.viewGenerationCount() == 1);
    REQUIRE(nativeSolver->nativeOrderingBuildCount() == 1);

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
        {"periodic-matrix-equivalence", testPeriodicEquivalence},
        {"antiperiodic-matrix-equivalence", testAntiperiodicEquivalence},
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
