#include "NativeForceTorqueEquivalence.h"

#include "AnalysisSession.h"
#include "FSolverAnalysisBackend.h"
#include "FemmReader.h"
#include "TangleMesherBackend.h"
#include "CBlockLabel.h"
#include "CCircuit.h"
#include "CMaterialProp.h"
#include "CPointProp.h"
#include "mesh/InstancedMesh.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using namespace femm;
using namespace femm::mesh;

const double pi = std::acos(-1.0);

void check(bool condition, const std::string &label)
{
    if (!condition) throw std::runtime_error(label);
}

std::unique_ptr<FemmProblem> coilProblem()
{
    auto p = std::make_unique<FemmProblem>(FileType::MagneticsFile);
    p->Depth = 0.15;
    p->Precision = 1e-9;

    const double b = 0.15;
    for (double sign : {-1., 1.}) {
        auto fixed = std::make_unique<CMPointProp>();
        fixed->PointName = sign > 0 ? "A_top" : "A_bottom";
        fixed->A = b * sign * 1.25;
        p->nodeproplist.push_back(std::move(fixed));
    }
    auto coil = std::make_unique<CMMaterialProp>();
    coil->BlockName = "coil";
    coil->mu_x = coil->mu_y = 1.;
    p->blockproplist.push_back(std::move(coil));
    auto air = std::make_unique<CMMaterialProp>();
    air->BlockName = "air";
    air->mu_x = air->mu_y = 1.;
    p->blockproplist.push_back(std::move(air));

    auto circuit = std::make_unique<CMCircuit>();
    circuit->CircName = "drive";
    circuit->CircType = 1;
    circuit->Amps = 3.;
    p->circproplist.push_back(std::move(circuit));
    auto label = std::make_unique<CMBlockLabel>();
    label->BlockType = 0;
    label->BlockTypeName = "coil";
    label->InCircuit = 0;
    label->Turns = 4;
    p->labellist.push_back(std::move(label));
    auto background = std::make_unique<CMBlockLabel>();
    background->BlockType = 1;
    background->BlockTypeName = "air";
    background->InCircuit = -1;
    background->Turns = 0;
    p->labellist.push_back(std::move(background));
    return p;
}

// A magnetically driven central coil, an annular air buffer, and an outer
// Dirichlet rectangle imposing A=B_x*y. The two air rings isolate the coil
// from the boundary and allow a valid FEMM weighted-stress mask.
// The coil is centred, so self-field Lorentz force cancels by symmetry;
// the imposed uniform field creates an independent nonzero force prediction.
InstancedMesh biasedCoil()
{
    SolverMesh mesh;
    for (double half : {0.25, 0.75, 1.25}) {
        const int markerBottom = half == 1.25 ? 2 : 0;
        const int markerTop = half == 1.25 ? 3 : 0;
        mesh.nodes.push_back({-half, -half, markerBottom});
        mesh.nodes.push_back({ half, -half, markerBottom});
        mesh.nodes.push_back({ half,  half, markerTop});
        mesh.nodes.push_back({-half,  half, markerTop});
    }
    mesh.elements.push_back({{{0, 1, 2}}, 1});
    mesh.elements.push_back({{{0, 2, 3}}, 1});
    for (int layer = 0; layer < 2; ++layer) {
        const int a = 4 * layer;
        const int b = a + 4;
        for (int i = 0; i < 4; ++i) {
            const int j = (i + 1) % 4;
            mesh.elements.push_back({{{a+i, b+i, b+j}}, 2});
            mesh.elements.push_back({{{a+i, b+j, a+j}}, 2});
        }
    }
    for (int i = 0; i < 4; ++i)
        mesh.edges.push_back({static_cast<MeshIndex>(8+i),
                              static_cast<MeshIndex>(8+(i+1)%4), 0});

    MeshTemplate templateMesh;
    templateMesh.localMesh = std::move(mesh);
    InstancedMesh result;
    result.templates.push_back(std::move(templateMesh));
    result.instances.push_back({});
    return result;
}

struct Pair {
    std::shared_ptr<FSolverAnalysisBackend> refBackend =
        std::make_shared<FSolverAnalysisBackend>();
    std::shared_ptr<FSolverAnalysisBackend> nativeBackend =
        std::make_shared<FSolverAnalysisBackend>();
    AnalysisSession reference;
    AnalysisSession native;

    Pair(std::unique_ptr<FemmProblem> refProblem,
         std::unique_ptr<FemmProblem> nativeProblem, const InstancedMesh &instanced)
        : reference(ModelDefinition(std::move(refProblem)), refBackend),
          native(ModelDefinition(std::move(nativeProblem)), nativeBackend)
    {
        reference.setInstancedMesh(instanced);
        native.setInstancedMesh(instanced);
        native.setNativeInstanced(true);
    }
};

void biasedCoilForce()
{
    auto mesh = biasedCoil();
    Pair fixture(coilProblem(), coilProblem(), mesh);

    xfemm_test::ForceTorqueOptions opt;
    opt.selectedLabels = {0};
    opt.compareLorentz = true;
    opt.compareWeightedStress = true;
    opt.requireNonzeroLorentz = true;
    opt.requireNonzeroWeightedStress = true;
    opt.baseline.field.maximumRelativeResidual = 1e-3;
    auto result = xfemm_test::compareForceTorqueSessions(
        fixture.reference, *fixture.refBackend,
        fixture.native, *fixture.nativeBackend, opt);
    check(result.passed, "biased coil: " + result.message);
    check(fixture.native.materializationCount() == 0,
          "biased coil: native solver materialized a mesh");
    check(fixture.refBackend->solveCount() == 1 &&
          fixture.nativeBackend->nativeSolveCount() == 1,
          "biased coil: both paths must solve exactly once");
    // The externally imposed uniform Bx=0.15 T acts on a 3A*4turn
    // current through a 0.15m-depth planar domain:
    // F_y = Bx * I * N * Depth = 0.27 N. The self-field part integrates
    // to zero on the symmetric coil (allow finite-mesh numerical error).
    const double expectedFy = 0.15 * 3. * 4. * 0.15;
    check(std::abs(result.reference.lorentzFy - expectedFy) < 0.03,
          "biased coil: wrong analytical Lorentz force: " +
          std::to_string(result.reference.lorentzFy));
    check(std::abs(result.reference.lorentzFx) < 0.03,
          "biased coil: unexpected lateral Lorentz force");
    check(std::abs(result.reference.stressFy) > 1e-5 ||
          std::abs(result.reference.stressFx) > 1e-5,
          "biased coil: weighted stress reported a vacuous zero");
    std::cout << "G5.4 biased coil " << result.message
              << " expected_Lorentz_Fy_N=" << expectedFy << '\n';

    // A distinct test-side field must change the measured Lorentz load;
    // otherwise two zero/borrowed-result comparisons would pass vacuously.
    const auto &rhs = fixture.refBackend->solvedSystem().rhs();
    std::vector<double> changed(rhs.begin(), rhs.end());
    for (double &a : changed) a *= 1.1;
    const auto altered = xfemm_test::measureForceTorque(
        fixture.reference.model().problem(),
        fixture.refBackend->solvedSolver(), fixture.refBackend->solvedSolver(),
        changed, opt);
    check(std::abs(altered.lorentzFy - result.reference.lorentzFy) > 1e-3,
          "biased coil: Lorentz integral ignores the supplied nodal field");

    opt.selectedLabels = {5000};
    bool rejected = false;
    try {
        (void)xfemm_test::measureForceTorque(
            fixture.reference.model().problem(),
            fixture.refBackend->solvedSolver(), fixture.refBackend->solvedSolver(),
            changed, opt);
    } catch (const std::invalid_argument &) {
        rejected = true;
    }
    check(rejected, "invalid force selection was not rejected");
}

std::unique_ptr<FemmProblem> loadBenchmark(const std::string &path)
{
    auto problem = std::make_unique<FemmProblem>(FileType::MagneticsFile);
    auto borrowed = std::shared_ptr<FemmProblem>(problem.get(), [](FemmProblem *) {});
    MagneticsReader reader(borrowed, std::cerr);
    if (reader.parse(path) != F_FILE_OK)
        throw std::runtime_error("could not load G5.4 AGE benchmark: " + path);
    return problem;
}

// The established FEMM anti-periodic AGE benchmark's torque is sin(theta)
// N m to within 0.02 N m, independent of the native/materialized comparison.
// Use one template of the real benchmark mesh here to avoid redraw effects;
// multi-tile machine coverage is reserved for G5.5.
void ageTorque(const std::string &path)
{
    auto source = std::make_shared<FSolverAnalysisBackend>();
    AnalysisSession meshing(ModelDefinition(loadBenchmark(path)),
                            std::make_shared<fmesher::TangleMesherBackend>(), source);
    auto mesh = meshing.ensureMesh();
    check(mesh && !mesh->nodes.empty() && !mesh->airGaps.empty(),
          "G5.4 benchmark mesh lacks an AGE");
    MeshTemplate prototype;
    prototype.localMesh = *mesh;
    InstancedMesh instanced;
    instanced.templates.push_back(std::move(prototype));
    instanced.instances.push_back({});

    Pair fixture(loadBenchmark(path), loadBenchmark(path), instanced);
    xfemm_test::ForceTorqueOptions opt;
    opt.airGapBoundary = "AGE";
    opt.requireNonzeroAirGapTorque = true;
    opt.absoluteTorqueTolerance = 1e-5;
    opt.relativeTorqueTolerance = 2e-3;
    opt.baseline.requireNonzeroLinkage = false; // permanent-magnet model, no circuits
    opt.baseline.field.nodal.absoluteTolerance = 1e-5;
    opt.baseline.field.nodal.relativeTolerance = 1e-5;
    opt.baseline.field.maximumRelativeResidual = 2e-2;
    opt.baseline.field.residualAgreementTolerance = 1e-4;
    opt.baseline.field.matrixRelativeTolerance = 1e-5;
    opt.baseline.field.fieldRelativeTolerance = 1e-4;
    opt.baseline.field.hRelativeTolerance = 1e-4;

    for (double angle : {15., 30., 70.}) {
        fixture.reference.setAirGapAngle(fixture.reference.model().airGap("AGE"), angle, 0.);
        fixture.native.setAirGapAngle(fixture.native.model().airGap("AGE"), angle, 0.);
        const auto result = xfemm_test::compareForceTorqueSessions(
            fixture.reference, *fixture.refBackend,
            fixture.native, *fixture.nativeBackend, opt);
        check(result.passed, "AGE at " + std::to_string(angle) + ": " + result.message);
        const double analytic = std::sin(angle * pi / 180.);
        check(std::abs(result.reference.airGapTorque - analytic) < 0.02,
              "AGE analytical torque: angle=" + std::to_string(angle) +
              " expected=" + std::to_string(analytic) +
              " actual=" + std::to_string(result.reference.airGapTorque));
        check(fixture.native.materializationCount() == 0,
              "G5.4 AGE native path unexpectedly materialized");
        std::cout << "G5.4 AGE angle=" << angle << ' ' << result.message
                  << " analytical_torque_Nm=" << analytic << '\n';
    }
    check(fixture.nativeBackend->nativeSolveCount() == 3 &&
          fixture.refBackend->solveCount() == 3,
          "G5.4 AGE expected exactly three solves on each path");
    check(fixture.native.viewGenerationCount() == 1,
          "G5.4 AGE angle changes rebuilt the native mesh");
}

} // namespace

int main(int argc, char **argv)
{
    if (argc != 2) {
        std::cerr << "usage: native_force_torque_equivalence_test AGE_fixture.fem\n";
        return 2;
    }
    try {
        biasedCoilForce();
        ageTorque(argv[1]);
        std::cout << "G5.4 native force and torque equivalence passed\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "G5.4 FAILED: " << e.what() << '\n';
        return 1;
    }
}
