#include "NativeCircuitEnergyEquivalence.h"

#include "AnalysisSession.h"
#include "FSolverAnalysisBackend.h"
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

bool require(bool ok, const std::string &label, const std::string &detail = {})
{
    if (!ok)
        std::cerr << "FAIL G5.3 " << label << ": " << detail << '\n';
    return ok;
}

std::unique_ptr<FemmProblem> makeProblem(bool nonlinear, bool opposing)
{
    auto problem = std::make_unique<FemmProblem>(FileType::MagneticsFile);
    problem->Depth = 0.15;
    problem->Precision = 1e-9;

    auto fixed = std::make_unique<CMPointProp>();
    fixed->PointName = "fixed";
    fixed->A = 0;
    problem->nodeproplist.push_back(std::move(fixed));

    auto steel = std::make_unique<CMMaterialProp>();
    steel->BlockName = "steel";
    steel->mu_x = steel->mu_y = 1000;
    steel->LamType = 0;
    steel->LamFill = 1;
    if (nonlinear) {
        steel->Bdata = {0., 0.5, 1., 1.5, 2.};
        steel->Hdata = {0., 100., 500., 2000., 10000.};
        steel->BHpoints = 5;
    }
    problem->blockproplist.push_back(std::move(steel));

    auto phase = std::make_unique<CMCircuit>();
    phase->CircName = "phase";
    phase->CircType = 1; // Series; expanded internally to per-label subcircuits.
    phase->Amps = nonlinear ? 50. : 3.;
    problem->circproplist.push_back(std::move(phase));

    if (opposing) {
        auto sense = std::make_unique<CMCircuit>();
        sense->CircName = "sense";
        sense->CircType = 1;
        sense->Amps = 0.; // Passive coil: flux linkage must still be measured.
        problem->circproplist.push_back(std::move(sense));
    }

    auto label = std::make_unique<CMBlockLabel>();
    label->BlockType = 0;
    label->BlockTypeName = "steel";
    label->InCircuit = 0;
    label->Turns = 2;
    problem->labellist.push_back(std::move(label));
    return problem;
}

InstancedMesh makeRing(bool opposing)
{
    const std::size_t count = 6;
    const double step = 360. / static_cast<double>(count);
    const double angle = step * pi / 180.;
    SolverMesh mesh;
    mesh.nodes = {{1., 0., 2}, {2., 0., 0},
                  {std::cos(angle), std::sin(angle), 0},
                  {2.*std::cos(angle), 2.*std::sin(angle), 0}};
    mesh.elements = {{{{0, 1, 3}}, 1}, {{{0, 3, 2}}, 1}};
    mesh.edges = {{0, 1, 0}, {1, 3, 0}, {2, 3, 0}, {0, 2, 0}, {0, 3, 0}};

    MeshTemplate prototype;
    prototype.localMesh = std::move(mesh);
    prototype.seams = {{"left", {0, 1}}, {"right", {2, 3}}};
    InstancedMesh result;
    result.templates.push_back(std::move(prototype));
    for (std::size_t i = 0; i < count; ++i) {
        MeshInstance instance;
        instance.templateIndex = 0;
        instance.transform.rotationDegrees = i * step;
        SeamConnection connection;
        connection.seam = 1;
        connection.otherInstance = (i + 1) % count;
        connection.otherSeam = 0;
        instance.seamConnections.push_back(connection);
        if (opposing) {
            InstanceRegionOverride override;
            override.sourceBlockLabel = 0;
            if (i == 2 || i == 4) override.circuit = 1;
            if (i == 3 || i == 4) override.currentScale = -1.;
            if (i == 5) override.currentScale = 2.;
            instance.regionOverrides.push_back(override);
        }
        result.instances.push_back(std::move(instance));
    }
    return result;
}

struct Fixture {
    std::shared_ptr<FSolverAnalysisBackend> referenceSolver =
        std::make_shared<FSolverAnalysisBackend>();
    std::shared_ptr<FSolverAnalysisBackend> nativeSolver =
        std::make_shared<FSolverAnalysisBackend>();
    AnalysisSession reference;
    AnalysisSession native;

    Fixture(bool nonlinear, bool opposing)
        : reference(ModelDefinition(makeProblem(nonlinear, opposing)), referenceSolver),
          native(ModelDefinition(makeProblem(nonlinear, opposing)), nativeSolver)
    {
        auto ring = makeRing(opposing);
        reference.setInstancedMesh(ring);
        native.setInstancedMesh(std::move(ring));
        native.setNativeInstanced(true);
    }
};

bool testLinearSeries()
{
    Fixture setup(false, false);
    xfemm_test::CircuitEnergyOptions opts;
    opts.field.maximumRelativeResidual = 1e-4;
    auto result = xfemm_test::compareCircuitEnergySessions(
        setup.reference, *setup.referenceSolver,
        setup.native, *setup.nativeSolver, opts);
    if (!require(result.passed, "linear series", result.message)) return false;
    if (!require(result.reference.circuits.size() == 1, "one series circuit")) return false;
    if (!require(result.reference.circuits[0].fluxLinkage != 0,
                 "series circuit does not silently report placeholder zero")) return false;
    if (!require(result.reference.energy > 0, "positive magnetic energy")) return false;
    const double gap = std::abs(result.reference.energy - result.reference.coenergy);
    if (!require(gap <= 1e-6 * std::max(result.reference.energy, 1e-9),
                 "linear energy equals coenergy", std::to_string(gap))) return false;
    if (!require(setup.native.materializationCount() == 0,
                 "native solver has no materialization")) return false;
    std::cout << "G5.3 linear " << result.message << '\n';
    return true;
}

bool testOpposingAndPassive()
{
    Fixture setup(false, true);
    xfemm_test::CircuitEnergyOptions opts;
    opts.field.maximumRelativeResidual = 1e-4;
    auto result = xfemm_test::compareCircuitEnergySessions(
        setup.reference, *setup.referenceSolver,
        setup.native, *setup.nativeSolver, opts);
    if (!require(result.passed, "opposing and passive circuits", result.message))
        return false;
    if (!require(result.reference.circuits.size() == 2,
                 "two original user-visible circuits")) return false;
    if (!require(result.reference.circuits[0].current == 3. &&
                 result.reference.circuits[1].current == 0.,
                 "prescribed versus passive currents")) return false;
    if (!require(std::isfinite(result.reference.circuits[1].fluxLinkage),
                 "passive flux linkage is finite")) return false;
    if (!require(std::abs(result.reference.circuits[1].fluxLinkage) > 1e-12,
                 "passive coil has induced linkage from the driven phase"))
        return false;

    // A different field must change measured energy; comparing two calls
    // with reference-field values would be a vacuous test.
    const auto &rhs = setup.referenceSolver->solvedSystem().rhs();
    std::vector<double> altered(rhs.begin(), rhs.end());
    for (double &a : altered) a *= 1.1;
    const auto changed = xfemm_test::measureCircuitEnergy(
        setup.reference.model().problem(), setup.referenceSolver->solvedSolver(),
        setup.referenceSolver->solvedSolver(), altered);
    if (!require(std::abs(changed.energy - result.reference.energy) >
                     1e-5 * result.reference.energy,
                 "postprocessor measures supplied field, not cached reference"))
        return false;

    std::cout << "G5.3 opposed " << result.message << '\n';
    return true;
}

bool testNonlinearEnergy()
{
    Fixture setup(true, false);
    xfemm_test::CircuitEnergyOptions opts;
    opts.field.nodal.absoluteTolerance = opts.field.nodal.relativeTolerance = 1e-4;
    opts.field.matrixRelativeTolerance = 1e-3;
    opts.field.rhsRelativeTolerance = 1e-3;
    opts.field.maximumRelativeResidual = 1e-2;
    opts.field.fieldRelativeTolerance = 1e-4;
    opts.field.hRelativeTolerance = 1e-3;
    opts.fluxRelativeTolerance = 1e-3;
    opts.energyRelativeTolerance = 1e-3;
    auto result = xfemm_test::compareCircuitEnergySessions(
        setup.reference, *setup.referenceSolver,
        setup.native, *setup.nativeSolver, opts);
    if (!require(result.passed, "nonlinear energy/coenergy", result.message))
        return false;
    if (!require(result.reference.energy > 0 && result.reference.coenergy > 0,
                 "nonlinear magnetic energy and coenergy are positive")) return false;
    std::cout << "G5.3 nonlinear " << result.message << '\n';
    return true;
}

bool testInputGuard()
{
    Fixture setup(false, false);
    // First prepare a real materialized solve; a vector of the wrong dimension
    // must never silently read the solver's previous RHS.
    setup.reference.solve();
    try {
        (void)xfemm_test::measureCircuitEnergy(
            setup.reference.model().problem(), setup.referenceSolver->solvedSolver(),
            setup.referenceSolver->solvedSolver(), {1.});
        return require(false, "reject incorrect nodal field dimensions");
    } catch (const std::invalid_argument &) {
        return true;
    }
}

} // namespace

int main()
{
    bool success = true;
    success = testLinearSeries() && success;
    success = testOpposingAndPassive() && success;
    success = testNonlinearEnergy() && success;
    success = testInputGuard() && success;
    if (success)
        std::cout << "G5.3 circuit, energy and coenergy equivalence tests passed\n";
    return success ? 0 : 1;
}
