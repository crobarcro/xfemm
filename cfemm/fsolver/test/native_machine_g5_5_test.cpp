/*
 * G5.5: checked-in RNFoundry machine, native vs materialized.
 *
 * Fast CI: one 60-degree sector (one rotor pole-pair and six stator slots),
 * at three relative rotor angles / two winding excitations.
 * Extended: --full retains the complete 6-rotor/36-stator model. Do not
 * assume a redraw model is pointwise comparable with an AGE motion model.
 * Fixture provenance: mfemm/testing/radial_machine/FIXTURE_PROVENANCE.md
 */
#include "AnalysisSession.h"
#include "FSolverAnalysisBackend.h"
#include "NativeForceTorqueEquivalence.h"
#include "TangleMesherBackend.h"
#include "TiledModel.h"
#include "TiledModelJson.h"
#include "TiledModelMesher.h"
#include "mesh/LogicalMeshView.h"

#include "CBoundaryProp.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using namespace femm;

void check(bool ok, const std::string &why)
{
    if (!ok) throw std::runtime_error(why);
}

std::string readFile(const std::string &path)
{
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot read pinned machine JSON: " + path);
    std::ostringstream text;
    text << input.rdbuf();
    return text.str();
}

void sectorise(tiled::TiledModel &model, unsigned divisor)
{
    if (divisor == 1) return;
    for (auto &tile : model.tiles) {
        check(tile.repeat.count % divisor == 0, "sector divisor does not divide tile instances");
        tile.repeat.count /= divisor;
        tile.repeat.closure = tiled::Closure::Periodic;
    }
    for (auto &override : model.overrides) {
        std::size_t count = 0;
        for (const auto &tile : model.tiles)
            if (tile.name == override.tile) count = tile.repeat.count;
        check(count != 0, "orphan machine region override");
        if (override.circuit.size() > count) override.circuit.resize(count);
        if (override.turnScale.size() > count) override.turnScale.resize(count);
        if (override.magDir.size() > count) override.magDir.resize(count);
    }
}

std::unique_ptr<FemmProblem> sessionProblem(const tiled::TiledModel &model)
{
    auto problem = tiled::buildTileProblem(model, 0);
    // The shared templates supply all mesh geometry and per-instance labels.
    problem->nodelist.clear();
    problem->linelist.clear();
    problem->arclist.clear();
    problem->labellist.clear();
    problem->pathName.clear();
    return problem;
}

std::string airGapName(const tiled::TiledModel &model)
{
    for (const auto &property : model.boundaryProps) {
        const auto *magnetic = dynamic_cast<const CMBoundaryProp *>(property.get());
        if (magnetic && (magnetic->BdryFormat == 6 || magnetic->BdryFormat == 7))
            return magnetic->BdryName;
    }
    throw std::runtime_error("pinned machine has no AGE boundary");
}

struct Machine {
    std::shared_ptr<FSolverAnalysisBackend> referenceBackend =
        std::make_shared<FSolverAnalysisBackend>();
    std::shared_ptr<FSolverAnalysisBackend> nativeBackend =
        std::make_shared<FSolverAnalysisBackend>();
    AnalysisSession reference;
    AnalysisSession native;

    Machine(const tiled::TiledModel &model, mesh::InstancedMesh instanced)
        : reference(ModelDefinition(sessionProblem(model)), referenceBackend),
          native(ModelDefinition(sessionProblem(model)), nativeBackend)
    {
        reference.setInstancedMesh(instanced);
        native.setInstancedMesh(std::move(instanced));
        native.setNativeInstanced(true);
    }
};

void currents(Machine &machine, const std::vector<double> &ampere)
{
    check(ampere.size() == 3, "machine fixture does not have three phases");
    for (std::size_t i = 0; i < ampere.size(); ++i) {
        const std::string name = std::to_string(i + 1);
        machine.reference.setCircuitCurrent(
            machine.reference.model().circuit(name), CComplex(ampere[i], 0.));
        machine.native.setCircuitCurrent(
            machine.native.model().circuit(name), CComplex(ampere[i], 0.));
    }
}

std::vector<std::size_t> indicesForBlock(const AnalysisSession &session, int block,
                                         std::size_t first, std::size_t last)
{
    const auto &labels = session.prepared().labels;
    check(last <= labels.size(), "invalid prepared machine region range");
    std::vector<std::size_t> selected;
    for (std::size_t i = first; i < last; ++i)
        if (labels[i].BlockType == block) selected.push_back(i);
    return selected;
}

void runMachine(const std::string &path, bool full, bool onePosition)
{
    tiled::TiledModel model;
    std::vector<tiled::TiledDiagnostic> diagnostics;
    check(tiled::loadTiledModelJson(readFile(path), model, diagnostics),
          "cannot parse pinned machine JSON" +
          (diagnostics.empty() ? std::string() : ": " + diagnostics.front().message));
    check(model.tiles.size() == 2, "expected rotor and stator machine tiles");
    check(model.circuitProps.size() == 3, "expected three machine phase circuits");
    // Same deterministic 60-degree reduction as instanced_solver_benchmark.
    sectorise(model, full ? 1u : 6u);
    const auto validation = tiled::validateTiledModel(model);
    check(validation.succeeded(), "invalid machine model" +
          (validation.diagnostics.empty() ? std::string() :
           ": " + validation.diagnostics.front().message));

    fmesher::TangleMesherBackend mesher;
    auto meshed = fmesher::meshTiledModel(model, mesher);
    check(meshed.ok, "machine tiling/meshing failed" +
          (meshed.diagnostics.empty() ? std::string() :
           ": " + meshed.diagnostics.front().message));
    check(meshed.instanced.templates.size() == 2, "machine must retain two mesh templates");
    const std::size_t expectedInstances =
        model.tiles[0].repeat.count + model.tiles[1].repeat.count;
    check(meshed.instanced.instances.size() == expectedInstances,
          "machine template instance count differs");

    Machine fixture(model, std::move(meshed.instanced));
    fixture.reference.synchronize();
    fixture.native.synchronize();
    const auto expanded = fixture.reference.mesh();
    const auto logical = fixture.native.logicalView();
    check(expanded && logical && logical->valid(), "missing machine reference/native mesh");
    check(expanded->nodes.size() == logical->nodeCount() &&
          expanded->elements.size() == logical->elementCount(),
          "machine topology counts mismatch before solve");
    check(expanded->airGaps.size() == 1 && logical->airGaps().size() == 1,
          "machine needs one coupled AGE");
    check(logical->storedByteCount() < logical->expandedByteCount(),
          "machine compressed view does not save stored mesh memory");
    check(fixture.native.materializationCount() == 0,
          "native machine mesh was materialized during setup");

    const auto initialMaterializations = fixture.reference.materializationCount();
    const auto initialViewGenerations = fixture.native.viewGenerationCount();
    const auto initialNativeOrdering = fixture.nativeBackend->nativeOrderingBuildCount();
    check(initialMaterializations == 1 && initialViewGenerations == 1,
          "expected one reference materialization and one native logical view");

    // Source labels are emitted per instance, in tile order. Select only
    // magnetic rotor regions for the mask, and stator coil regions for Lorentz.
    const std::size_t rotorEnd = model.tiles[0].geometry.labels.size() *
                                 model.tiles[0].repeat.count;
    const auto rotorMagnet = indicesForBlock(fixture.reference, 1, 0, rotorEnd);
    const auto rotorSteel = indicesForBlock(fixture.reference, 2, 0, rotorEnd);
    auto rotor = rotorMagnet;
    rotor.insert(rotor.end(), rotorSteel.begin(), rotorSteel.end());
    const auto coils = indicesForBlock(fixture.reference, 3, rotorEnd,
                                      fixture.reference.prepared().labels.size());
    check(!rotorMagnet.empty() && !rotorSteel.empty() && !coils.empty(),
          "machine fixture lost permanent magnets, rotor iron or phase coils");
    const std::string gap = airGapName(model);

    struct State { double angle; double ia, ib, ic; };
    const std::vector<State> states = onePosition ?
        std::vector<State>{{5., 0., 0., 0.}} :
        std::vector<State>{{0., 0., 0., 0.},
                           {5., 1.25, -0.625, -0.625},
                           {15., -1.25, 0.625, 0.625}};

    double maximumTorque = 0.;
    double maximumFlux = 0.;
    std::size_t solvedElements = 0;
    for (std::size_t position = 0; position < states.size(); ++position) {
        const auto &state = states[position];
        currents(fixture, {state.ia, state.ib, state.ic});
        fixture.reference.setAirGapAngle(
            fixture.reference.model().airGap(gap), state.angle, 0.);
        fixture.native.setAirGapAngle(
            fixture.native.model().airGap(gap), state.angle, 0.);

        xfemm_test::ForceTorqueOptions options;
        options.airGapBoundary = gap;
        // The fixture deliberately uses Lua magnetisation directions theta and
        // theta+180. G5.2's constitutive-H comparator currently rejects these
        // rather than evaluating Lua. B and all other physical checks remain ON.
        options.baseline.field.compareH = false;
        options.baseline.field.nodal.absoluteTolerance = 1e-5;
        options.baseline.field.nodal.relativeTolerance = 1e-5;
        options.baseline.field.matrixRelativeTolerance = 1e-5;
        options.baseline.field.rhsRelativeTolerance = 1e-5;
        options.baseline.field.maximumRelativeResidual = 2e-2;
        options.baseline.field.residualAgreementTolerance = 1e-4;
        options.baseline.field.fieldRelativeTolerance = 1e-4;
        options.baseline.fluxAbsoluteTolerance = 1e-8;
        options.baseline.fluxRelativeTolerance = 1e-4;
        options.baseline.energyRelativeTolerance = 1e-4;
        options.relativeTorqueTolerance = 5e-3;
        options.absoluteTorqueTolerance = 1e-5;
        // Each physical load is tested on the correct kind of region.
        // The full-machine test also exercises the rotor's independent mask.
        if (position == 0) {
            options.compareWeightedStress = true;
            options.selectedLabels = rotor;
        } else {
            options.compareLorentz = true;
            options.selectedLabels = coils;
        }

        const auto result = xfemm_test::compareForceTorqueSessions(
            fixture.reference, *fixture.referenceBackend,
            fixture.native, *fixture.nativeBackend, options);
        check(result.passed, "machine G5.5 angle=" + std::to_string(state.angle) +
                             " degrees: " + result.message);

        const auto &field = result.baseline.baseline;
        check(field.comparedElements == expanded->elements.size() &&
              field.matrix.numericNonzeros > 0 &&
              field.nodal.matchedNodes == expanded->nodes.size(),
              "machine comparison omitted elements, matrix entries or nodes");
        solvedElements = field.comparedElements;
        check(result.baseline.reference.circuits.size() == 3,
              "not all original three phase windings compared");
        for (const auto &circuit : result.baseline.reference.circuits)
            maximumFlux = std::max(maximumFlux, std::abs(circuit.fluxLinkage));
        maximumTorque = std::max(maximumTorque, std::abs(result.reference.airGapTorque));
        check(fixture.native.materializationCount() == 0,
              "G5.5 native path unexpectedly allocated expanded SolverMesh");
        check(fixture.native.viewGenerationCount() == initialViewGenerations &&
              fixture.nativeBackend->nativeOrderingBuildCount() == initialNativeOrdering,
              "machine angle/current change rebuilt native topology or ordering");
        check(fixture.reference.materializationCount() == initialMaterializations,
              "machine angle/current change rebuilt materialized topology");
        std::cout << "G5.5 " << (full ? "full" : "sector")
                  << " angle_deg=" << state.angle
                  << " currents_A=(" << state.ia << "," << state.ib << "," << state.ic
                  << ") A_max_delta=" << field.nodal.maxAbsoluteError
                  << " B_max_delta_T=" << field.maximumBError
                  << " nnz=" << field.matrix.numericNonzeros
                  << " max_flux_delta_Wb=" << result.baseline.maximumFluxDifference
                  << " max_force_delta_N=" << result.maximumForceDifference
                  << " max_torque_delta_Nm=" << result.maximumTorqueDifference
                  << " AGE_torque_Nm=" << result.reference.airGapTorque << "\n";
    }
    check(solvedElements > 0 && maximumFlux > 1e-12,
          "machine fixture compared a trivial field/circuit linkage");
    check(maximumTorque > 1e-8, "machine AGE torque was trivial across all positions");
    check(fixture.nativeBackend->nativeSolveCount() == states.size() &&
          fixture.referenceBackend->solveCount() == states.size(),
          "G5.5 performed extra or missing physics solves");
    std::cout << "G5.5 PASS fixture=RNFoundry-f4d42805"
              << " sector_divisor=" << (full ? 1 : 6)
              << " instances=" << expectedInstances
              << " nodes=" << logical->nodeCount()
              << " elements=" << logical->elementCount()
              << " positions=" << states.size()
              << " native_stored_bytes=" << logical->storedByteCount()
              << " materialized_equivalent_bytes=" << logical->expandedByteCount()
              << "\n";
}
} // namespace

int main(int argc, char **argv)
{
    if (argc < 2 || argc > 4) {
        std::cerr << "usage: native_machine_g5_5_test tiled_machine.json [--full] [--single]\n";
        return 2;
    }
    bool full = false, single = false;
    for (int i = 2; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--full") full = true;
        else if (arg == "--single") single = true;
        else {
            std::cerr << "invalid G5.5 option: " << arg << "\n";
            return 2;
        }
    }
    try {
        runMachine(argv[1], full, single);
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "G5.5 FAILED: " << e.what() << "\n";
        return 1;
    }
}
