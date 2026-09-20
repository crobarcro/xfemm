#include "AnalysisSession.h"
#include "FSolverAnalysisBackend.h"

#include "TangleMesherBackend.h"
#include "TiledModel.h"
#include "TiledModelMesher.h"

#include "mesh/LogicalMeshView.h"
#include "mesh/SolverMeshValidator.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace {

using namespace femm::mesh;

constexpr double Pi = 3.141592653589793238462643383279502884;

int fail(const std::string &message)
{
    std::cerr << message << '\n';
    return 1;
}

femm::tiled::Tile annularTile(const std::string &name, double innerRadius,
                              double outerRadius, double sectorDegrees,
                              std::size_t count, femm::tiled::Closure closure,
                              bool airGapOnOuter)
{
    const double sector = sectorDegrees * Pi / 180.0;
    femm::tiled::Tile tile;
    tile.name = name;
    tile.repeat.count = count;
    tile.repeat.closure = closure;
    tile.seamBoundary = "seam";
    tile.airGapBoundary = "gap";
    tile.geometry.nodes = {
        {innerRadius, 0.0, -1, 0},
        {outerRadius, 0.0, -1, 0},
        {innerRadius * std::cos(sector), innerRadius * std::sin(sector), -1, 0},
        {outerRadius * std::cos(sector), outerRadius * std::sin(sector), -1, 0},
    };
    tile.geometry.segments = {{0, 1, -1, 0, false, 0}, {2, 3, -1, 0, false, 0}};
    // Boundary 1 = "outer" (non-physical), boundary 2 = "gap" (AGE).
    const int outerArcBoundary = airGapOnOuter ? 2 : 1;
    const int innerArcBoundary = airGapOnOuter ? 1 : 2;
    tile.geometry.arcs = {
        {1, 3, sectorDegrees, 2.0, outerArcBoundary, false, 0},
        {2, 0, sectorDegrees, 2.0, innerArcBoundary, false, 0},
    };
    tile.geometry.labels = {{"coil", 0.5 * (innerRadius + outerRadius),
                             0.5 * (innerRadius + outerRadius) * 0.1, 0, 0, 1, 0.0, "",
                             1e-4, 0, false}};
    return tile;
}

femm::tiled::TiledModel makeTwoTileModel()
{
    femm::tiled::TiledModel model;
    model.depth = 0.1;

    auto seam = std::make_unique<femm::CMBoundaryProp>();
    seam->BdryName = "seam";
    seam->BdryFormat = 4; // periodic
    model.boundaryProps.push_back(std::move(seam));

    auto outer = std::make_unique<femm::CMBoundaryProp>();
    outer->BdryName = "outer";
    outer->BdryFormat = 0;
    model.boundaryProps.push_back(std::move(outer));

    auto gap = std::make_unique<femm::CMBoundaryProp>();
    gap->BdryName = "gap";
    gap->BdryFormat = 6; // periodic AGE
    model.boundaryProps.push_back(std::move(gap));

    auto material = std::make_unique<femm::CMMaterialProp>();
    material->BlockName = "steel";
    material->mu_x = material->mu_y = 1000;
    model.materialProps.push_back(std::move(material));

    auto circuit = std::make_unique<femm::CMCircuit>();
    circuit->CircName = "phase";
    circuit->CircType = 1;
    circuit->Amps = 1;
    model.circuitProps.push_back(std::move(circuit));

    // Rotor pole (30 deg, air gap on its outer arc) and stator slot (10 deg,
    // air gap on its inner arc); both cover the 60-degree sector.
    model.tiles.push_back(annularTile("rotor", 0.04, 0.05, 30.0, 2,
                                      femm::tiled::Closure::Periodic, true));
    model.tiles.push_back(annularTile("stator", 0.06, 0.08, 10.0, 6,
                                      femm::tiled::Closure::Periodic, false));
    model.couplings.push_back({"rotor", "stator", "gap", 0.0, 0.0, 0.05, 0.06});
    return model;
}

int runPipelineChecks(femm::mesh::InstancedMesh &instanced)
{
    if (instanced.templates.size() != 2)
        return fail("expected two templates");
    if (instanced.instances.size() != 8)
        return fail("expected eight instances (2 rotor + 6 stator)");
    if (instanced.airGapCouplings.size() != 1)
        return fail("expected one air-gap coupling");
    if (instanced.periodicClosures.size() != 2)
        return fail("expected one periodic closure per open tile");

    const auto materialized = instanced.materialize();
    if (!materialized.succeeded()) {
        for (const auto &diagnostic : materialized.diagnostics)
            std::cerr << "  materialize: " << diagnostic.message << '\n';
        return fail("materialisation failed");
    }
    if (!validateSolverMesh(materialized.mesh).valid())
        return fail("materialised mesh is invalid");
    if (materialized.mesh.airGaps.size() != 1)
        return fail("expected one materialised AGE");
    const auto &age = materialized.mesh.airGaps.front();
    if (age.boundaryName != "gap")
        return fail("AGE has the wrong boundary name");
    if (age.innerRing.empty() || age.innerRing.size() != age.outerRing.size())
        return fail("AGE rings have incompatible cardinality");
    // The two domains are coupled, not welded: no rotor air-gap node is a
    // stator air-gap node.
    for (const auto &inner : age.innerRing)
        for (const auto &outer : age.outerRing)
            if (inner.node == outer.node)
                return fail("AGE rings share a node (domains were welded)");
    return 0;
}

/**
 * G4/G5: the native compressed assembly must reproduce the materialised AGE
 * solve, including after the AGE is repositioned without remeshing.
 */
int compareNative(femm::AnalysisSession &session,
                  const std::shared_ptr<femm::FSolverAnalysisBackend> &solver,
                  const std::string &label)
{
    const auto &materializedSystem = solver->solvedSystem();
    const std::size_t materializedNodes = solver->solvedSolver().NumNodes;
    const std::vector<double> reference(materializedSystem.rhs().begin(),
                                        materializedSystem.rhs().end());
    // Tangle's synchronised seam splitting can place a welded node a few
    // nanometres from the representative, so match coordinates at a
    // 1e-3 cm (1e-5 m) tolerance rather than exact equality.
    std::map<std::pair<long long, long long>, std::size_t> lookup;
    for (std::size_t i = 0; i < solver->solvedSolver().meshnode.size(); ++i) {
        lookup[{std::llround(solver->solvedSolver().meshnode[i].x * 1e3),
                std::llround(solver->solvedSolver().meshnode[i].y * 1e3)}] = i;
    }

    std::vector<MaterializationDiagnostic> diagnostics;
    const LogicalMeshView view = LogicalMeshView::build(*session.instancedMesh(), diagnostics);
    if (!view.valid())
        return fail(label + ": logical view is invalid");
    if (view.nodeCount() != materializedNodes)
        return fail(label + ": node count differs");
    if (!solver->solveNative(view, session.instanceLabelBases(), session.airGapPositioning()))
        return fail(label + ": native solve failed");

    const auto canonical = session.mesh();
    for (std::size_t node = 0; node < view.nodeCount(); ++node) {
        double x = 0.0, y = 0.0;
        view.nodeCoordinates(node, x, y);
        if (std::abs(x - canonical->nodes[node].x) > 1e-6 ||
            std::abs(y - canonical->nodes[node].y) > 1e-6) {
            std::cerr << "view/materialised coordinate mismatch at " << node << ": view=" << x
                      << "," << y << " mesh=" << canonical->nodes[node].x << ","
                      << canonical->nodes[node].y << '\n';
            return fail(label + ": view coordinates differ from the materialised mesh");
        }
    }

    const auto &native = solver->nativeSystem();
    double largest = 0.0;
    for (std::size_t node = 0; node < view.nodeCount(); ++node) {
        double x = 0.0, y = 0.0;
        view.nodeCoordinates(node, x, y);
        const auto found =
            lookup.find({std::llround(x * 100.0 * 1e3), std::llround(y * 100.0 * 1e3)});
        if (found == lookup.end())
            return fail(label + ": node is missing from the materialised mesh");
        const double nativeValue = native.rhs()[node];
        const double referenceValue = reference[found->second];
        largest = std::max(largest, std::abs(referenceValue));
        if (!(std::abs(nativeValue - referenceValue) < 1e-6 * (1.0 + std::abs(referenceValue))))
            return fail(label + ": native field differs from the materialised field");
    }
    if (largest <= 1e-9)
        return fail(label + ": trivial field");
    return 0;
}

} // namespace

int main()
{
    femm::tiled::TiledModel model = makeTwoTileModel();
    fmesher::TangleMesherBackend backend;
    fmesher::TiledMeshResult meshed = fmesher::meshTiledModel(model, backend);
    if (!meshed.ok) {
        for (const auto &diagnostic : meshed.diagnostics)
            std::cerr << "  mesh: " << diagnostic.message << '\n';
        return fail("meshTiledModel failed");
    }
    if (const int status = runPipelineChecks(meshed.instanced))
        return status;

    // The session accepts the multi-tile instanced mesh, materialises it, and
    // reuses topology across AGE positions.
    std::unique_ptr<femm::FemmProblem> sessionModel = femm::tiled::buildTileProblem(model, 0);
    sessionModel->nodelist.clear();
    sessionModel->linelist.clear();
    sessionModel->arclist.clear();
    sessionModel->labellist.clear();
    sessionModel->pathName.clear();

    // Keep a copy for the session-native AGE-angle check below; the original is
    // moved into the materialised session.
    femm::mesh::InstancedMesh nativeInstanced = meshed.instanced;

    auto solver = std::make_shared<femm::FSolverAnalysisBackend>();
    femm::AnalysisSession session(femm::ModelDefinition(std::move(sessionModel)), solver);
    session.setInstancedMesh(std::move(meshed.instanced));

    session.synchronize();
    const auto mesh = session.mesh();
    if (!mesh || mesh->airGaps.size() != 1)
        return fail("session did not materialise the coupled AGE");
    if (session.prepared().labels.size() != 8)
        return fail("session did not build one label set per instance");
    const std::size_t materializations = session.materializationCount();

    session.solve();
    if (const int status = compareNative(session, solver, "native-default"))
        return status;

    // A relative rotor position updates the AGE coupling without remeshing the
    // stored templates; the native path must track it.
    session.setAirGapAngle(session.model().airGap("gap"), 15.0, 0.0);
    session.solve();
    if (const int status = compareNative(session, solver, "native-angle"))
        return status;
    if (session.materializationCount() != materializations)
        return fail("changing the AGE angle rematerialised the mesh");

    // Session-native mode must reuse the topology-derived ordering across a
    // physics-only AGE angle change.
    {
        std::unique_ptr<femm::FemmProblem> nativeModel = femm::tiled::buildTileProblem(model, 0);
        nativeModel->nodelist.clear();
        nativeModel->linelist.clear();
        nativeModel->arclist.clear();
        nativeModel->labellist.clear();
        nativeModel->pathName.clear();
        auto nativeSolver = std::make_shared<femm::FSolverAnalysisBackend>();
        femm::AnalysisSession nativeSession(femm::ModelDefinition(std::move(nativeModel)),
                                            nativeSolver);
        nativeSession.setInstancedMesh(std::move(nativeInstanced));
        nativeSession.setNativeInstanced(true);
        nativeSession.solve();
        if (nativeSolver->nativeOrderingBuildCount() != 1)
            return fail("native ordering was not built exactly once");
        nativeSession.setAirGapAngle(nativeSession.model().airGap("gap"), 15.0, 0.0);
        nativeSession.solve();
        if (nativeSolver->nativeOrderingBuildCount() != 1)
            return fail("AGE angle change rebuilt the native ordering");
        if (nativeSession.viewGenerationCount() != 1)
            return fail("AGE angle change rebuilt the native view");
    }

    for (double angle : {15.0, 30.0}) {
        session.setAirGapAngle(session.model().airGap("gap"), angle, 0.0);
        session.ensureMesh();
    }
    if (session.materializationCount() != materializations)
        return fail("changing the AGE angle rematerialised the mesh");

    std::cout << "tiled model air-gap coupling tests passed\n";
    return 0;
}
