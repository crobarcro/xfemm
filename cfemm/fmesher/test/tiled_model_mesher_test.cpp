#include "TiledModelMesher.h"
#include "TangleMesherBackend.h"

#include "TiledModel.h"
#include "mesh/SolverMeshValidator.h"

#include <tangle_mesh.h>

#include <cmath>
#include <iostream>
#include <memory>
#include <string>

namespace {

using namespace femm::mesh;

constexpr double Pi = 3.141592653589793238462643383279502884;

int fail(const std::string &message)
{
    std::cerr << message << '\n';
    return 1;
}

femm::tiled::TiledModel makeModel(double sectorDegrees, std::size_t count,
                                  femm::tiled::Closure closure, bool withAirGap)
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

    if (withAirGap) {
        auto gap = std::make_unique<femm::CMBoundaryProp>();
        gap->BdryName = "sliding-gap";
        gap->BdryFormat = 6; // periodic AGE
        model.boundaryProps.push_back(std::move(gap));
    }

    auto material = std::make_unique<femm::CMMaterialProp>();
    material->BlockName = "steel";
    material->mu_x = material->mu_y = 1000;
    model.materialProps.push_back(std::move(material));

    auto circuit = std::make_unique<femm::CMCircuit>();
    circuit->CircName = "phase";
    circuit->CircType = 1;
    model.circuitProps.push_back(std::move(circuit));

    const double sector = sectorDegrees * Pi / 180.0;
    const double ri = 0.05;
    const double ro = 0.06;
    femm::tiled::Tile tile;
    tile.name = "slot";
    tile.repeat.count = count;
    tile.repeat.closure = closure;
    tile.seamBoundary = "seam";
    if (withAirGap)
        tile.airGapBoundary = "sliding-gap";
    tile.geometry.nodes = {
        {ri, 0.0, -1, 0},
        {ro, 0.0, -1, 0},
        {ri * std::cos(sector), ri * std::sin(sector), -1, 0},
        {ro * std::cos(sector), ro * std::sin(sector), -1, 0},
    };
    tile.geometry.segments = {{0, 1, -1, 0, false, 0}, {2, 3, -1, 0, false, 0}};
    tile.geometry.arcs = {
        {1, 3, sectorDegrees, 2.0, 1, false, 0},
        {2, 0, sectorDegrees, 2.0, 1, false, 0},
    };
    tile.geometry.labels = {
        {"coil", 0.5 * (ri + ro), 0.004, 0, 0, 1, 0.0, "", 1e-4, 0, false}};
    model.tiles.push_back(tile);
    return model;
}

int runCase(const femm::tiled::TiledModel &model, const std::string &label,
            int expectedTemplates, int expectedInstances, bool expectClosure)
{
    int inMemoryCalls = 0;
    fmesher::TangleMesherBackend backend(
        {},
        [&](const ::FemProblem &problem, const ::MeshOptions &options, ::Mesh &mesh) {
            ++inMemoryCalls;
            return tangle_mesh_fem(problem, options, mesh);
        });
    fmesher::TiledMeshResult result = fmesher::meshTiledModel(model, backend);
    if (!result.ok) {
        for (const auto &diagnostic : result.diagnostics)
            std::cerr << "  " << label << ": " << diagnostic.message << '\n';
        return fail(label + ": meshTiledModel failed");
    }
    if (inMemoryCalls != expectedTemplates)
        return fail(label + ": wrong number of tile meshing calls");
    if (static_cast<int>(result.instanced.templates.size()) != expectedTemplates)
        return fail(label + ": wrong template count");
    if (static_cast<int>(result.instanced.instances.size()) != expectedInstances)
        return fail(label + ": wrong instance count");
    if (expectClosure != !result.instanced.periodicClosures.empty())
        return fail(label + ": unexpected periodic closure state");

    const auto materialized = result.instanced.materialize();
    if (!materialized.succeeded()) {
        for (const auto &diagnostic : materialized.diagnostics)
            std::cerr << "  " << label << ": " << diagnostic.message << '\n';
        return fail(label + ": materialisation failed");
    }
    if (!validateSolverMesh(materialized.mesh).valid())
        return fail(label + ": materialised mesh is invalid");
    if (expectClosure && materialized.mesh.periodicConstraints.empty())
        return fail(label + ": open sector has no periodic constraints");
    return 0;
}

int testOverridesApplied()
{
    femm::tiled::TiledModel model = makeModel(60.0, 6, femm::tiled::Closure::Closed, false);
    femm::tiled::TileLabelOverride override;
    override.tile = "slot";
    override.label = "coil";
    override.circuit = {-1, -1, 0, -1, -1, -1};
    override.turnScale = {1.0, 1.0, -2.0, 1.0, 1.0, 1.0};
    override.magDir = {0.0, 0.0, 45.0, 0.0, 0.0, 0.0};
    model.overrides.push_back(override);

    fmesher::TangleMesherBackend backend;
    fmesher::TiledMeshResult result = fmesher::meshTiledModel(model, backend);
    if (!result.ok) {
        for (const auto &diagnostic : result.diagnostics)
            std::cerr << "  override case: " << diagnostic.message << '\n';
        return fail("override case: meshTiledModel failed");
    }
    const auto &instances = result.instanced.instances;
    if (instances.size() != 6)
        return fail("override case: wrong instance count");
    for (std::size_t k = 0; k < instances.size(); ++k) {
        const auto &overrides = instances[k].regionOverrides;
        if (overrides.size() != 1)
            return fail("override case: instance has wrong override count");
        const auto &region = overrides.front();
        if (region.sourceBlockLabel != 0)
            return fail("override case: wrong source label");
        if (k == 2) {
            if (!region.circuit || *region.circuit != 0)
                return fail("override case: circuit override not applied");
            if (!region.currentScale || std::abs(*region.currentScale + 2.0) > 1e-12)
                return fail("override case: turn scale override not applied");
            if (!region.magnetisationRotationDegrees ||
                std::abs(*region.magnetisationRotationDegrees - 45.0) > 1e-12)
                return fail("override case: magnetisation override not applied");
        } else {
            if (region.circuit)
                return fail("override case: unused circuit entry became active");
            if (region.currentScale && std::abs(*region.currentScale - 1.0) > 1e-12)
                return fail("override case: unused turn scale entry changed");
            if (region.magnetisationRotationDegrees &&
                std::abs(*region.magnetisationRotationDegrees) > 1e-12)
                return fail("override case: unused magnetisation entry changed");
        }
    }
    return 0;
}

} // namespace

int main()
{
    if (const int status = runCase(makeModel(10.0, 36, femm::tiled::Closure::Closed, false),
                                   "closed ring", 1, 36, false))
        return status;
    if (const int status = runCase(makeModel(10.0, 6, femm::tiled::Closure::Periodic, false),
                                   "open sector", 1, 6, true))
        return status;
    if (const int status = testOverridesApplied())
        return status;
    std::cout << "tiled model mesher tests passed\n";
    return 0;
}
