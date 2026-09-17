#include "TiledModelMesher.h"

#include "MesherBackend.h"

#include "FemmProblem.h"
#include "mesh/SolverMesh.h"

#include <algorithm>
#include <cmath>
#include <set>
#include <string>

namespace fmesher {
namespace {

using femm::mesh::InstancedMesh;
using femm::mesh::MeshDiagnostic;
using femm::mesh::MeshDiagnosticSeverity;
using femm::mesh::MeshIndex;
using femm::mesh::MeshTemplate;
using femm::mesh::SolverMesh;

constexpr double Pi = 3.141592653589793238462643383279502884;

MeshDiagnostic error(const std::string &message)
{
    return {MeshDiagnosticSeverity::Error, message, "TiledModelMesher", 0};
}

void addDiagnostics(TiledMeshResult &result,
                    const std::vector<femm::tiled::TiledDiagnostic> &diagnostics)
{
    for (const auto &diagnostic : diagnostics)
        result.diagnostics.push_back(
            {MeshDiagnosticSeverity::Error, diagnostic.message, "TiledModel", 0});
}

bool chainsMatch(const SolverMesh &mesh, const std::vector<MeshIndex> &from,
                 const std::vector<MeshIndex> &to, double centerX, double centerY,
                 double angleDegrees)
{
    if (from.empty() || from.size() != to.size())
        return false;
    const femm::mesh::RigidTransform2D transform =
        femm::mesh::RigidTransform2D::rotationAbout(centerX, centerY, angleDegrees);
    for (std::size_t i = 0; i < from.size(); ++i) {
        if (from[i] >= mesh.nodes.size() || to[i] >= mesh.nodes.size())
            return false;
        double x = 0.0, y = 0.0;
        transform.applyPoint(mesh.nodes[from[i]].x, mesh.nodes[from[i]].y, x, y);
        const auto &target = mesh.nodes[to[i]];
        if (std::abs(x - target.x) > femm::mesh::InstancedMeshWeldToleranceMetres ||
            std::abs(y - target.y) > femm::mesh::InstancedMeshWeldToleranceMetres)
            return false;
    }
    return true;
}

/** Order a radial seam's nodes from the rotation centre outwards. */
std::vector<MeshIndex> orderByRadius(const SolverMesh &mesh,
                                     const std::vector<MeshIndex> &nodes, double centerX,
                                     double centerY)
{
    std::vector<MeshIndex> ordered = nodes;
    const auto radiusSquared = [&](MeshIndex node) {
        const double dx = mesh.nodes[node].x - centerX;
        const double dy = mesh.nodes[node].y - centerY;
        return dx * dx + dy * dy;
    };
    std::stable_sort(ordered.begin(), ordered.end(),
                     [&](MeshIndex a, MeshIndex b) { return radiusSquared(a) < radiusSquared(b); });
    return ordered;
}

/**
 * Ordered nodes on the tile's air-gap boundary, by angle about the rotation
 * centre. The boundary is a radial arc, so angular order is also edge order.
 */
std::vector<MeshIndex> airGapChain(const SolverMesh &mesh, int boundaryMarker,
                                   double centerX, double centerY)
{
    // Tangle encodes a boundary property index b as the segment marker
    // -(b + 2); the converter preserves that marker on SolverMesh edges.
    const std::int32_t marker = -static_cast<std::int32_t>(boundaryMarker + 2);
    std::set<MeshIndex> nodes;
    for (const auto &edge : mesh.edges)
        if (edge.boundaryMarker == marker) {
            nodes.insert(edge.first);
            nodes.insert(edge.second);
        }
    std::vector<MeshIndex> chain(nodes.begin(), nodes.end());
    std::stable_sort(chain.begin(), chain.end(), [&](MeshIndex a, MeshIndex b) {
        const auto angle = [&](MeshIndex node) {
            double value = std::atan2(mesh.nodes[node].y - centerY,
                                      mesh.nodes[node].x - centerX) *
                           180.0 / Pi;
            if (value < 0.0)
                value += 360.0;
            return value;
        };
        return angle(a) < angle(b);
    });
    return chain;
}

} // namespace

TiledMeshResult meshTiledModel(const femm::tiled::TiledModel &model, MesherBackend &mesher,
                               const femm::mesh::MeshingOptions &options)
{
    TiledMeshResult result;

    const femm::tiled::TiledValidationResult validation =
        femm::tiled::validateTiledModel(model);
    if (!validation.succeeded()) {
        addDiagnostics(result, validation.diagnostics);
        return result;
    }

    // Mesh each tile once and build its template.
    for (std::size_t tileIndex = 0; tileIndex < model.tiles.size(); ++tileIndex) {
        const femm::tiled::Tile &tile = model.tiles[tileIndex];
        const std::unique_ptr<femm::FemmProblem> problem =
            femm::tiled::buildTileProblem(model, tileIndex);
        if (!problem) {
            result.diagnostics.push_back(error("could not build tile problem"));
            return result;
        }
        femm::mesh::MeshingRequest request;
        request.options = options;
        request.createPeriodicFieldConstraints = false;
        femm::mesh::MeshResult meshed = mesher.mesh(*problem, request);
        if (!meshed.succeeded()) {
            for (const auto &diagnostic : meshed.diagnostics)
                result.diagnostics.push_back(diagnostic);
            return result;
        }

        const int seamIndex = femm::tiled::boundaryIndexByName(model, tile.seamBoundary);
        const femm::mesh::MeshBoundaryMatch *match = nullptr;
        for (const auto &candidate : meshed.boundaryMatches)
            if (static_cast<int>(candidate.boundaryProperty) == seamIndex) {
                match = &candidate;
                break;
            }
        if (!match || match->firstNodes.empty() ||
            match->firstNodes.size() != match->secondNodes.size()) {
            result.diagnostics.push_back(
                error("tile '" + tile.name + "' has no usable matched seam pair"));
            return result;
        }

        MeshTemplate meshTemplate;
        meshTemplate.localMesh = std::move(meshed.mesh);
        meshTemplate.localMesh.periodicConstraints.clear();
        meshTemplate.localMesh.airGaps.clear();

        const double step = validation.tiles[tileIndex].totalAngleDegrees /
                            static_cast<double>(tile.repeat.count);
        const std::vector<MeshIndex> first =
            orderByRadius(meshTemplate.localMesh, match->firstNodes,
                          tile.repeat.centerXMetres, tile.repeat.centerYMetres);
        const std::vector<MeshIndex> second =
            orderByRadius(meshTemplate.localMesh, match->secondNodes,
                          tile.repeat.centerXMetres, tile.repeat.centerYMetres);
        std::vector<MeshIndex> leadNodes;
        std::vector<MeshIndex> trailNodes;
        if (chainsMatch(meshTemplate.localMesh, first, second, tile.repeat.centerXMetres,
                        tile.repeat.centerYMetres, step)) {
            leadNodes = first;
            trailNodes = second;
        } else if (chainsMatch(meshTemplate.localMesh, second, first,
                               tile.repeat.centerXMetres, tile.repeat.centerYMetres, step)) {
            leadNodes = second;
            trailNodes = first;
        } else {
            result.diagnostics.push_back(
                error("tile '" + tile.name + "' seams are not related by the pitch"));
            return result;
        }
        meshTemplate.seams.push_back({"lead", std::move(leadNodes)});
        meshTemplate.seams.push_back({"trail", std::move(trailNodes)});

        const int gapIndex = femm::tiled::boundaryIndexByName(model, tile.airGapBoundary);
        if (gapIndex >= 0) {
            const std::vector<MeshIndex> chain = airGapChain(
                meshTemplate.localMesh, gapIndex, tile.repeat.centerXMetres,
                tile.repeat.centerYMetres);
            meshTemplate.seams.push_back({"airgap", chain});
        }

        result.instanced.templates.push_back(std::move(meshTemplate));
    }

    // Place the instances, weld internal joins, and close open sectors.
    std::vector<std::size_t> instanceBase(model.tiles.size(), 0);
    for (std::size_t tileIndex = 0; tileIndex < model.tiles.size(); ++tileIndex) {
        const femm::tiled::Tile &tile = model.tiles[tileIndex];
        const std::size_t base = result.instanced.instances.size();
        const double step = validation.tiles[tileIndex].totalAngleDegrees /
                            static_cast<double>(tile.repeat.count);
        for (std::size_t k = 0; k < tile.repeat.count; ++k) {
            femm::mesh::MeshInstance instance;
            instance.templateIndex = tileIndex;
            instance.transform = femm::mesh::RigidTransform2D::rotationAbout(
                tile.repeat.centerXMetres, tile.repeat.centerYMetres,
                static_cast<double>(k) * step);
            if (k + 1 < tile.repeat.count)
                instance.seamConnections.push_back(
                    {1, base + k + 1, 0, femm::mesh::SeamOrientation::Forward});
            result.instanced.instances.push_back(std::move(instance));
        }
        if (tile.repeat.closure == femm::tiled::Closure::Closed) {
            result.instanced.instances[base + tile.repeat.count - 1]
                .seamConnections.push_back(
                    {1, base, 0, femm::mesh::SeamOrientation::Forward});
        } else {
            femm::mesh::PeriodicClosure closure;
            closure.firstInstance = base + tile.repeat.count - 1;
            closure.firstSeam = 1; // trail
            closure.secondInstance = base;
            closure.secondSeam = 0; // lead
            closure.orientation = femm::mesh::SeamOrientation::Forward;
            closure.periodicity =
                tile.repeat.closure == femm::tiled::Closure::Antiperiodic
                    ? SolverMesh::Periodicity::Antiperiodic
                    : SolverMesh::Periodicity::Periodic;
            result.instanced.periodicClosures.push_back(closure);
        }
    }

    // Cross-tile air-gap couplings.
    for (const femm::tiled::TileCoupling &coupling : model.couplings) {
        const auto tileIndexOf = [&](const std::string &name) -> std::size_t {
            for (std::size_t i = 0; i < model.tiles.size(); ++i)
                if (model.tiles[i].name == name)
                    return i;
            return model.tiles.size();
        };
        const std::size_t inner = tileIndexOf(coupling.innerTile);
        const std::size_t outer = tileIndexOf(coupling.outerTile);
        if (inner >= model.tiles.size() || outer >= model.tiles.size()) {
            result.diagnostics.push_back(error("coupling references an unknown tile"));
            return result;
        }
        femm::mesh::AirGapCoupling converted;
        converted.innerTemplate = inner;
        converted.innerSeam = 2; // the air-gap seam appended above
        converted.outerTemplate = outer;
        converted.outerSeam = 2;
        converted.boundaryName = coupling.boundary;
        converted.centerXMetres = coupling.centerXMetres;
        converted.centerYMetres = coupling.centerYMetres;
        converted.innerRadiusMetres = coupling.innerRadiusMetres;
        converted.outerRadiusMetres = coupling.outerRadiusMetres;
        converted.totalArcLengthDegrees =
            validation.tiles[outer].totalAngleDegrees;
        result.instanced.airGapCouplings.push_back(converted);
    }

    result.ok = true;
    return result;
}

} // namespace fmesher
