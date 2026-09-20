#ifndef FEMM_FMESHER_TILEDMODELMESHER_H
#define FEMM_FMESHER_TILEDMODELMESHER_H

#include "TiledModel.h"
#include "mesh/InstancedMesh.h"
#include "mesh/Meshing.h"

#include <vector>

namespace fmesher {

class MesherBackend;

/**
 * Result of turning a validated TiledModel into an InstancedMesh: one meshed
 * template per tile, the instances produced by the repeat rule, the internal
 * welds, the open-sector periodic closure, and the cross-tile air-gap
 * couplings. The instanced mesh is not yet materialised; call materialize().
 */
struct TiledMeshResult {
    bool ok = false;
    femm::mesh::InstancedMesh instanced;
    std::vector<femm::mesh::MeshDiagnostic> diagnostics;
};

/**
 * Mesh every tile in \p model exactly once and assemble the instanced mesh.
 * Each tile is meshed from its in-memory FemmProblem through \p mesher; the
 * tile's declared periodic seam is used for topology-only boundary matching,
 * so no field constraint is emitted.
 */
TiledMeshResult meshTiledModel(const femm::tiled::TiledModel &model, MesherBackend &mesher,
                               const femm::mesh::MeshingOptions &options = {});

} // namespace fmesher

#endif
