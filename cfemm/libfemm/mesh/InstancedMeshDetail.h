#ifndef FEMM_MESH_INSTANCEDMESHDETAIL_H
#define FEMM_MESH_INSTANCEDMESHDETAIL_H

#include "InstancedMesh.h"

#include <cstddef>
#include <functional>
#include <vector>

namespace femm {
namespace mesh {
namespace detail {

/**
 * Result of validating an InstancedMesh and resolving its declared seam welds.
 *
 * The plan is the shared basis of both InstancedMesh::materialize() and
 * LogicalMeshView so the two cannot disagree about which logical nodes are
 * welded or in which order global nodes are first encountered.
 */
struct WeldPlan {
    /** Flattened index of each instance's first local node. */
    std::vector<std::size_t> instanceBase;
    /** Total number of logical nodes before welding. */
    std::size_t totalNodes = 0;
    /**
     * Representative flattened index for every flattened logical node. Welded
     * nodes share a representative; the representative is the first flattened
     * index in instance-then-local order that belongs to the weld class.
     */
    std::vector<std::size_t> root;
};

/**
 * Validate the instanced mesh, apply every unique declared seam connection, and
 * verify that welded nodes actually coincide. Returns false and fills
 * diagnostics when any invariant fails; the plan is only meaningful on success.
 */
bool buildWeldPlan(const InstancedMesh &instanced, WeldPlan &plan,
                   std::vector<MaterializationDiagnostic> &diagnostics);

/** Remapped ordinary periodic constraints and AGE topology of a materialisation. */
struct TopologyRemap {
    std::vector<SolverMesh::PeriodicConstraint> periodicConstraints;
    std::vector<SolverMesh::AirGap> airGaps;
};

/**
 * Remap every local periodic constraint, periodic closure, per-instance AGE,
 * and cross-template AGE coupling through the supplied node and coordinate
 * accessors. The accessors must already honour welded seams. Returns false and
 * fills diagnostics when any referenced local node is invalid or an AGE
 * structure is inconsistent.
 */
bool remapTopology(const InstancedMesh &instanced,
                   const std::function<MeshIndex(std::size_t, MeshIndex)> &nodeFor,
                   const std::function<void(MeshIndex, double &, double &)> &coordinates,
                   TopologyRemap &out,
                   std::vector<MaterializationDiagnostic> &diagnostics);

} // namespace detail
} // namespace mesh
} // namespace femm

#endif
