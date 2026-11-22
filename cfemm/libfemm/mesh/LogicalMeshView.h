#ifndef FEMM_MESH_LOGICALMESHVIEW_H
#define FEMM_MESH_LOGICALMESHVIEW_H

#include "InstancedMesh.h"
#include "SolverMesh.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <map>
#include <utility>
#include <vector>

namespace femm {
namespace mesh {

/**
 * A compact, solver-facing view of an InstancedMesh.
 *
 * The view presents the same logical mesh that InstancedMesh::materialize()
 * would produce, but it stores template geometry once plus per-instance
 * metadata instead of a fully expanded element list. Logical element
 * connectivity is computed from the template and the instance transform, and
 * welded seam nodes resolve to a shared global node. Unconnected instance nodes
 * keep independent global indices.
 *
 * The global node numbering is deliberately identical to materialize(): a
 * LogicalMeshView can be compared directly against the materialised oracle.
 *
 * Physics metadata (circuits, magnetisation, per-instance overrides) is not
 * part of the view; a consumer resolves region attributes per instance using
 * elementProvenance().
 */
class LogicalMeshView {
public:
    LogicalMeshView() = default;

    /**
     * Validate and build the view. On failure `diagnostics` holds every
     * rejection and valid() is false; the view is otherwise unusable.
     */
    static LogicalMeshView build(const InstancedMesh &instanced,
                                 std::vector<MaterializationDiagnostic> &diagnostics);

    bool valid() const { return m_valid; }
    const std::vector<MaterializationDiagnostic> &diagnostics() const { return m_diagnostics; }

    std::size_t nodeCount() const { return m_nodeCount; }
    std::size_t elementCount() const { return m_elementCount; }
    std::size_t instanceCount() const { return m_instances.size(); }
    /** Number of unique undirected edges; computed from element connectivity. */
    std::size_t edgeCount() const;

    void nodeCoordinates(MeshIndex node, double &x, double &y) const;
    /** Fill x/y with every global node coordinate in global node order. */
    void allNodeCoordinates(std::vector<double> &x, std::vector<double> &y) const;
    /** Raw template boundary marker of a global node (see SolverMesh::Node). */
    std::int32_t nodeBoundaryMarker(MeshIndex node) const;
    std::array<MeshIndex, 3> elementNodes(MeshIndex element) const;
    /** Element corner coordinates in metres, computed without reverse provenance. */
    void elementCoordinates(MeshIndex element, double x[3], double y[3]) const;
    /** Raw template boundary markers of the three element sides. */
    std::array<std::int32_t, 3> elementEdgeMarkers(MeshIndex element) const;
    std::int32_t elementRegionAttribute(MeshIndex element) const;
    std::size_t elementInstance(MeshIndex element) const;
    ElementProvenance elementProvenance(MeshIndex element) const;
    NodeProvenance nodeProvenance(MeshIndex node) const;

    /** Global node of one instance-local template node, honouring welds. */
    MeshIndex nodeFor(std::size_t instance, MeshIndex localNode) const;
    /** Global element of one instance-local template element. */
    MeshIndex elementFor(std::size_t instance, MeshIndex localElement) const;

    /** Canonical boundary marker of an edge, or zero when the edge is unmarked. */
    std::int32_t boundaryMarkerForEdge(MeshIndex first, MeshIndex second) const;

    const std::vector<SolverMesh::PeriodicConstraint> &periodicConstraints() const
    {
        return m_periodicConstraints;
    }
    const std::vector<SolverMesh::AirGap> &airGaps() const { return m_airGaps; }

    /** Symmetric adjacency in compressed sparse-row form. */
    struct Adjacency {
        std::vector<std::size_t> offsets;
        std::vector<MeshIndex> neighbors;
        std::size_t degree(std::size_t node) const { return offsets[node + 1] - offsets[node]; }
    };

    Adjacency buildAdjacency() const;

    /**
     * Deterministic Cuthill--McKee ordering of the logical nodes. The returned
     * vector maps an old global node to its new position. The overload reuses
     * an already-built adjacency instead of rebuilding it.
     */
    std::vector<MeshIndex> cuthillMcKeeOrdering() const;
    std::vector<MeshIndex> cuthillMcKeeOrdering(const Adjacency &adjacency) const;

    /** Matrix bandwidth of an ordering: max |new(a) - new(b)| + 1 over edges. */
    static std::size_t bandwidth(const Adjacency &adjacency,
                                 const std::vector<MeshIndex> &oldToNew);
    /** Profile of an ordering: sum of (new - minimum neighbour new index). */
    static std::size_t profile(const Adjacency &adjacency,
                               const std::vector<MeshIndex> &oldToNew);

    /**
     * Approximate bytes held by the compact representation (template geometry
     * plus instance metadata). It does not grow with logical element count.
     */
    std::size_t storedByteCount() const;
    /** Approximate bytes of the equivalent fully materialised SolverMesh. */
    std::size_t expandedByteCount() const;

private:
    std::vector<MeshTemplate> m_templates;
    std::vector<MeshInstance> m_instances;
    bool m_valid = false;
    std::vector<MaterializationDiagnostic> m_diagnostics;

    std::size_t m_nodeCount = 0;
    std::size_t m_elementCount = 0;
    /** First new global node of each instance; size instanceCount + 1. */
    std::vector<std::size_t> m_instanceNodeBase;
    /** Per instance, sorted (local node, global node) pairs welded to an earlier node. */
    std::vector<std::vector<std::pair<MeshIndex, MeshIndex>>> m_weldedNodes;
    /** Sorted (global node, instance, local node) for welded global nodes. */
    std::vector<std::pair<MeshIndex, std::pair<std::size_t, MeshIndex>>> m_weldedGlobalToLocal;
    /** First global element of each instance; size instanceCount + 1. */
    std::vector<std::size_t> m_instanceElementBase;
    /** Sorted canonical edge -> nonzero boundary marker. */
    std::vector<std::pair<std::pair<MeshIndex, MeshIndex>, std::int32_t>> m_markedEdges;
    /** Per template, canonical local edge -> marker for fast element-side lookup. */
    std::vector<std::map<std::pair<MeshIndex, MeshIndex>, std::int32_t>> m_templateEdgeMarkers;
    std::vector<SolverMesh::PeriodicConstraint> m_periodicConstraints;
    std::vector<SolverMesh::AirGap> m_airGaps;
};

} // namespace mesh
} // namespace femm

#endif
