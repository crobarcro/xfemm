#ifndef FEMM_MESH_INSTANCEDMESH_H
#define FEMM_MESH_INSTANCEDMESH_H

#include "SolverMesh.h"

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

namespace femm {
namespace mesh {

/**
 * Orientation-preserving rigid transform in the xy plane.
 *
 * The type deliberately cannot express a reflection: reflections reverse
 * element orientation and magnetisation handedness and need separate semantics.
 * A rotation followed by a translation preserves signed element area.
 */
struct RigidTransform2D {
    double rotationDegrees = 0.0;
    double translationXMetres = 0.0;
    double translationYMetres = 0.0;

    bool isFinite() const;
    /** Transform a point: rotate about the origin, then translate. */
    void applyPoint(double x, double y, double &outX, double &outY) const;
    /** Transform a direction: rotate only. */
    void applyDirection(double x, double y, double &outX, double &outY) const;
    /** Build a transform that rotates about a declared centre. */
    static RigidTransform2D rotationAbout(double centerXMetres, double centerYMetres,
                                          double angleDegrees);
};

/**
 * A named ordered chain of local nodes forming one side of a template seam.
 * The order is the correspondence used for welding, not a geometric search.
 */
struct TemplateSeam {
    std::string name;
    std::vector<MeshIndex> orderedNodes;
};

/** A template mesh and the seams that may be welded to other instances. */
struct MeshTemplate {
    SolverMesh localMesh;
    std::vector<TemplateSeam> seams;
};

/** How a connection aligns the other seam's ordered nodes with this one. */
enum class SeamOrientation { Forward, Reverse };

/**
 * Declare that one seam of this instance is welded to a seam of another
 * instance. Connections are undirected: declaring the same pair from either
 * side is equivalent, and a conflicting orientation is rejected.
 */
struct SeamConnection {
    std::size_t seam = 0;
    std::size_t otherInstance = 0;
    std::size_t otherSeam = 0;
    SeamOrientation orientation = SeamOrientation::Forward;
};

/**
 * Per-instance physics override for one source block label of the template.
 *
 * Physics metadata deliberately lives on the instance, never on the shared
 * template. A rigid transform rotates directional quantities; it must not
 * mutate the template. The initial implementation supports isotropic
 * materials only.
 */
struct InstanceRegionOverride {
    /** Zero-based index into the template problem's block-label list. */
    std::size_t sourceBlockLabel = 0;
    /** Replace the label's circuit membership. */
    std::optional<std::size_t> circuit;
    /** Add to the label's magnetisation direction, in degrees. */
    std::optional<double> magnetisationRotationDegrees;
    /** Scale the label's signed turns (a negative value reverses current). */
    std::optional<double> currentScale;
};

/** One placed occurrence of a template. */
struct MeshInstance {
    std::size_t templateIndex = 0;
    RigidTransform2D transform;
    std::vector<SeamConnection> seamConnections;
    /** Per-instance physics; not part of the topology or layout identity. */
    std::vector<InstanceRegionOverride> regionOverrides;
};

/** Stable diagnostic category for a rejected instanced mesh. */
enum class MaterializationDiagnosticCategory {
    InvalidTemplateIndex,
    NonFiniteTransform,
    InvalidInstanceSeam,
    InvalidConnectionInstance,
    InvalidConnectionSeam,
    SelfConnection,
    ConflictingConnectionOrientation,
    EmptySeam,
    InvalidSeamNode,
    DuplicateSeamNode,
    DisconnectedSeamChain,
    SeamCardinalityMismatch,
    SeamCoordinateMismatch,
    DegenerateElement,
    ReversedElement,
    InvalidElementNode,
    InvalidEdgeNode,
    EdgeMarkerConflict,
    InvalidPeriodicNode,
    InvalidAirGapNode,
    InvalidAirGapStructure,
    IndexOverflow
};

/** One deterministic, machine-readable reason a materialisation was rejected. */
struct MaterializationDiagnostic {
    MaterializationDiagnosticCategory category;
    std::size_t instanceIndex = 0;
    std::size_t objectIndex = 0;
    std::size_t localIndex = 0;
    std::string message;
};

/** Provenance of a materialised global node or element. */
struct NodeProvenance {
    std::size_t templateIndex = 0;
    std::size_t instanceIndex = 0;
    MeshIndex localNode = InvalidMeshIndex;
};

struct ElementProvenance {
    std::size_t templateIndex = 0;
    std::size_t instanceIndex = 0;
    MeshIndex localElement = InvalidMeshIndex;
};

/** Local-to-global maps exposed alongside a materialised mesh. */
struct InstancingProvenance {
    /** [instance][local node] -> global node. */
    std::vector<std::vector<MeshIndex>> nodeMap;
    /** [instance][local element] -> global element. */
    std::vector<std::vector<MeshIndex>> elementMap;
    /** [global node] -> template/instance/local node (a weld representative). */
    std::vector<NodeProvenance> nodeProvenance;
    /** [global element] -> template/instance/local element. */
    std::vector<ElementProvenance> elementProvenance;
};

/**
 * A materialised SolverMesh plus the maps needed to translate between local and
 * global indices. The mesh is only meaningful when succeeded() is true.
 */
struct MaterializationResult {
    SolverMesh mesh;
    InstancingProvenance provenance;
    std::vector<MaterializationDiagnostic> diagnostics;

    bool succeeded() const { return diagnostics.empty(); }
};

/**
 * A set of templates and their placed instances.
 *
 * Materialisation expands each instance with its transform, welds only the
 * explicitly connected seam nodes, and remaps periodic and AGE topology. It
 * never identifies fields across a seam: that remains an explicit solver
 * constraint.
 */
struct InstancedMesh {
    std::vector<MeshTemplate> templates;
    std::vector<MeshInstance> instances;

    /** Validate template-local topology, transforms, seams, and connections. */
    std::vector<MaterializationDiagnostic> validate() const;

    /** Deterministically expand the instances into one SolverMesh. */
    MaterializationResult materialize() const;
};

/** Metre-scale tolerance used only to reject an inconsistent seam weld. */
constexpr double InstancedMeshWeldToleranceMetres = 1e-9;

/**
 * Stable 64-bit identities for cache invalidation. They depend only on value
 * data (never pointers or addresses) and are identical across repeated runs.
 * templateTopologyIdentity covers template meshes and seams;
 * instanceLayoutIdentity covers instance transforms and seam connections;
 * materializedTopologyIdentity covers a completed SolverMesh.
 */
std::uint64_t templateTopologyIdentity(const InstancedMesh &instanced);
std::uint64_t instanceLayoutIdentity(const InstancedMesh &instanced);
std::uint64_t materializedTopologyIdentity(const SolverMesh &mesh);

} // namespace mesh
} // namespace femm

#endif
