#ifndef FEMM_MESH_MESHING_H
#define FEMM_MESH_MESHING_H

#include "InstancedMesh.h"
#include "SolverMesh.h"

#include <optional>
#include <string>
#include <vector>

namespace femm {
namespace mesh {

/** Backend-neutral controls for creating a mesh. */
struct MeshingOptions {
    double minimumAngleDegrees = 0.0; ///< Zero requests the backend default.
    double defaultElementSize = 0.0; ///< Source-problem units; zero requests automatic sizing.
    bool forceMaximumElementArea = false;
    bool suppressExteriorSteinerPoints = false;
    bool suppressUnusedVertices = false;
    bool verbose = false;
};

enum class MeshDiagnosticSeverity { Information, Warning, Error };

/**
 * A diagnostic safe to display or log without interpreting backend errors.
 * backendErrorCode is opaque to the solver; only the named backend defines it.
 */
struct MeshDiagnostic {
    MeshDiagnosticSeverity severity = MeshDiagnosticSeverity::Information;
    std::string message;
    std::string backendName;
    int backendErrorCode = 0;
};

enum class MeshStatus {
    Success,
    SuccessWithWarnings,
    InvalidInput,
    /** The selected backend cannot perform part of the request. */
    Unsupported,
    BackendFailure,
    ConversionFailure
};

/**
 * Stable reference to one source boundary entity in a FemmProblem.
 *
 * Boundary matching deliberately references source entities by kind and index
 * rather than by transient mesh indices: adding or removing unrelated geometry
 * keeps a reference valid, while a removed or reordered entity is detected as
 * out of range before the backend runs.
 */
struct SourceEntityReference {
    enum class Kind { Segment, Arc };
    Kind kind = Kind::Segment;
    /** Zero-based index into FemmProblem::linelist or arclist. */
    std::size_t index = 0;
};

/**
 * Ask the mesher to give two source boundary entities the same ordered
 * discretisation. Matching topology is not a physical boundary condition, so a
 * match only emits a periodic or antiperiodic field constraint when its own
 * createFieldConstraint flag says so.
 *
 * The two references must name the same periodic boundary property (Tangle
 * pairs FEMM boundaries by property), must be the same geometry kind, and that
 * property must name exactly these two entities. The request's periodicity must
 * match the property's declared periodic or antiperiodic type.
 */
struct BoundaryMatch {
    SourceEntityReference first;
    SourceEntityReference second;
    bool createFieldConstraint = false;
    SolverMesh::Periodicity periodicity = SolverMesh::Periodicity::Periodic;
};

/**
 * Mesh one rotational tile once and repeat it around a centre.
 *
 * The initial implementation treats the loaded problem as the tile: its
 * matched periodic boundaries are used as topology-only seams, so no temporary
 * problem file or upstream in-memory PSLG API is required. Selecting a tile
 * out of a larger model is deferred until that API exists.
 */
struct TemplateRequest {
    /** Centre of rotation in metres. */
    double centerXMetres = 0.0;
    double centerYMetres = 0.0;
    /** Number of rotational instances; must be at least two. */
    std::size_t instanceCount = 1;
    /** Total angle covered by all instances; a closed ring requires 360. */
    double totalAngleDegrees = 360.0;
    /**
     * Boundary properties whose matched chains become the template seams.
     * Empty selects every match the mesher reports (which must be exactly one
     * pair for a rotational tile).
     */
    std::vector<std::size_t> seamBoundaryProperties;
};

/** Backend-neutral description of one meshing operation. */
struct MeshingRequest {
    MeshingOptions options;
    /**
     * Compatibility flag for the pre-request API, used only when
     * boundaryMatches is empty: emit ordinary periodic or antiperiodic
     * constraints for every boundary the problem declares. When
     * boundaryMatches is non-empty it selects exactly which matches emit
     * constraints.
     */
    bool createPeriodicFieldConstraints = false;
    std::vector<BoundaryMatch> boundaryMatches;
    /** At most one rotational template is supported initially. */
    std::vector<TemplateRequest> templates;
};

/** Direction of the second matched chain relative to the first. */
enum class BoundaryMatchOrientation { Forward, Reverse };

/**
 * Ordered node correspondence produced for one matched boundary pair.
 * firstNodes[i] and secondNodes[i] are one-to-one; the correspondence is
 * topology only and is never itself a periodic field constraint.
 */
struct MeshBoundaryMatch {
    /** Zero-based index into FemmProblem::lineproplist of the shared boundary. */
    std::size_t boundaryProperty = 0;
    SolverMesh::Periodicity periodicity = SolverMesh::Periodicity::Periodic;
    BoundaryMatchOrientation orientation = BoundaryMatchOrientation::Forward;
    std::vector<MeshIndex> firstNodes;
    std::vector<MeshIndex> secondNodes;
};

/** Result of the complete backend-to-solver mesh operation. */
struct MeshResult {
    MeshStatus status = MeshStatus::BackendFailure;
    SolverMesh mesh;
    /** Ordered seam correspondences; additive and independent of SolverMesh. */
    std::vector<MeshBoundaryMatch> boundaryMatches;
    /**
     * Canonical instanced templates when the mesh was produced by materialising
     * a template, before materialisation. Additive: a session can cache this
     * and re-materialise on a layout change without re-meshing the tile.
     */
    std::optional<InstancedMesh> instancedTemplates;
    /**
     * Local-to-global provenance when the mesh was produced by materialising a
     * template. Additive: consumers that do not instance ignore it.
     */
    std::optional<InstancingProvenance> instancing;
    std::vector<MeshDiagnostic> diagnostics;

    bool succeeded() const
    {
        return status == MeshStatus::Success || status == MeshStatus::SuccessWithWarnings;
    }
};

} // namespace mesh
} // namespace femm

#endif
