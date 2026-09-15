#ifndef FEMM_MESH_MESHING_H
#define FEMM_MESH_MESHING_H

#include "SolverMesh.h"

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
    std::vector<MeshDiagnostic> diagnostics;

    bool succeeded() const
    {
        return status == MeshStatus::Success || status == MeshStatus::SuccessWithWarnings;
    }
};

} // namespace mesh
} // namespace femm

#endif
