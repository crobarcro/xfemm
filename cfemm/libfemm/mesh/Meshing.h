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
    BackendFailure,
    ConversionFailure
};

/** Result of the complete backend-to-solver mesh operation. */
struct MeshResult {
    MeshStatus status = MeshStatus::BackendFailure;
    SolverMesh mesh;
    std::vector<MeshDiagnostic> diagnostics;

    bool succeeded() const
    {
        return status == MeshStatus::Success || status == MeshStatus::SuccessWithWarnings;
    }
};

} // namespace mesh
} // namespace femm

#endif
