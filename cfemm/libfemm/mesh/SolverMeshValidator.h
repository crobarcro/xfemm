#ifndef FEMM_MESH_SOLVERMESHVALIDATOR_H
#define FEMM_MESH_SOLVERMESHVALIDATOR_H

#include "SolverMesh.h"

#include <cstddef>
#include <string>
#include <vector>

namespace femm {
namespace mesh {

enum class SolverMeshValidationCategory {
    NonFiniteNodeCoordinate,
    InvalidElementNode,
    RepeatedElementNode,
    NonFiniteElementArea,
    NonPositiveElementArea,
    InvalidEdgeNode,
    SelfEdge,
    EdgeNotOwnedByElement,
    InvalidPeriodicNode,
    SelfPeriodicConstraint,
    InvalidPeriodicity,
    NonFiniteAirGapValue,
    InvalidAirGapGeometry,
    InvalidAirGapStructure,
    InvalidAirGapQuadratureNode,
    NonFiniteAirGapQuadratureWeight,
    InvalidAirGapRingNode,
    NonFiniteAirGapRingValue,
    InvalidAirGapNodeIndex
};

/** One deterministic, machine-readable failure found in a SolverMesh. */
struct SolverMeshValidationDiagnostic {
    SolverMeshValidationCategory category;
    std::size_t objectIndex = 0;
    std::size_t localIndex = 0;
    std::string context;
};

struct SolverMeshValidationResult {
    std::vector<SolverMeshValidationDiagnostic> diagnostics;
    bool valid() const { return diagnostics.empty(); }
};

/** Validate a completed SolverMesh without making assumptions about its backend. */
SolverMeshValidationResult validateSolverMesh(const SolverMesh &mesh);

} // namespace mesh
} // namespace femm

#endif
