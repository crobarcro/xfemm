#ifndef FEMM_FMESHER_BOUNDARYMATCHVALIDATOR_H
#define FEMM_FMESHER_BOUNDARYMATCHVALIDATOR_H

#include "mesh/Meshing.h"

#include <vector>

namespace femm { class FemmProblem; }

namespace fmesher {

struct BoundaryMatchValidation {
    std::vector<femm::mesh::MeshDiagnostic> diagnostics;
    bool valid() const { return diagnostics.empty(); }
};

/**
 * Validate a request's boundary matches against a problem before any backend
 * runs. Rejects stale (out-of-range), mixed line/arc, duplicated, and missing
 * references, and references to boundaries the problem does not declare as a
 * pair of periodic entities.
 */
BoundaryMatchValidation validateBoundaryMatches(
    const femm::FemmProblem &problem,
    const femm::mesh::MeshingRequest &request);

} // namespace fmesher

#endif
