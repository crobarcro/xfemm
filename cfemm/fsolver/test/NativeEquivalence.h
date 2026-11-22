#ifndef XFEMM_TEST_NATIVE_EQUIVALENCE_H
#define XFEMM_TEST_NATIVE_EQUIVALENCE_H

#include "AnalysisSession.h"

#include <cstddef>
#include <string>
#include <vector>

namespace xfemm_test {

/** G5 test-only comparison controls. Values are in metres for coordinates and
 * webers/metre for planar magnetic vector potential, as in TrialSolution.
 * The relative scale is the maximum absolute reference A, not a per-node 1+A.
 */
struct EquivalenceOptions {
    double coordinateTolerance = 1e-7;
    double absoluteTolerance = 1e-7;
    double relativeTolerance = 1e-7;
    bool removeGaugeOffset = false;
    bool requireNontrivialField = true;
    double minimumFieldMagnitude = 1e-9;
};

struct EquivalenceResult {
    bool passed = false;
    std::string message;
    std::size_t matchedNodes = 0;
    std::size_t worstReferenceNode = 0;
    double referenceScale = 0;
    double maxAbsoluteError = 0;
    double rmsError = 0;
    /** Added to native A before comparison with materialized A. */
    double appliedGaugeOffset = 0;
};

/** Compare arbitrary solver-ordered TrialSolutions by a one-to-one geometric
 * node association. Ambiguous coincident nodes are rejected, never overwritten
 * in a coordinate map or paired by accident. Future support for truly
 * coincident independent DOFs should use an explicit provenance mapping.
 */
/** Reuse the G5.1 unambiguous coordinate association for comparing solver
 * DOF orderings and canonical mesh numbering. Mapping[trial] = reference.
 * Returns false with a diagnostic on duplicate, missing or invalid nodes.
 */
bool associateNodalCoordinates(const femm::RealNodalSolution &reference,
                               const femm::RealNodalSolution &trial,
                               double tolerance, std::vector<std::size_t> &mapping,
                               std::string &message);

EquivalenceResult compareNodalSolutions(const femm::TrialSolution &materialized,
                                         const femm::TrialSolution &native,
                                         const EquivalenceOptions &options = {});

/** Solve two independent, identically configured instanced sessions and compare.
 * The materialized session must use setNativeInstanced(false) and the other
 * setNativeInstanced(true). Validates logical/materialized topology and ensures
 * the native solve never materializes a SolverMesh.
 *
 * This intentionally uses the public AnalysisSession execution paths; no
 * direct solveNative() call may reconfigure an already-solved reference backend.
 */
EquivalenceResult compareSessionSolutions(femm::AnalysisSession &materialized,
                                           femm::AnalysisSession &native,
                                           const EquivalenceOptions &options = {});

} // namespace xfemm_test

#endif
