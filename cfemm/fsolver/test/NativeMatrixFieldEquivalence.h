#ifndef XFEMM_TEST_NATIVE_MATRIX_FIELD_EQUIVALENCE_H
#define XFEMM_TEST_NATIVE_MATRIX_FIELD_EQUIVALENCE_H

#include "FSolverAnalysisBackend.h"
#include "NativeEquivalence.h"

#include <cstddef>
#include <string>
#include <vector>

namespace xfemm_test {

/** G5.2 tolerances are independent of those for nodal A. Linear-system
 * quantities use the internal FEMM solver scaling, while B/H use SI units.
 */
struct MatrixFieldOptions {
    EquivalenceOptions nodal;
    double matrixAbsoluteTolerance = 1e-8;
    double matrixRelativeTolerance = 1e-7;
    double rhsAbsoluteTolerance = 1e-8;
    double rhsRelativeTolerance = 1e-7;
    /** Filter numerical zeros before comparing the *numerical* nonzero pattern.
     * Exact storage allocation patterns may include zero-valued entries.
     */
    double sparsityThreshold = 1e-13;
    double maximumRelativeResidual = 1e-5;
    /** Difference between reference and native normalized equation residuals,
     * even when both satisfy an AGE-specific absolute residual cap.
     */
    double residualAgreementTolerance = 1e-5;
    double fieldAbsoluteTolerance = 1e-6; // tesla
    double fieldRelativeTolerance = 1e-5;
    double hAbsoluteTolerance = 1e-3;     // ampere/metre
    double hRelativeTolerance = 1e-5;
    bool compareMatricesAndRhs = true;
    bool compareFields = true;
    bool compareH = true;
};

struct MatrixComparisonResult {
    bool passed = false;
    std::string message;
    std::size_t numericNonzeros = 0;
    std::size_t sparsityMismatches = 0;
    std::size_t coefficientMismatches = 0;
    double maximumMatrixDifference = 0;
    double maximumRelativeMatrixDifference = 0;
    double maximumRhsDifference = 0;
    double referenceRelativeResidual = 0;
    double nativeRelativeResidual = 0;
};

/** Compare snapshots in different DOF orderings, including the actual
 * pre-solve RHS and the equation residual K*x-b. The mapping is native
 * solver DOF -> materialized solver DOF (not logical global-node index).
 */
MatrixComparisonResult compareAssemblySnapshots(
    const femm::FSolverAnalysisBackend::StaticAssemblySnapshot &reference,
    const femm::FSolverAnalysisBackend::StaticAssemblySnapshot &native,
    const std::vector<std::size_t> &nativeToReference,
    const std::vector<double> &referenceSolution,
    const std::vector<double> &nativeSolution,
    const MatrixFieldOptions &options = {});

struct MatrixFieldResult {
    bool passed = false;
    std::string message;
    EquivalenceResult nodal;
    MatrixComparisonResult matrix;
    std::size_t comparedElements = 0;
    double maximumBError = 0;
    double maximumHError = 0;
    double referenceBScale = 0;
    double referenceHScale = 0;
};

/** Execute both independent sessions, compare nodal A and the final assembled
 * system, then compute constant-per-triangle B and constitutive H from the
 * canonical materialized mesh using each solution.
 *
 * Requires matching instancing and static planar physics. Matrix snapshots
 * are opt-in and never used by production solvers. Only materialized tests
 * allocate global elements; a native solve must remain unmaterialized.
 */
MatrixFieldResult compareMatrixFieldSessions(
    femm::AnalysisSession &referenceSession,
    femm::FSolverAnalysisBackend &referenceSolver,
    femm::AnalysisSession &nativeSession,
    femm::FSolverAnalysisBackend &nativeSolver,
    const MatrixFieldOptions &options = {});

} // namespace xfemm_test

#endif
