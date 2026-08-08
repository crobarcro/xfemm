/*
   PETSc-backed linear system backend implementation.

   The assembled system lives in PETSc sequential sparse matrices (complex
   scalar build) and is solved with a PETSc KSP.  The legacy array views
   exposed to the solvers (rhs, solution, scratch, node flags) are mirrored
   std::vector buffers that are copied to/from the PETSc vectors around the
   solve; this keeps the exact same access pattern the solvers use with the
   legacy linked-list implementation.
*/

#include "PetscLinearSystemBackend.h"

#include <petsc.h>

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <set>
#include <type_traits>
#include <vector>

namespace femm {

namespace {

/** One-time PETSc initialization for the CLI solvers (single process). */
void ensure_petsc_initialized()
{
    PetscBool initialized = PETSC_FALSE;
    PetscInitialized(&initialized);
    if (!initialized)
    {
        PetscInitializeNoArguments();
        std::atexit([]() { PetscFinalize(); });
    }
}

template <typename Scalar> PetscScalar to_petsc(const Scalar &v);
template <> PetscScalar to_petsc<double>(const double &v) { return PetscCMPLX(v, 0.0); }
template <> PetscScalar to_petsc<CComplex>(const CComplex &v)
{
    return PetscCMPLX(v.re, v.im);
}

template <typename Scalar> Scalar from_petsc(const PetscScalar &s);
template <> double from_petsc<double>(const PetscScalar &s) { return PetscRealPart(s); }
template <> CComplex from_petsc<CComplex>(const PetscScalar &s)
{
    return CComplex(PetscRealPart(s), PetscImaginaryPart(s));
}

/** Value stored at the transposed position of a matrix in matrix slot k.
 *  M/Ms are complex-symmetric, Mh is Hermitian, Ma is anti-Hermitian. */
template <typename Scalar> PetscScalar symmetric_partner(const Scalar &v, int k);
template <> PetscScalar symmetric_partner<double>(const double &v, int /*k*/)
{
    return to_petsc(v);
}
template <> PetscScalar symmetric_partner<CComplex>(const CComplex &v, int k)
{
    switch (k)
    {
    case 1: return PetscCMPLX(v.re, -v.im); // conj
    case 3: return PetscCMPLX(-v.re, v.im); // -conj
    default: return PetscCMPLX(v.re, v.im);
    }
}

} // namespace

struct PetscImpl
{
    ~PetscImpl()
    {
        if (ksp) KSPDestroy(&ksp);
        if (rhs) VecDestroy(&rhs);
        if (sol) VecDestroy(&sol);
        if (A) MatDestroy(&A);
        if (Ah) MatDestroy(&Ah);
        if (As) MatDestroy(&As);
        if (Aa) MatDestroy(&Aa);
    }

    Mat A = nullptr;   // main matrix
    Mat Ah = nullptr;  // auxiliary Mh (harmonic Newton)
    Mat As = nullptr;  // auxiliary Ms
    Mat Aa = nullptr;  // auxiliary Ma
    Vec rhs = nullptr; // RHS vector (backs rhsBuf)
    Vec sol = nullptr; // solution vector (backs solBuf)
    KSP ksp = nullptr;

    std::vector<PetscScalar> rhsBuf;
    std::vector<PetscScalar> solBuf;

    int n = 0;
    int bandwidth = 0;
    double precision = 1e-8;
    bool assembled = false;
    bool auxTouched = false;
};

namespace {

/** Preallocation estimate for a row of the sparse matrix. */
PetscInt nnz_estimate(const PetscImpl &impl)
{
    const PetscInt n = impl.n;
    if (n <= 0) return 0;
    PetscInt est = impl.bandwidth > 0 ? 2 * impl.bandwidth + 1 : n;
    if (est > n) est = n;
    return est;
}

void create_matrix(PetscImpl &impl, Mat &mat)
{
    std::vector<PetscInt> nnz(impl.n, nnz_estimate(impl));
    MatCreateSeqAIJ(PETSC_COMM_SELF, impl.n, impl.n, PETSC_DEFAULT,
                    impl.n > 0 ? nnz.data() : nullptr, &mat);
    MatSetOption(mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
}

Mat &select_matrix(PetscImpl &impl, int matrix)
{
    switch (matrix)
    {
    case 1:
        if (!impl.Ah) create_matrix(impl, impl.Ah);
        return impl.Ah;
    case 2:
        if (!impl.As) create_matrix(impl, impl.As);
        return impl.As;
    case 3:
        if (!impl.Aa) create_matrix(impl, impl.Aa);
        return impl.Aa;
    default:
        return impl.A;
    }
}

/** Flush any pending MatSetValue buffers. */
void assemble(PetscImpl &impl)
{
    if (impl.assembled) return;
    Mat mats[] = {impl.A, impl.Ah, impl.As, impl.Aa};
    for (Mat m : mats)
        if (m)
        {
            MatAssemblyBegin(m, MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(m, MAT_FINAL_ASSEMBLY);
        }
    impl.assembled = true;
}

} // namespace

// ---------------------------------------------------------------------------
// common (scalar-agnostic) interface
// ---------------------------------------------------------------------------

template <typename Scalar>
PetscLinearSystemBackend<Scalar>::PetscLinearSystemBackend()
    : m_impl(new PetscImpl)
{
    ensure_petsc_initialized();
}

template <typename Scalar>
PetscLinearSystemBackend<Scalar>::~PetscLinearSystemBackend() = default;

template <typename Scalar>
ScalarType PetscLinearSystemBackend<Scalar>::scalar_type() const
{
    return std::is_same<Scalar, double>::value ? ScalarType::Real : ScalarType::Complex;
}

template <typename Scalar>
bool PetscLinearSystemBackend<Scalar>::create(int dimension, int bandwidth,
                                              int /*node_count*/)
{
    PetscImpl &impl = *m_impl;
    impl.n = dimension;
    impl.bandwidth = bandwidth;
    impl.precision = 1e-8;
    impl.assembled = false;
    impl.auxTouched = false;

    m_rhsVec.assign(dimension, Scalar(0));
    m_solVec.assign(dimension, Scalar(0));
    m_scratchVec.assign(dimension, Scalar(0));
    m_nodeFlagVec.assign(dimension, -2);
    m_rhsView.assign(m_rhsVec.data(), static_cast<std::size_t>(dimension));
    m_solView.assign(m_solVec.data(), static_cast<std::size_t>(dimension));
    m_scratchView.assign(m_scratchVec.data(), static_cast<std::size_t>(dimension));
    m_nodeFlagView.assign(m_nodeFlagVec.data(), static_cast<std::size_t>(dimension));

    impl.rhsBuf.assign(dimension, PetscScalar(0));
    impl.solBuf.assign(dimension, PetscScalar(0));

    if (!impl.A)
        create_matrix(impl, impl.A);
    else
        MatZeroEntries(impl.A);

    if (!impl.rhs)
        VecCreateSeqWithArray(PETSC_COMM_SELF, 1, dimension,
                              impl.rhsBuf.data(), &impl.rhs);
    if (!impl.sol)
        VecCreateSeqWithArray(PETSC_COMM_SELF, 1, dimension,
                              impl.solBuf.data(), &impl.sol);

    if (!impl.ksp)
    {
        KSPCreate(PETSC_COMM_SELF, &impl.ksp);
    }
    // Bind the matrix to the KSP.  The 4-argument form is used by PETSc 3.16+;
    // earlier versions take only the two matrices.
#if PETSC_VERSION_GE(3, 16, 0)
    KSPSetOperators(impl.ksp, impl.A, impl.A, SAME_NONZERO_PATTERN);
#else
    KSPSetOperators(impl.ksp, impl.A, impl.A);
#endif

    // KSP/PC defaults.  Apply PETSc option-database overrides first, then
    // install the backend defaults for anything the user did not specify.
    // Real problems are symmetric positive definite (CG + ILU; Jacobi is too
    // weak for large machines with high permeability contrast); harmonic
    // problems are complex-symmetric non-Hermitian (BiCGStab + ILU).
    KSPSetFromOptions(impl.ksp);
    PetscBool hasKsp = PETSC_FALSE, hasPc = PETSC_FALSE;
    PetscOptionsHasName(nullptr, nullptr, "-ksp_type", &hasKsp);
    PetscOptionsHasName(nullptr, nullptr, "-pc_type", &hasPc);
    if (!hasKsp)
        KSPSetType(impl.ksp, std::is_same<Scalar, double>::value ? KSPCG : KSPBCGS);
    if (!hasPc)
    {
        PC pc;
        KSPGetPC(impl.ksp, &pc);
        PCSetType(pc, PCILU);
        // Moderate fill reduces the number of Krylov iterations on large
        // machines (ILU(0) is markedly slower); shift zero pivots so that
        // degenerate (e.g. unconstrained) matrices are handled gracefully
        // instead of aborting the factorization.
        PCFactorSetLevels(pc, 3);
        PCFactorSetShiftType(pc, MAT_SHIFT_NONZERO);
    }

    return true;
}

template <typename Scalar>
int PetscLinearSystemBackend<Scalar>::dimension() const
{
    return m_impl->n;
}

template <typename Scalar>
void PetscLinearSystemBackend<Scalar>::wipe()
{
    PetscImpl &impl = *m_impl;
    if (impl.A) MatZeroEntries(impl.A);
    if (impl.Ah) MatZeroEntries(impl.Ah);
    if (impl.As) MatZeroEntries(impl.As);
    if (impl.Aa) MatZeroEntries(impl.Aa);
    std::fill(m_rhsVec.begin(), m_rhsVec.end(), Scalar(0));
    impl.assembled = false;
}

template <typename Scalar>
void PetscLinearSystemBackend<Scalar>::put(Scalar value, int row, int col,
                                           int matrix)
{
    PetscImpl &impl = *m_impl;
    if (matrix > 0) impl.auxTouched = true;
    Mat &mat = select_matrix(impl, matrix);
    MatSetValue(mat, row, col, to_petsc<Scalar>(value), INSERT_VALUES);
    if (row != col)
        MatSetValue(mat, col, row, symmetric_partner<Scalar>(value, matrix), INSERT_VALUES);
    impl.assembled = false;
}

template <typename Scalar>
void PetscLinearSystemBackend<Scalar>::add_to(Scalar value, int row, int col,
                                              int matrix)
{
    PetscImpl &impl = *m_impl;
    if (matrix > 0) impl.auxTouched = true;
    Mat &mat = select_matrix(impl, matrix);
    MatSetValue(mat, row, col, to_petsc<Scalar>(value), ADD_VALUES);
    if (row != col)
        MatSetValue(mat, col, row, symmetric_partner<Scalar>(value, matrix), ADD_VALUES);
    impl.assembled = false;
}

template <typename Scalar>
Scalar PetscLinearSystemBackend<Scalar>::get(int row, int col, int matrix)
{
    PetscImpl &impl = *m_impl;
    if (matrix > 0 && !impl.auxTouched) return Scalar(0);
    if (matrix > 0) impl.auxTouched = true;
    Mat &mat = select_matrix(impl, matrix);
    assemble(impl);
    PetscScalar v(0);
    MatGetValue(mat, row, col, &v);
    return from_petsc<Scalar>(v);
}

template <typename Scalar>
void PetscLinearSystemBackend<Scalar>::set_value(int i, Scalar x)
{
    PetscImpl &impl = *m_impl;
    assemble(impl);

    // Capture the row, then restore before mutating the matrix.
    const PetscInt *cols = nullptr;
    const PetscScalar *vals = nullptr;
    PetscInt nnz = 0;
    MatGetRow(impl.A, i, &nnz, &cols, &vals);
    std::vector<std::pair<PetscInt, PetscScalar>> row;
    row.reserve(nnz);
    for (PetscInt j = 0; j < nnz; ++j) row.emplace_back(cols[j], vals[j]);
    MatRestoreRow(impl.A, i, &nnz, &cols, &vals);

    PetscScalar diag(0);
    for (const auto &e : row)
    {
        if (e.first == i)
        {
            diag = e.second;
            continue;
        }
        // b[k] -= A[k][i]*x  (A[k][i] == A[i][k] for the symmetric main matrix)
        m_rhsVec[e.first] -= from_petsc<Scalar>(e.second) * x;
        MatSetValue(impl.A, i, e.first, PetscScalar(0), INSERT_VALUES);
        MatSetValue(impl.A, e.first, i, PetscScalar(0), INSERT_VALUES);
    }
    m_rhsVec[i] = from_petsc<Scalar>(diag) * x;
    impl.assembled = false;
}

template <typename Scalar>
void PetscLinearSystemBackend<Scalar>::constrain_periodic(int a, int b,
                                                          bool antiperiodic)
{
    PetscImpl &impl = *m_impl;
    assemble(impl);

    const PetscInt rows[2] = {a, b};
    const PetscInt *cols[2] = {nullptr, nullptr};
    const PetscScalar *vals[2] = {nullptr, nullptr};
    PetscInt nnz[2] = {0, 0};
    std::vector<std::pair<PetscInt, PetscScalar>> row[2];
    for (int r = 0; r < 2; ++r)
    {
        MatGetRow(impl.A, rows[r], &nnz[r], &cols[r], &vals[r]);
        row[r].reserve(nnz[r]);
        for (PetscInt j = 0; j < nnz[r]; ++j) row[r].emplace_back(cols[r][j], vals[r][j]);
        MatRestoreRow(impl.A, rows[r], &nnz[r], &cols[r], &vals[r]);
    }

    // Collect every column that couples to a or b.
    std::set<PetscInt> affected;
    for (int r = 0; r < 2; ++r)
        for (const auto &e : row[r])
            if (e.first != a && e.first != b) affected.insert(e.first);

    // Average the coupling rows/columns, mirroring the legacy transform.
    for (PetscInt k : affected)
    {
        PetscScalar v1(0), v2(0);
        for (const auto &e : row[0]) if (e.first == k) { v1 = e.second; break; }
        for (const auto &e : row[1]) if (e.first == k) { v2 = e.second; break; }
        const PetscScalar c = antiperiodic ? (v1 - v2) * 0.5 : (v1 + v2) * 0.5;
        const PetscScalar c2 = antiperiodic ? -c : c;
        MatSetValue(impl.A, a, k, c, INSERT_VALUES);
        MatSetValue(impl.A, k, a, c, INSERT_VALUES);
        MatSetValue(impl.A, b, k, c2, INSERT_VALUES);
        MatSetValue(impl.A, k, b, c2, INSERT_VALUES);
    }

    // Average the diagonal entries.
    PetscScalar da(0), db(0);
    for (const auto &e : row[0]) if (e.first == a) { da = e.second; break; }
    for (const auto &e : row[1]) if (e.first == b) { db = e.second; break; }
    const PetscScalar dc = (da + db) * 0.5;
    MatSetValue(impl.A, a, a, dc, INSERT_VALUES);
    MatSetValue(impl.A, b, b, dc, INSERT_VALUES);

    // Adjust the RHS.
    const Scalar ba = m_rhsVec[a];
    const Scalar bb = m_rhsVec[b];
    if (antiperiodic)
    {
        const Scalar c = (ba - bb) * 0.5;
        m_rhsVec[a] = c;
        m_rhsVec[b] = -c;
    }
    else
    {
        const Scalar c = (ba + bb) * 0.5;
        m_rhsVec[a] = c;
        m_rhsVec[b] = c;
    }
    impl.assembled = false;
}

template <typename Scalar>
ScalarView<Scalar> &PetscLinearSystemBackend<Scalar>::rhs() { return m_rhsView; }
template <typename Scalar>
const ScalarView<Scalar> &PetscLinearSystemBackend<Scalar>::rhs() const { return m_rhsView; }
template <typename Scalar>
ScalarView<Scalar> &PetscLinearSystemBackend<Scalar>::solution() { return m_solView; }
template <typename Scalar>
const ScalarView<Scalar> &PetscLinearSystemBackend<Scalar>::solution() const { return m_solView; }

template <typename Scalar>
ScalarView<int> &PetscLinearSystemBackend<Scalar>::node_flag() { return m_nodeFlagView; }
template <typename Scalar>
const ScalarView<int> &PetscLinearSystemBackend<Scalar>::node_flag() const { return m_nodeFlagView; }
template <typename Scalar>
ScalarView<Scalar> &PetscLinearSystemBackend<Scalar>::scratch() { return m_scratchView; }

template <typename Scalar>
bool PetscLinearSystemBackend<Scalar>::newton() const { return m_impl->auxTouched; }
template <typename Scalar>
void PetscLinearSystemBackend<Scalar>::set_newton(bool b) { m_impl->auxTouched = b; }

template <typename Scalar>
double PetscLinearSystemBackend<Scalar>::precision() const { return m_impl->precision; }
template <typename Scalar>
void PetscLinearSystemBackend<Scalar>::set_precision(double p) { m_impl->precision = p; }

template <typename Scalar>
SolveReport PetscLinearSystemBackend<Scalar>::solve(const SolveOptions &options)
{
    PetscImpl &impl = *m_impl;
    SolveReport report;
    report.solver = "petsc";

    if (impl.auxTouched)
    {
        std::fprintf(stderr,
                     "xfemm PETSc backend: harmonic Newton problems (ACSolver=1) are "
                     "not yet supported; use XFEMM_SOLVER_BACKEND=legacy.\n");
        return report;
    }

    assemble(impl);

    if (options.tolerance > 0.0)
        impl.precision = options.tolerance;

    // Upload the RHS and set up the initial guess.  The initial guess is only
    // honoured when the caller requests a warm start; -ksp_type preonly
    // (direct solve) ignores the guess and rejects a nonzero-initial-guess
    // flag, so it is never advertised for preonly.
    for (int i = 0; i < impl.n; ++i) impl.rhsBuf[i] = to_petsc<Scalar>(m_rhsVec[i]);
    if (options.warm_start)
    {
        for (int i = 0; i < impl.n; ++i) impl.solBuf[i] = to_petsc<Scalar>(m_solVec[i]);
        KSPType ktype = nullptr;
        KSPGetType(impl.ksp, &ktype);
        const PetscBool useGuess =
            (ktype && std::strcmp(ktype, KSPPREONLY) != 0) ? PETSC_TRUE : PETSC_FALSE;
        KSPSetInitialGuessNonzero(impl.ksp, useGuess);
    }
    else
    {
        std::fill(impl.solBuf.begin(), impl.solBuf.end(), PetscScalar(0));
        KSPSetInitialGuessNonzero(impl.ksp, PETSC_FALSE);
    }

    // Tolerances: honour explicit PETSc options, otherwise use the caller's.
    PetscBool hasRtol = PETSC_FALSE, hasMaxit = PETSC_FALSE;
    PetscOptionsHasName(nullptr, nullptr, "-ksp_rtol", &hasRtol);
    PetscOptionsHasName(nullptr, nullptr, "-ksp_max_it", &hasMaxit);
    PetscReal curRtol = 1e-5, curAtol = 1e-50, curDtol = 1e5;
    PetscInt curMaxit = 10000;
    KSPGetTolerances(impl.ksp, &curRtol, &curAtol, &curDtol, &curMaxit);
    const PetscReal rtol = hasRtol ? curRtol
                          : (options.tolerance > 0.0 ? options.tolerance : impl.precision);
    const PetscInt maxits = hasMaxit ? curMaxit
                          : (options.max_iterations > 0 ? options.max_iterations : 1000);
    KSPSetFromOptions(impl.ksp);
    KSPSetTolerances(impl.ksp, rtol, curAtol, curDtol, maxits);

    // The matrix values change between nonlinear iterations, so force the
    // preconditioner/factorization to be rebuilt for this solve.
    KSPSetReusePreconditioner(impl.ksp, PETSC_FALSE);

    const PetscErrorCode ierr = KSPSolve(impl.ksp, impl.rhs, impl.sol);
    if (ierr != 0)
    {
        report.converged = false;
        return report;
    }

    KSPConvergedReason reason = KSP_CONVERGED_ITERATING;
    KSPGetConvergedReason(impl.ksp, &reason);
    PetscInt its = 0;
    KSPGetIterationNumber(impl.ksp, &its);
    PetscReal rnorm = 0.0;
    KSPGetResidualNorm(impl.ksp, &rnorm);
    PetscReal bnorm = 0.0;
    VecNorm(impl.rhs, NORM_2, &bnorm);

    report.converged = (reason > 0);
    report.iterations = static_cast<int>(its);
    report.relative_residual = (bnorm > 0.0) ? rnorm / bnorm : rnorm;

    // Download the solution.
    for (int i = 0; i < impl.n; ++i) m_solVec[i] = from_petsc<Scalar>(impl.solBuf[i]);

    return report;
}

template class PetscLinearSystemBackend<double>;
template class PetscLinearSystemBackend<CComplex>;

} // namespace femm
