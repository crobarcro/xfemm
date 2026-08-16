/*
   PETSc-backed linear system backend.

   This backend stores the assembled system in PETSc sequential sparse
   matrices and solves it with a PETSc KSP, giving access to the PETSc
   solver/preconditioner catalogue (CG, GMRES, BiCGStab, ... ; Jacobi,
   SSOR, ILU, Hypre AMG, MUMPS/SuperLU direct solves, ...).

   The implementation lives in PetscLinearSystemBackend.cpp; this header
   only declares the interface so that consumers never need the PETSc
   headers.  PETSc must be built with complex scalars (the Debian/Ubuntu
   "libpetsc-complex-dev" package); the real-valued backend stores real
   values with a zero imaginary part.

   Limitations:
     - The harmonic Newton-Raphson path (magnetic problems with ACSolver=1)
       is not yet supported: solve() reports failure in that case.  Use the
       legacy backend for those problems.
*/

#ifndef XFEMM_LINSOLVE_PETSC_LINEAR_SYSTEM_BACKEND_H
#define XFEMM_LINSOLVE_PETSC_LINEAR_SYSTEM_BACKEND_H

#include "LinearSystemBackend.h"

#include <memory>
#include <vector>

namespace femm {

struct PetscImpl;

/** Linear system backend backed by PETSc (complex scalar build). */
template <typename Scalar>
class PetscLinearSystemBackend : public LinearSystemBackend<Scalar>
{
public:
    PetscLinearSystemBackend();
    ~PetscLinearSystemBackend() override;

    ScalarType scalar_type() const override;
    bool create(int dimension, int bandwidth, int node_count = -1) override;
    int dimension() const override;

    void wipe() override;

    void put(Scalar value, int row, int col, int matrix = 0) override;
    void add_to(Scalar value, int row, int col, int matrix = 0) override;
    Scalar get(int row, int col, int matrix = 0) override;

    void set_value(int i, Scalar x) override;
    void constrain_periodic(int a, int b, bool antiperiodic) override;

    ScalarView<Scalar> &rhs() override;
    const ScalarView<Scalar> &rhs() const override;
    ScalarView<Scalar> &solution() override;
    const ScalarView<Scalar> &solution() const override;

    ScalarView<int> &node_flag() override;
    const ScalarView<int> &node_flag() const override;
    ScalarView<Scalar> &scratch() override;

    bool newton() const override;
    void set_newton(bool b) override;

    double precision() const override;
    void set_precision(double p) override;

    SolveReport solve(const SolveOptions &options) override;

private:
    std::unique_ptr<PetscImpl> m_impl;
    std::vector<Scalar> m_rhsVec;
    std::vector<Scalar> m_solVec;
    std::vector<Scalar> m_scratchVec;
    std::vector<int> m_nodeFlagVec;
    ScalarView<Scalar> m_rhsView;
    ScalarView<Scalar> m_solView;
    ScalarView<Scalar> m_scratchView;
    ScalarView<int> m_nodeFlagView;
};

} // namespace femm

#endif // XFEMM_LINSOLVE_PETSC_LINEAR_SYSTEM_BACKEND_H
