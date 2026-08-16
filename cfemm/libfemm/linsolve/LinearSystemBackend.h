/*
   Abstraction layer for linear system assembly and solution.

   This file introduces the LinearSystemBackend interface that the
   solvers (fsolver, esolver, hsolver) use instead of reaching directly
   for the legacy CBigLinProb / CBigComplexLinProb classes.  The intended
   implementations are:

     - LegacyLinearSystemBackend (the current linked-list solver, kept as
       the validation oracle),
     - a future PETSc backend (and possibly an Eigen backend).

   Element assembly and FEMM conventions stay in the solvers; only the
   linear algebra storage and solve primitives move behind this interface.
*/

#ifndef XFEMM_LINSOLVE_LINEAR_SYSTEM_BACKEND_H
#define XFEMM_LINSOLVE_LINEAR_SYSTEM_BACKEND_H

#include "femmcomplex.h"

#include <cstddef>
#include <string>

namespace femm {

enum class ScalarType { Real, Complex };

/** Options controlling a single solve. */
struct SolveOptions {
    /** Requested relative tolerance; values <= 0 leave the backend's
     *  stored precision untouched. */
    double tolerance = -1.0;
    /** Maximum number of iterations; 0 means "use the backend default". */
    int max_iterations = 0;
    /** Reuse the current solution() as an initial guess (nonlinear
     *  / incremental iterations). */
    bool warm_start = false;
    bool verbose = false;
};

/** Outcome of a solve call. */
struct SolveReport {
    bool converged = false;
    int iterations = 0;
    double relative_residual = -1.0;
    std::string solver;
};

/** Non-owning view over a contiguous array of scalars.
 *
 *  The legacy solvers allocate their work arrays with raw new[]/calloc and
 *  expose them as bare pointers.  The backend implementations hand those
 *  arrays out through this lightweight span-like view so callers can keep
 *  using subscript syntax without owning the storage.
 */
template <typename T>
class ScalarView
{
public:
    ScalarView() = default;
    ScalarView(T *data, std::size_t size)
        : m_data(data), m_size(size) {}

    void assign(T *data, std::size_t size)
    {
        m_data = data;
        m_size = size;
    }

    std::size_t size() const { return m_size; }
    bool empty() const { return m_size == 0; }

    T &operator[](std::size_t index) { return m_data[index]; }
    const T &operator[](std::size_t index) const { return m_data[index]; }

    T *data() { return m_data; }
    const T *data() const { return m_data; }

    T *begin() { return m_data; }
    T *end() { return m_data + m_size; }
    const T *begin() const { return m_data; }
    const T *end() const { return m_data + m_size; }

private:
    T *m_data = nullptr;
    std::size_t m_size = 0;
};

/** Abstract interface for a sparse linear system backend.
 *
 *  The life cycle is:
 *   - create(dimension, bandwidth) once;
 *   - wipe() + assemble entries (put/add_to) + impose boundary conditions
 *     (set_value, constrain_periodic);
 *   - solve() one or more times;
 *   - read the results from solution()/rhs().
 *
 *  \tparam Scalar the value type: double for real problems, CComplex for
 *          harmonic (frequency-domain) problems.
 */
template <typename Scalar>
class LinearSystemBackend
{
public:
    virtual ~LinearSystemBackend() = default;

    virtual ScalarType scalar_type() const = 0;

    /** Allocate storage.  bandwidth is a speed hint for the legacy backend
     *  only; node_count is the number of field (non-circuit) unknowns, which
     *  the legacy complex backend needs to bound certain transformations.
     *  Returns false if the allocation failed. */
    virtual bool create(int dimension, int bandwidth, int node_count = -1) = 0;
    virtual int dimension() const = 0;

    /** Reset the matrix values and the right-hand side to zero, keeping the
     *  allocated storage and the current solution() intact. */
    virtual void wipe() = 0;

    /** Matrix assembly.  matrix selects the matrix to operate on:
     *  0 = M, 1 = Mh, 2 = Ms, 3 = Ma (the last three are the auxiliary
     *  Newton-Raphson matrices used by the harmonic solvers; real backends
     *  ignore the selector). */
    virtual void put(Scalar value, int row, int col, int matrix = 0) = 0;
    virtual void add_to(Scalar value, int row, int col, int matrix = 0) = 0;
    virtual Scalar get(int row, int col, int matrix = 0) = 0;

    /** Impose a Dirichlet value at node i, eliminating the corresponding
     *  row and column (and adjusting the RHS). */
    virtual void set_value(int i, Scalar x) = 0;

    /** Constrain nodes a and b to share a variable: equal for periodic
     *  boundaries, opposite for antiperiodic boundaries. */
    virtual void constrain_periodic(int a, int b, bool antiperiodic) = 0;

    /** Right-hand side and solution vectors.  Backends may reuse these
     *  vectors between solves. */
    virtual ScalarView<Scalar> &rhs() = 0;
    virtual const ScalarView<Scalar> &rhs() const = 0;
    virtual ScalarView<Scalar> &solution() = 0;
    virtual const ScalarView<Scalar> &solution() const = 0;

    /** Legacy bookkeeping: an integer flag per unknown (used by the
     *  electrostatic and heat-flow solvers) and a scratch vector (reused by
     *  the conductor-charge computations). */
    virtual ScalarView<int> &node_flag() = 0;
    virtual const ScalarView<int> &node_flag() const = 0;
    virtual ScalarView<Scalar> &scratch() = 0;

    /** Newton-Raphson flag for the harmonic solvers (true once the
     *  auxiliary matrices have been populated). */
    virtual bool newton() const = 0;
    virtual void set_newton(bool b) = 0;

    /** Stored solve precision. */
    virtual double precision() const = 0;
    virtual void set_precision(double p) = 0;

    /** Solve A x = b. */
    virtual SolveReport solve(const SolveOptions &options) = 0;
};

} // namespace femm

#endif // XFEMM_LINSOLVE_LINEAR_SYSTEM_BACKEND_H
