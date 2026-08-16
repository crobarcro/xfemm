/*
   Legacy linear system backend: wraps the original linked-list solvers
   (CBigLinProb / CBigComplexLinProb) behind the LinearSystemBackend
   interface.

   This backend preserves the exact numerical behaviour of the original
   solver, which makes it the validation oracle for any alternative
   backend (PETSc, Eigen, ...).
*/

#ifndef XFEMM_LINSOLVE_LEGACY_LINEAR_SYSTEM_BACKEND_H
#define XFEMM_LINSOLVE_LEGACY_LINEAR_SYSTEM_BACKEND_H

#include "LinearSystemBackend.h"
#include "cspars.h"
#include "spars.h"

namespace femm {

template <typename Scalar>
class LegacyLinearSystemBackend;

/** Real-valued legacy backend (CBigLinProb). */
template <>
class LegacyLinearSystemBackend<double> : public LinearSystemBackend<double>
{
public:
    ScalarType scalar_type() const override { return ScalarType::Real; }

    bool create(int dimension, int bandwidth, int /*node_count*/) override
    {
        if (m_system.Create(dimension, bandwidth) == 0)
            return false;
        const int n = m_system.n;
        m_rhs.assign(m_system.b, static_cast<std::size_t>(n));
        m_solution.assign(m_system.V, static_cast<std::size_t>(n));
        m_scratch.assign(m_system.P, static_cast<std::size_t>(n));
        m_nodeFlag.assign(m_system.Q, static_cast<std::size_t>(n));
        return true;
    }

    int dimension() const override { return m_system.n; }

    void wipe() override { m_system.Wipe(); }

    void put(double value, int row, int col, int /*matrix*/) override
    {
        m_system.Put(value, row, col);
    }
    void add_to(double value, int row, int col, int /*matrix*/) override
    {
        m_system.AddTo(value, row, col);
    }
    double get(int row, int col, int /*matrix*/) override
    {
        return m_system.Get(row, col);
    }

    void set_value(int i, double x) override { m_system.SetValue(i, x); }
    void constrain_periodic(int a, int b, bool antiperiodic) override
    {
        if (antiperiodic)
            m_system.AntiPeriodicity(a, b);
        else
            m_system.Periodicity(a, b);
    }

    ScalarView<double> &rhs() override { return m_rhs; }
    const ScalarView<double> &rhs() const override { return m_rhs; }
    ScalarView<double> &solution() override { return m_solution; }
    const ScalarView<double> &solution() const override { return m_solution; }

    ScalarView<int> &node_flag() override { return m_nodeFlag; }
    const ScalarView<int> &node_flag() const override { return m_nodeFlag; }
    ScalarView<double> &scratch() override { return m_scratch; }

    bool newton() const override { return false; }
    void set_newton(bool) override {}

    double precision() const override { return m_system.Precision; }
    void set_precision(double p) override { m_system.Precision = p; }

    SolveReport solve(const SolveOptions &options) override
    {
        if (options.tolerance > 0.0)
            m_system.Precision = options.tolerance;
        SolveReport report;
        report.solver = "legacy-pcg";
        report.converged = m_system.PCGSolve(options.warm_start ? 1 : 0);
        return report;
    }

private:
    CBigLinProb m_system;
    ScalarView<double> m_rhs;
    ScalarView<double> m_solution;
    ScalarView<double> m_scratch;
    ScalarView<int> m_nodeFlag;
};

/** Complex-valued legacy backend (CBigComplexLinProb). */
template <>
class LegacyLinearSystemBackend<CComplex> : public LinearSystemBackend<CComplex>
{
public:
    ScalarType scalar_type() const override { return ScalarType::Complex; }

    bool create(int dimension, int bandwidth, int node_count) override
    {
        const int nodes = (node_count < 0) ? dimension : node_count;
        if (m_system.Create(dimension, bandwidth, nodes) == 0)
            return false;
        const int n = m_system.n;
        m_rhs.assign(m_system.b, static_cast<std::size_t>(n));
        m_solution.assign(m_system.V, static_cast<std::size_t>(n));
        m_scratch.assign(m_system.P, static_cast<std::size_t>(n));
        m_nodeFlag.assign(nullptr, 0);
        return true;
    }

    int dimension() const override { return m_system.n; }

    void wipe() override { m_system.Wipe(); }

    void put(CComplex value, int row, int col, int matrix) override
    {
        m_system.Put(value, row, col, matrix);
    }
    void add_to(CComplex value, int row, int col, int /*matrix*/) override
    {
        m_system.AddTo(value, row, col);
    }
    CComplex get(int row, int col, int matrix) override
    {
        return m_system.Get(row, col, matrix);
    }

    void set_value(int i, CComplex x) override { m_system.SetValue(i, x); }
    void constrain_periodic(int a, int b, bool antiperiodic) override
    {
        if (antiperiodic)
            m_system.AntiPeriodicity(a, b);
        else
            m_system.Periodicity(a, b);
    }

    ScalarView<CComplex> &rhs() override { return m_rhs; }
    const ScalarView<CComplex> &rhs() const override { return m_rhs; }
    ScalarView<CComplex> &solution() override { return m_solution; }
    const ScalarView<CComplex> &solution() const override { return m_solution; }

    ScalarView<int> &node_flag() override { return m_nodeFlag; }
    const ScalarView<int> &node_flag() const override { return m_nodeFlag; }
    ScalarView<CComplex> &scratch() override { return m_scratch; }

    bool newton() const override { return m_system.bNewton; }
    void set_newton(bool b) override { m_system.bNewton = b; }

    double precision() const override { return m_system.Precision; }
    void set_precision(double p) override { m_system.Precision = p; }

    SolveReport solve(const SolveOptions &options) override
    {
        if (options.tolerance > 0.0)
            m_system.Precision = options.tolerance;
        SolveReport report;
        report.solver = "legacy-pbcg";
        report.converged =
            m_system.PBCGSolveMod(options.warm_start ? 1 : 0, options.verbose) != 0;
        return report;
    }

private:
    CBigComplexLinProb m_system;
    ScalarView<CComplex> m_rhs;
    ScalarView<CComplex> m_solution;
    ScalarView<CComplex> m_scratch;
    ScalarView<int> m_nodeFlag;
};

} // namespace femm

#endif // XFEMM_LINSOLVE_LEGACY_LINEAR_SYSTEM_BACKEND_H
