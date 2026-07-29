#ifndef XFEMM_FSOLVERANALYSISBACKEND_H
#define XFEMM_FSOLVERANALYSISBACKEND_H

#include "AnalysisSession.h"
#include "fsolver.h"

#include <cstddef>
#include <memory>

namespace femm {

/** AnalysisSession adapter for the legacy magnetic FSolver engine.
 *
 * The imported mesh and its Cuthill--McKee numbering are owned by this
 * backend and replaced only when meshTopologyIdentity changes.  Linear-system
 * factorization is deliberately not cached yet.
 */
class FSolverAnalysisBackend final : public AnalysisSolverBackend {
public:
    FSolverAnalysisBackend();
    ~FSolverAnalysisBackend() override;

    void synchronize(const ModelDefinition &, const SolveParameters &,
                     const PreparedAnalysis &, std::shared_ptr<const mesh::SolverMesh>,
                     std::uint64_t meshTopologyIdentity, Dirty rebuilt) override;
    TrialSolution solve(const ModelDefinition &, const SolveParameters &,
                        const PreparedAnalysis &) override;

    std::size_t topologyImportCount() const { return m_topologyImports; }
    std::size_t orderingCount() const { return m_orderings; }
    std::size_t couplingRegenerationCount() const { return m_couplingRegenerations; }

private:
    void configure(const ModelDefinition &, const SolveParameters &,
                   const PreparedAnalysis &);
    void positionAirGaps(const PreparedAnalysis &);
    std::unique_ptr<FSolver> m_solver;
    std::shared_ptr<const mesh::SolverMesh> m_mesh;
    std::uint64_t m_topologyIdentity = 0;
    std::size_t m_topologyImports = 0;
    std::size_t m_orderings = 0;
    std::size_t m_couplingRegenerations = 0;
};

} // namespace femm

#endif
