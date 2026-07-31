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

    /** Write the most recently solved field in the legacy .ans format. */
    void writeSolution(const std::string &ansPath);

    std::size_t topologyImportCount() const { return m_topologyImports; }
    std::size_t orderingCount() const { return m_orderings; }
    std::size_t couplingRegenerationCount() const { return m_couplingRegenerations; }
    std::size_t operatorAssemblyCount() const { return m_operatorAssemblies; }
    std::size_t rightHandSideAssemblyCount() const { return m_rightHandSideAssemblies; }
    std::size_t solveCount() const { return m_solves; }
    std::size_t meshFileReadCount() const { return m_meshFileReads; }
    std::size_t meshFileWriteCount() const { return m_meshFileWrites; }

private:
    void configure(const ModelDefinition &, const SolveParameters &,
                   const PreparedAnalysis &);
    void positionAirGaps(const PreparedAnalysis &);
    std::unique_ptr<FSolver> m_solver;
    std::unique_ptr<CBigLinProb> m_lastSystem;
    std::shared_ptr<const mesh::SolverMesh> m_mesh;
    std::uint64_t m_topologyIdentity = 0;
    std::size_t m_topologyImports = 0;
    std::size_t m_orderings = 0;
    std::size_t m_couplingRegenerations = 0;
    std::size_t m_operatorAssemblies = 0;
    std::size_t m_rightHandSideAssemblies = 0;
    std::size_t m_solves = 0;
    // The session path is in-memory. These counters make accidental legacy
    // mesh-file I/O observable to callers and regression tests.
    std::size_t m_meshFileReads = 0;
    std::size_t m_meshFileWrites = 0;
};

} // namespace femm

#endif
