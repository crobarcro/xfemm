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
    void synchronizeNative(const ModelDefinition &, const SolveParameters &,
                           const PreparedAnalysis &,
                           std::shared_ptr<const mesh::LogicalMeshView>,
                           const std::vector<std::size_t> &,
                           const std::map<std::string, std::pair<double, double>> &,
                           Dirty rebuilt) override;
    TrialSolution solve(const ModelDefinition &, const SolveParameters &,
                        const PreparedAnalysis &) override;

    /** Write the most recently solved field in the legacy .ans format. */
    void writeSolution(const std::string &ansPath);

    /**
     * Solve the current prepared problem through the native compressed path,
     * iterating a logical view instead of the imported materialised mesh.
     * Returns false when the native path does not support the configured
     * problem (nonlinear materials, functional magnetisation, non-planar
     * coordinates, incremental previous solutions, or an unsupported AGE).
     * On success nativeSystem() holds the result, indexed by logical global
     * node rather than by the materialised Cuthill ordering.
     */
    bool solveNative(const mesh::LogicalMeshView &view,
                     const std::vector<std::size_t> &instanceLabelBase,
                     const std::map<std::string, std::pair<double, double>>
                         &airGapPositions = {},
                     const std::vector<std::size_t> &nodePermutation = {},
                     int bandwidth = -1);
    const femm::LinearSystemBackend<double> &nativeSystem() const;
    std::size_t nativeSolveCount() const { return m_nativeSolves; }
    /** System bandwidth of the last native setup. */
    int nativeBandwidth() const { return m_nativeBandwidth; }

    /** Native solved state used to construct an in-memory post-processor view. */
    const FSolver &solvedSolver() const;
    const femm::LinearSystemBackend<double> &solvedSystem() const;

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
    std::unique_ptr<femm::LinearSystemBackend<double>> m_lastSystem;
    std::shared_ptr<const mesh::SolverMesh> m_mesh;
    bool m_nativeMode = false;
    std::shared_ptr<const mesh::LogicalMeshView> m_nativeView;
    std::vector<std::size_t> m_nativeLabelBases;
    std::map<std::string, std::pair<double, double>> m_nativeAirGapPositions;
    /** Cuthill-McKee ordering (old -> new) cached with the native view. */
    std::vector<std::size_t> m_nativeOrdering;
    int m_nativeBandwidth = 0;
    /** new -> old global node map of the last native assembly. */
    std::vector<std::size_t> m_nativeNewToOld;
    std::uint64_t m_topologyIdentity = 0;
    std::size_t m_topologyImports = 0;
    std::size_t m_orderings = 0;
    std::size_t m_couplingRegenerations = 0;
    std::size_t m_operatorAssemblies = 0;
    std::size_t m_rightHandSideAssemblies = 0;
    std::size_t m_solves = 0;
    std::size_t m_nativeSolves = 0;
    // The session path is in-memory. These counters make accidental legacy
    // mesh-file I/O observable to callers and regression tests.
    std::size_t m_meshFileReads = 0;
    std::size_t m_meshFileWrites = 0;
};

} // namespace femm

#endif
