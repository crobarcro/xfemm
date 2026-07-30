#include "FSolverAnalysisBackend.h"

#include "CBoundaryProp.h"
#include "CBlockLabel.h"
#include "CCircuit.h"
#include "CMaterialProp.h"
#include "CPointProp.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace femm {
namespace {

template<class Target, class Source>
Target magneticCopy(const std::unique_ptr<Source> &source, const char *description)
{
    const auto *typed = dynamic_cast<const Target *>(source.get());
    if (!typed)
        throw std::invalid_argument(std::string("magnetic model contains a non-magnetic ") + description);
    return *typed;
}

} // namespace

FSolverAnalysisBackend::FSolverAnalysisBackend() : m_solver(new FSolver) {}
FSolverAnalysisBackend::~FSolverAnalysisBackend() = default;

void FSolverAnalysisBackend::configure(const ModelDefinition &model,
                                       const SolveParameters &parameters,
                                       const PreparedAnalysis &prepared)
{
    const auto &problem = model.problem();
    m_solver->FileFormat = problem.FileFormat;
    m_solver->Frequency = parameters.frequency;
    m_solver->Precision = problem.Precision;
    m_solver->MinAngle = problem.MinAngle;
    m_solver->Depth = problem.Depth;
    m_solver->LengthUnits = problem.LengthUnits;
    m_solver->Coords = problem.Coords;
    m_solver->ProblemType = problem.problemType;
    m_solver->extZo = problem.extZo; m_solver->extRo = problem.extRo; m_solver->extRi = problem.extRi;
    m_solver->ACSolver = problem.ACSolver;
    m_solver->PrevType = problem.PrevType;
    m_solver->previousSolutionFile = problem.previousSolutionFile;
    m_solver->Relax = 1.0;

    m_solver->nodeproplist.clear();
    for (const auto &item : problem.nodeproplist)
        m_solver->nodeproplist.push_back(magneticCopy<CMPointProp>(item, "point property"));
    m_solver->lineproplist.clear();
    for (const auto &item : problem.lineproplist)
        m_solver->lineproplist.push_back(magneticCopy<CMBoundaryProp>(item, "boundary property"));
    m_solver->blockproplist.clear();
    for (const auto &item : prepared.materials) {
        CMSolverMaterialProp solverMaterial;
        static_cast<CMMaterialProp &>(solverMaterial) = item;
        m_solver->blockproplist.push_back(std::move(solverMaterial));
    }
    m_solver->labellist.clear();
    for (const auto &item : problem.labellist)
        m_solver->labellist.push_back(magneticCopy<CMBlockLabel>(item, "block label"));
    m_solver->circproplist.clear();
    for (std::size_t i = 0; i < problem.circproplist.size(); ++i) {
        auto circuit = magneticCopy<CMCircuit>(problem.circproplist[i], "circuit");
        const auto constraint = parameters.circuitConstraints.at(CircuitId{i});
        circuit.Amps = constraint.value;
        circuit.Case = constraint.kind == CircuitConstraintKind::PrescribedCurrent ? 0 : 2;
        circuit.dVolts = constraint.kind == CircuitConstraintKind::PrescribedVoltage
                       ? constraint.value : CComplex();
        m_solver->circproplist.push_back(circuit);
    }
    m_solver->NumPointProps = static_cast<int>(m_solver->nodeproplist.size());
    m_solver->NumLineProps = static_cast<int>(m_solver->lineproplist.size());
    m_solver->NumBlockProps = static_cast<int>(m_solver->blockproplist.size());
    m_solver->NumBlockLabels = static_cast<int>(m_solver->labellist.size());
    m_solver->NumCircProps = m_solver->NumCircPropsOrig = static_cast<int>(m_solver->circproplist.size());
}

void FSolverAnalysisBackend::positionAirGaps(const PreparedAnalysis &prepared)
{
    for (const auto &entry : prepared.airGapPositions)
        for (auto &gap : m_solver->agelist) {
            if (gap.BdryName != m_solver->lineproplist[entry.first.value].BdryName) continue;
            if (gap.innerRingTopology.empty() || gap.outerRingTopology.empty()) {
                gap.InnerAngle = entry.second.innerAngle;
                gap.OuterAngle = entry.second.outerAngle;
                continue;
            }
            const double step = gap.totalArcLength / gap.totalArcElements;
            auto positioned = [step](const std::vector<CQuadPoint> &topology, double angle) {
                auto ring = topology;
                for (auto &point : ring) {
                    point.w0 = std::fmod(point.w0 * step + angle, 360.0);
                    if (point.w0 < 0) point.w0 += 360.0;
                    point.w0 /= step;
                }
                std::stable_sort(ring.begin(), ring.end(), [](const CQuadPoint &a,
                                                              const CQuadPoint &b) {
                    return a.w0 < b.w0;
                });
                return ring;
            };
            const auto inner = positioned(gap.innerRingTopology, entry.second.innerAngle);
            const auto outer = positioned(gap.outerRingTopology, entry.second.outerAngle);
            const int fullCount = static_cast<int>(inner.size());
            gap.InnerShift = inner.front().w0;
            gap.OuterShift = outer.front().w0;
            gap.quadNode.clear();
            gap.quadNode.reserve(gap.totalArcElements + 1);
            for (int i = 0; i <= gap.totalArcElements; ++i) {
                const int p1 = i == fullCount ? 0 : i;
                const int p0 = p1 == 0 ? fullCount - 1 : p1 - 1;
                CQuadPoint q;
                q.n0=inner[p0].n0; q.n1=inner[p1].n0;
                q.n2=outer[p0].n0; q.n3=outer[p1].n0;
                q.w0=inner[p0].w1; q.w1=inner[p1].w1;
                q.w2=outer[p0].w1; q.w3=outer[p1].w1;
                gap.quadNode.push_back(q);
            }
            gap.InnerAngle = entry.second.innerAngle;
            gap.OuterAngle = entry.second.outerAngle;
            ++m_couplingRegenerations;
        }
}

void FSolverAnalysisBackend::synchronize(const ModelDefinition &model,
                                         const SolveParameters &parameters,
                                         const PreparedAnalysis &prepared,
                                         std::shared_ptr<const mesh::SolverMesh> mesh,
                                         std::uint64_t topologyIdentity, Dirty)
{
    if (!mesh)
        throw std::invalid_argument("FSolver backend requires a mesh");
    configure(model, parameters, prepared);
    if (topologyIdentity != m_topologyIdentity) {
        if (m_solver->LoadMesh(*mesh) != NOERROR)
            throw std::runtime_error("FSolver could not import the session mesh");
        std::vector<std::pair<std::size_t, std::size_t>> connectivity;
        connectivity.reserve(mesh->edges.size());
        for (const auto &edge : mesh->edges)
            connectivity.emplace_back(edge.first, edge.second);
        if (!m_solver->Cuthill(connectivity))
            throw std::runtime_error("FSolver Cuthill-McKee ordering failed");
        m_topologyIdentity = topologyIdentity;
        m_mesh = std::move(mesh);
        ++m_topologyImports;
        ++m_orderings;
    }
    positionAirGaps(prepared);
}

TrialSolution FSolverAnalysisBackend::solve(const ModelDefinition &model,
                                            const SolveParameters &parameters,
                                            const PreparedAnalysis &prepared)
{
    // Static2D/StaticAxisymmetric currently assemble both parts of a fresh
    // linear system on every evaluation.  Keep the two counters separate so
    // future dirty-flag based reuse can be introduced without changing the API.
    ++m_operatorAssemblies;
    ++m_rightHandSideAssemblies;
    configure(model, parameters, prepared);
    if (parameters.frequency != 0)
        throw std::invalid_argument("FSolverAnalysisBackend currently returns real (zero-frequency) solutions only");
    CBigLinProb system;
    system.Precision = m_solver->Precision;
    if (!system.Create(m_solver->NumNodes, m_solver->BandWidth))
        throw std::runtime_error("FSolver could not allocate the linear system");
    const bool solved = m_solver->ProblemType == PLANAR
                      ? m_solver->Static2D(system) : m_solver->StaticAxisymmetric(system);
    if (!solved)
        throw std::runtime_error("FSolver failed to solve the analysis");
    ++m_solves;

    TrialSolution result;
    result.real.emplace();
    result.real->nodal.magneticVectorPotential.reserve(m_solver->NumNodes);
    for (int i = 0; i < m_solver->NumNodes; ++i)
        result.real->nodal.magneticVectorPotential.push_back(system.b[i]);
    result.real->nodal.x.reserve(m_solver->NumNodes);
    result.real->nodal.y.reserve(m_solver->NumNodes);
    for (const auto &node : m_solver->meshnode) {
        // FSolver stores imported coordinates internally in centimetres.
        result.real->nodal.x.push_back(node.x / 100.0);
        result.real->nodal.y.push_back(node.y / 100.0);
    }
    for (std::size_t i = 0; i < m_solver->circproplist.size(); ++i) {
        const auto &constraint = parameters.circuitConstraints.at(CircuitId{i});
        CComplex current = constraint.kind == CircuitConstraintKind::PrescribedCurrent
                         ? constraint.value : m_solver->circproplist[i].Amps;
        std::optional<CComplex> voltage;
        if (constraint.kind == CircuitConstraintKind::PrescribedVoltage)
            voltage = constraint.value;
        result.circuits.push_back({CircuitId{i}, current, CComplex(), voltage});
        std::optional<double> realVoltage;
        if (voltage)
            realVoltage = voltage->re;
        result.real->circuits.push_back({CircuitId{i}, current.re, 0.0, realVoltage});
    }
    return result;
}

} // namespace femm
