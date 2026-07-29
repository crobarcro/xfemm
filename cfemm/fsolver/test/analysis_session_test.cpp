#include "AnalysisSession.h"
#include "FSolverAnalysisBackend.h"

#include "CBoundaryProp.h"
#include "CBlockLabel.h"
#include "CCircuit.h"
#include "CMaterialProp.h"
#include "MesherBackend.h"

#include <cassert>
#include <memory>
#include <stdexcept>

namespace {

class RecordingBackend final : public femm::AnalysisSolverBackend {
public:
    void synchronize(const femm::ModelDefinition &, const femm::SolveParameters &,
                     const femm::PreparedAnalysis &,
                     std::shared_ptr<const femm::mesh::SolverMesh> mesh,
                     std::uint64_t topology, femm::Dirty rebuilt) override
    {
        lastMesh = std::move(mesh);
        lastTopology = topology;
        lastRebuilt = rebuilt;
        ++synchronizations;
    }

    femm::TrialSolution solve(const femm::ModelDefinition &,
                              const femm::SolveParameters &parameters,
                              const femm::PreparedAnalysis &) override
    {
        ++solves;
        femm::TrialSolution result;
        result.backendState = "accepted-field-state";
        for (const auto &entry : parameters.circuitConstraints) {
            if (entry.second.kind == femm::CircuitConstraintKind::PrescribedCurrent)
                result.circuits.push_back({entry.first, entry.second.value,
                                           entry.second.value * 2.0, std::nullopt});
        }
        return result;
    }

    femm::Dirty lastRebuilt = femm::Dirty::None;
    int synchronizations = 0;
    int solves = 0;
    std::shared_ptr<const femm::mesh::SolverMesh> lastMesh;
    std::uint64_t lastTopology = 0;
};

class RecordingMesher final : public fmesher::MesherBackend {
public:
    femm::mesh::MeshResult mesh(femm::FemmProblem &, bool periodic,
                                const femm::mesh::MeshingOptions &) override
    {
        ++calls;
        lastPeriodic = periodic;
        femm::mesh::MeshResult result;
        result.status = femm::mesh::MeshStatus::Success;
        result.mesh.nodes = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}};
        result.mesh.elements.push_back({{{0, 1, 2}}, 1});
        result.mesh.edges = {{0, 1, 0}, {1, 2, 0}, {2, 0, 0}};
        return result;
    }
    int calls = 0;
    bool lastPeriodic = false;
};

std::unique_ptr<femm::FemmProblem> makeProblem()
{
    auto problem = std::make_unique<femm::FemmProblem>(femm::FileType::MagneticsFile);

    auto material = std::make_unique<femm::CMMaterialProp>();
    material->BlockName = "steel";
    material->mu_x = 100;
    material->mu_y = 100;
    problem->blockproplist.push_back(std::move(material));

    auto circuit = std::make_unique<femm::CMCircuit>();
    circuit->CircName = "phase-a";
    circuit->CircType = 1;
    circuit->Amps = 3;
    problem->circproplist.push_back(std::move(circuit));

    auto first = std::make_unique<femm::CMBlockLabel>();
    first->InCircuit = 0;
    first->Turns = 10;
    problem->labellist.push_back(std::move(first));
    auto second = std::make_unique<femm::CMBlockLabel>();
    second->InCircuit = 0;
    second->Turns = -5;
    problem->labellist.push_back(std::move(second));

    auto gap = std::make_unique<femm::CMBoundaryProp>();
    gap->BdryName = "rotor-gap";
    gap->BdryFormat = 6;
    gap->InnerAngle = 1;
    gap->OuterAngle = 2;
    problem->lineproplist.push_back(std::move(gap));
    return problem;
}

} // namespace

int main()
{
    auto backend = std::make_shared<RecordingBackend>();
    auto mesher = std::make_shared<RecordingMesher>();
    femm::AnalysisSession session(femm::ModelDefinition(makeProblem()), mesher, backend);

    const auto circuit = session.model().circuit("phase-a");
    const auto material = session.model().material("steel");
    const auto gap = session.model().airGap("rotor-gap");

    session.synchronize();
    assert(backend->synchronizations == 1);
    assert(mesher->calls == 1);
    assert(backend->lastMesh == session.mesh());
    assert(backend->lastTopology == session.meshTopologyIdentity());
    assert(session.prepared().circuits.size() == 2);
    assert(session.prepared().circuits[0].signedTurns == 10);
    assert(session.prepared().circuits[1].signedTurns == -5);
    // Preparing series windings must not rewrite the authoritative labels or circuit.
    assert(session.model().problem().labellist[0]->InCircuit == 0);
    assert(dynamic_cast<const femm::CMCircuit *>(
               session.model().problem().circproplist[0].get())->CircType == 1);

    session.setCircuitCurrent(circuit, CComplex(7, 0));
    session.synchronize();
    assert(backend->synchronizations == 2);
    assert(mesher->calls == 1); // circuit changes do not affect topology.
    assert((session.dirty() & femm::Dirty::SolveState) != femm::Dirty::None);
    assert(session.prepared().circuits[0].constraint.value == CComplex(7, 0));

    session.setAirGapAngle(gap, 12, 2);
    session.setFrequency(50);
    session.setMaterialProperty(material, femm::MaterialProperty::Conductivity, 5.8);
    session.synchronize();
    assert(session.prepared().frequency == 50);
    assert(session.prepared().airGapPositions.at(gap).innerAngle == 12);
    assert(session.model().problem().Frequency == 0); // not copied to the model.

    const auto trial = session.solve();
    assert(backend->solves == 1);
    assert(trial.id != 0);
    assert(trial.circuits[0].fluxLinkage == CComplex(14, 0));
    const auto accepted = session.acceptSolution(trial);
    assert(accepted->backendState == "accepted-field-state");
    assert(session.solveParameters().initialState == accepted);

    session.setTime(0.01);
    bool staleRejected = false;
    try {
        session.acceptSolution(trial);
    } catch (const std::invalid_argument &) {
        staleRejected = true;
    }
    assert(staleRejected);

    const auto before = session.solveParameters();
    bool invalidBatchRejected = false;
    try {
        session.updateSolveParameters([](femm::SolveParameters &parameters) {
            parameters.frequency = -1;
        });
    } catch (const std::invalid_argument &) {
        invalidBatchRejected = true;
    }
    assert(invalidBatchRejected);
    assert(session.solveParameters().frequency == before.frequency);

    auto concrete = std::make_shared<femm::FSolverAnalysisBackend>();
    auto concreteMesher = std::make_shared<RecordingMesher>();
    femm::AnalysisSession concreteSession(femm::ModelDefinition(makeProblem()),
                                          concreteMesher, concrete);
    concreteSession.synchronize();
    assert(concreteSession.meshGenerationCount() == 1);
    assert(concrete->topologyImportCount() == 1);
    assert(concrete->orderingCount() == 1);
    constexpr int evaluations = 3;
    for (int step = 0; step < evaluations; ++step) {
        concreteSession.setCircuitCurrent(concreteSession.model().circuit("phase-a"),
                                          CComplex(9 + step, 0));
        concreteSession.setAirGapAngle(concreteSession.model().airGap("rotor-gap"),
                                       4 + step, 5);
        concreteSession.setTime(0.01 * step);
        concreteSession.solve();
    }
    assert(concreteMesher->calls == 1);
    assert(concreteSession.meshGenerationCount() == 1);
    assert(concrete->topologyImportCount() == 1);
    assert(concrete->orderingCount() == 1);
    assert(concrete->solveCount() == evaluations);
    // Assembly reuse is intentionally not part of this contract yet.
    assert(concrete->operatorAssemblyCount() >= 1);
    assert(concrete->operatorAssemblyCount() <= evaluations);
    assert(concrete->rightHandSideAssemblyCount() >= 1);
    assert(concrete->rightHandSideAssemblyCount() <= evaluations);
    assert(concrete->meshFileWriteCount() == 0);
    assert(concrete->meshFileReadCount() == 0);
}
