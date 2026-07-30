#ifndef XFEMM_ANALYSISSESSION_H
#define XFEMM_ANALYSISSESSION_H

#include "FemmProblem.h"
#include "mesh/Meshing.h"

#include <cstdint>
#include <functional>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <variant>
#include <vector>

namespace femm {

struct CircuitId {
    std::size_t value;
    friend bool operator<(CircuitId a, CircuitId b) { return a.value < b.value; }
    friend bool operator==(CircuitId a, CircuitId b) { return a.value == b.value; }
};

struct MaterialId {
    std::size_t value;
};

struct AirGapId {
    std::size_t value;
    friend bool operator<(AirGapId a, AirGapId b) { return a.value < b.value; }
};

/** User-authored magnetic model. Mutable access is deliberately limited to AnalysisSession. */
class ModelDefinition {
public:
    /** Transfers exclusive ownership so no caller can bypass session dirty tracking. */
    explicit ModelDefinition(std::unique_ptr<FemmProblem> problem);
    const FemmProblem &problem() const { return *m_problem; }

    CircuitId circuit(const std::string &name) const;
    MaterialId material(const std::string &name) const;
    AirGapId airGap(const std::string &name) const;

private:
    friend class AnalysisSession;
    std::unique_ptr<FemmProblem> m_problem;
};

enum class CircuitConstraintKind {
    PrescribedCurrent,
    PrescribedVoltage,
    OpenCircuit,
    CoupledUnknown
};

struct CircuitConstraint {
    CircuitConstraintKind kind = CircuitConstraintKind::OpenCircuit;
    CComplex value;
};

struct AirGapPosition {
    double innerAngle = 0;
    double outerAngle = 0;
};

struct AcceptedState {
    std::uint64_t modelRevision = 0;
    std::uint64_t parameterRevision = 0;
    double time = 0;
    std::string backendState;
};

/** Inputs which may change between evaluations without changing the physical model. */
struct SolveParameters {
    double frequency = 0;
    double time = 0;
    std::map<CircuitId, CircuitConstraint> circuitConstraints;
    std::map<AirGapId, AirGapPosition> airGapPositions;
    std::shared_ptr<const AcceptedState> initialState;
};

enum class MaterialProperty {
    RelativePermeabilityX,
    RelativePermeabilityY,
    Coercivity,
    Conductivity,
    LaminationType,
    LaminationFill,
    LaminationThickness,
    AppliedCurrentDensity,
    BHCurve
};

struct BHCurve {
    std::vector<double> fluxDensity;
    std::vector<CComplex> fieldStrength;
};

using MaterialPropertyValue = std::variant<double, int, CComplex, BHCurve>;

enum class Dirty : std::uint32_t {
    None = 0,
    ModelIndex = 1u << 0,
    Mesh = 1u << 1,
    PreparedMaterials = 1u << 2,
    PreparedCircuits = 1u << 3,
    AirGapCoupling = 1u << 4,
    Operator = 1u << 5,
    RightHandSide = 1u << 6,
    InitialState = 1u << 7,
    SolveState = 1u << 8,
    All = (1u << 9) - 1
};

Dirty operator|(Dirty lhs, Dirty rhs);
Dirty operator&(Dirty lhs, Dirty rhs);

struct PreparedCircuit {
    CircuitId source;
    std::size_t labelIndex;
    int signedTurns;
    CircuitConstraint constraint;
};

/** Disposable solver input. It is recreated from authoritative user data. */
struct PreparedAnalysis {
    std::vector<CMMaterialProp> materials;
    std::vector<PreparedCircuit> circuits;
    std::map<AirGapId, AirGapPosition> airGapPositions;
    double frequency = 0;
};

struct CircuitPortResult {
    CircuitId id;
    CComplex current;
    CComplex fluxLinkage;
    std::optional<CComplex> terminalVoltage;
};

/** Real-valued magnetic vector potential in solver node order. */
struct RealNodalSolution {
    std::vector<double> magneticVectorPotential;
    /** Mesh coordinates in the same (possibly reordered) solver-node order. */
    std::vector<double> x;
    std::vector<double> y;
};

struct RealCircuitPortResult {
    CircuitId id;
    double current = 0;
    double fluxLinkage = 0;
    std::optional<double> terminalVoltage;
};

/** Strongly typed result payload for magnetostatic analyses. */
struct RealTrialSolution {
    RealNodalSolution nodal;
    std::vector<RealCircuitPortResult> circuits;
};

struct TrialSolution {
    std::uint64_t id = 0;
    std::uint64_t sessionId = 0;
    std::uint64_t modelRevision = 0;
    std::uint64_t parameterRevision = 0;
    double time = 0;
    std::vector<CircuitPortResult> circuits;
    std::optional<RealTrialSolution> real;
    std::string backendState;
};

/** Adapter implemented by a concrete field solver; it never receives mutable model state. */
class AnalysisSolverBackend {
public:
    virtual ~AnalysisSolverBackend() = default;
    virtual void synchronize(const ModelDefinition &model,
                             const SolveParameters &parameters,
                             const PreparedAnalysis &prepared,
                             std::shared_ptr<const mesh::SolverMesh> mesh,
                             std::uint64_t meshTopologyIdentity,
                             Dirty rebuilt) = 0;
    virtual TrialSolution solve(const ModelDefinition &model,
                                const SolveParameters &parameters,
                                const PreparedAnalysis &prepared) = 0;
};

} // namespace femm

namespace fmesher { class MesherBackend; }

namespace femm {

class AnalysisSession {
public:
    AnalysisSession(ModelDefinition model, std::shared_ptr<AnalysisSolverBackend> backend);
    AnalysisSession(ModelDefinition model, std::shared_ptr<fmesher::MesherBackend> mesher,
                    std::shared_ptr<AnalysisSolverBackend> backend);

    const ModelDefinition &model() const { return m_model; }
    const SolveParameters &solveParameters() const { return m_parameters; }
    const PreparedAnalysis &prepared() const { return m_prepared; }
    Dirty dirty() const { return m_dirty; }
    const mesh::MeshingOptions &meshingOptions() const { return m_meshingOptions; }
    const std::vector<mesh::MeshDiagnostic> &meshDiagnostics() const { return m_meshDiagnostics; }
    std::shared_ptr<const mesh::SolverMesh> mesh() const { return m_mesh; }
    std::uint64_t meshTopologyIdentity() const { return m_meshTopologyIdentity; }
    std::size_t meshGenerationCount() const { return m_meshGenerations; }

    /** Select a mesher. The currently owned mesh is discarded. */
    void setMesher(std::shared_ptr<fmesher::MesherBackend> mesher);
    /** Change controls used for the next mesh. */
    void setMeshingOptions(const mesh::MeshingOptions &options);
    /** Return the current immutable mesh, creating it when necessary. */
    std::shared_ptr<const mesh::SolverMesh> ensureMesh();

    void setCircuitCurrent(CircuitId id, CComplex current);
    void setCircuitVoltage(CircuitId id, CComplex voltage);
    void setCircuitOpen(CircuitId id);
    void setCircuitCoupled(CircuitId id);
    void setAirGapAngle(AirGapId id, double innerAngle, double outerAngle);
    void setFrequency(double frequency);
    void setTime(double time);
    void setInitialState(std::shared_ptr<const AcceptedState> state);
    void setMaterialProperty(MaterialId id, MaterialProperty property,
                             const MaterialPropertyValue &value);

    /** Apply several changes atomically and coalesce their invalidations. */
    void updateSolveParameters(const std::function<void(SolveParameters &)> &update);

    void synchronize();
    TrialSolution solve();
    std::shared_ptr<const AcceptedState> acceptSolution(const TrialSolution &solution);

private:
    void requireCircuit(CircuitId id) const;
    void requireAirGap(AirGapId id) const;
    void invalidate(Dirty dirty);
    void rebuildMaterials(PreparedAnalysis &candidate) const;
    void rebuildCircuits(PreparedAnalysis &candidate) const;

    ModelDefinition m_model;
    SolveParameters m_parameters;
    PreparedAnalysis m_prepared;
    std::shared_ptr<AnalysisSolverBackend> m_backend;
    std::shared_ptr<fmesher::MesherBackend> m_mesher;
    mesh::MeshingOptions m_meshingOptions;
    std::shared_ptr<const mesh::SolverMesh> m_mesh;
    std::vector<mesh::MeshDiagnostic> m_meshDiagnostics;
    std::uint64_t m_meshTopologyIdentity = 0;
    std::size_t m_meshGenerations = 0;
    Dirty m_dirty = Dirty::All;
    std::uint64_t m_modelRevision = 1;
    std::uint64_t m_parameterRevision = 1;
    std::uint64_t m_nextSolutionId = 1;
    std::uint64_t m_sessionId = 0;
};

} // namespace femm

#endif
