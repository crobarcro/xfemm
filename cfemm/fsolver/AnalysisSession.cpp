#include "AnalysisSession.h"

#include "CBoundaryProp.h"
#include "CBlockLabel.h"
#include "CCircuit.h"
#include "CMaterialProp.h"
#include "MesherBackend.h"
#include "TangleMesherBackend.h"
#include "TriangleMesherBackend.h"

#include <cmath>
#include <atomic>
#include <stdexcept>

namespace femm {
namespace {

std::atomic<std::uint64_t> nextSessionId{1};

std::shared_ptr<fmesher::MesherBackend> makeDefaultMesher()
{
    // Tangle is the tested production default for sessions. Triangle remains
    // available through explicit injection for compatibility and baselines.
    return std::make_shared<fmesher::TangleMesherBackend>();
}

std::uint32_t bits(Dirty value) { return static_cast<std::uint32_t>(value); }
bool has(Dirty value, Dirty flag) { return bits(value & flag) != 0; }

bool sameComplex(const CComplex &lhs, const CComplex &rhs)
{
    return lhs.re == rhs.re && lhs.im == rhs.im;
}

bool sameConstraints(const std::map<CircuitId, CircuitConstraint> &lhs,
                     const std::map<CircuitId, CircuitConstraint> &rhs)
{
    if (lhs.size() != rhs.size())
        return false;
    auto lit = lhs.begin();
    auto rit = rhs.begin();
    for (; lit != lhs.end(); ++lit, ++rit) {
        if (!(lit->first == rit->first) || lit->second.kind != rit->second.kind ||
            !sameComplex(lit->second.value, rit->second.value))
            return false;
    }
    return true;
}

bool sameAirGaps(const std::map<AirGapId, AirGapPosition> &lhs,
                 const std::map<AirGapId, AirGapPosition> &rhs)
{
    if (lhs.size() != rhs.size())
        return false;
    auto lit = lhs.begin();
    auto rit = rhs.begin();
    for (; lit != lhs.end(); ++lit, ++rit) {
        if (lit->first.value != rit->first.value ||
            lit->second.innerAngle != rit->second.innerAngle ||
            lit->second.outerAngle != rit->second.outerAngle)
            return false;
    }
    return true;
}

template<typename T>
const T &as(const MaterialPropertyValue &value, const char *description)
{
    const auto *typed = std::get_if<T>(&value);
    if (!typed)
        throw std::invalid_argument(std::string("wrong value type for ") + description);
    return *typed;
}

} // namespace

Dirty operator|(Dirty lhs, Dirty rhs)
{
    return static_cast<Dirty>(bits(lhs) | bits(rhs));
}

Dirty operator&(Dirty lhs, Dirty rhs)
{
    return static_cast<Dirty>(bits(lhs) & bits(rhs));
}

ModelDefinition::ModelDefinition(std::unique_ptr<FemmProblem> problem)
    : m_problem(std::move(problem))
{
    if (!m_problem || m_problem->filetype != FileType::MagneticsFile)
        throw std::invalid_argument("ModelDefinition requires a magnetic FemmProblem");
    m_problem->updateBlockMap();
    m_problem->updateCircuitMap();
    m_problem->updateLineMap();
}

CircuitId ModelDefinition::circuit(const std::string &name) const
{
    const auto found = m_problem->circuitMap.find(name);
    if (found == m_problem->circuitMap.end())
        throw std::out_of_range("unknown circuit: " + name);
    return {static_cast<std::size_t>(found->second)};
}

MaterialId ModelDefinition::material(const std::string &name) const
{
    const auto found = m_problem->blockMap.find(name);
    if (found == m_problem->blockMap.end())
        throw std::out_of_range("unknown material: " + name);
    return {static_cast<std::size_t>(found->second)};
}

AirGapId ModelDefinition::airGap(const std::string &name) const
{
    const auto found = m_problem->lineMap.find(name);
    if (found == m_problem->lineMap.end())
        throw std::out_of_range("unknown boundary: " + name);
    const auto id = static_cast<std::size_t>(found->second);
    const auto *boundary = dynamic_cast<const CMBoundaryProp *>(m_problem->lineproplist[id].get());
    if (!boundary || (boundary->BdryFormat != 6 && boundary->BdryFormat != 7))
        throw std::invalid_argument("boundary is not an air-gap element: " + name);
    return {id};
}

AnalysisSession::AnalysisSession(ModelDefinition model, std::shared_ptr<AnalysisSolverBackend> backend)
    : m_model(std::move(model)), m_backend(std::move(backend)),
      m_mesher(makeDefaultMesher()),
      m_sessionId(nextSessionId.fetch_add(1))
{
    if (!m_backend)
        throw std::invalid_argument("AnalysisSession requires a solver backend");
    m_parameters.frequency = m_model.problem().Frequency;

    for (std::size_t i = 0; i < m_model.problem().circproplist.size(); ++i) {
        const auto *circuit = dynamic_cast<const CMCircuit *>(m_model.problem().circproplist[i].get());
        if (!circuit)
            throw std::invalid_argument("magnetic model contains a non-magnetic circuit");
        m_parameters.circuitConstraints[{i}] = {CircuitConstraintKind::PrescribedCurrent,
                                                circuit->Amps};
    }
    for (std::size_t i = 0; i < m_model.problem().lineproplist.size(); ++i) {
        const auto *boundary = dynamic_cast<const CMBoundaryProp *>(m_model.problem().lineproplist[i].get());
        if (boundary && (boundary->BdryFormat == 6 || boundary->BdryFormat == 7))
            m_parameters.airGapPositions[{i}] = {boundary->InnerAngle, boundary->OuterAngle};
    }
}

AnalysisSession::AnalysisSession(ModelDefinition model,
                                 std::shared_ptr<fmesher::MesherBackend> mesher,
                                 std::shared_ptr<AnalysisSolverBackend> backend)
    : AnalysisSession(std::move(model), std::move(backend))
{
    setMesher(std::move(mesher));
}

const char *AnalysisSession::mesherBackendName() const
{
    return m_mesher ? m_mesher->name() : "none";
}

void AnalysisSession::setMesher(std::shared_ptr<fmesher::MesherBackend> mesher)
{
    if (!mesher)
        throw std::invalid_argument("AnalysisSession requires a mesher backend");
    if (m_mesher == mesher)
        return;
    m_mesher = std::move(mesher);
    m_mesh.reset();
    m_meshDiagnostics.clear();
    invalidate(Dirty::Mesh | Dirty::Operator | Dirty::RightHandSide);
}

void AnalysisSession::setMeshingOptions(const mesh::MeshingOptions &options)
{
    if (!std::isfinite(options.minimumAngleDegrees) || options.minimumAngleDegrees < 0 ||
        !std::isfinite(options.defaultElementSize) || options.defaultElementSize < 0)
        throw std::invalid_argument("meshing sizes and angles must be finite and nonnegative");
    m_meshingOptions = options;
    m_mesh.reset();
    m_meshDiagnostics.clear();
    invalidate(Dirty::Mesh | Dirty::Operator | Dirty::RightHandSide);
}

std::shared_ptr<const mesh::SolverMesh> AnalysisSession::ensureMesh()
{
    if (m_instanced)
        return ensureInstancedMesh();
    if (m_mesh && !has(m_dirty, Dirty::Mesh))
        return m_mesh;
    if (!m_mesher)
        throw std::logic_error("no mesher backend selected");

    bool periodic = false;
    for (const auto &property : m_model.problem().lineproplist) {
        if (property && property->isPeriodic()) {
            periodic = true;
            break;
        }
    }
    auto result = m_mesher->mesh(*m_model.m_problem, periodic, m_meshingOptions);
    ++m_meshGenerations;
    m_meshDiagnostics = std::move(result.diagnostics);
    if (!result.succeeded())
        throw std::runtime_error("meshing failed");

    m_mesh = std::make_shared<const mesh::SolverMesh>(std::move(result.mesh));
    ++m_meshTopologyIdentity;
    m_dirty = static_cast<Dirty>(bits(m_dirty) & ~bits(Dirty::Mesh));
    return m_mesh;
}

void AnalysisSession::setInstancedMesh(mesh::InstancedMesh instanced)
{
    m_instanced = std::move(instanced);
    m_canonicalMesh.reset();
    m_instancingProvenance = {};
    m_templateTopologyIdentity = 0;
    m_instanceLayoutIdentity = 0;
    m_mesh.reset();
    m_meshDiagnostics.clear();
    invalidate(Dirty::Mesh | Dirty::PreparedCircuits | Dirty::Operator | Dirty::RightHandSide);
}

void AnalysisSession::clearInstancedMesh()
{
    if (!m_instanced)
        return;
    m_instanced.reset();
    m_canonicalMesh.reset();
    m_instancingProvenance = {};
    m_templateTopologyIdentity = 0;
    m_instanceLayoutIdentity = 0;
    m_mesh.reset();
    m_meshDiagnostics.clear();
    invalidate(Dirty::Mesh | Dirty::PreparedCircuits | Dirty::Operator | Dirty::RightHandSide);
}

void AnalysisSession::setInstanceTransform(std::size_t instance,
                                           const mesh::RigidTransform2D &transform)
{
    if (!m_instanced || instance >= m_instanced->instances.size())
        throw std::out_of_range("invalid instance index");
    if (!transform.isFinite())
        throw std::invalid_argument("instance transform must be finite");
    m_instanced->instances[instance].transform = transform;
    m_mesh.reset();
    // The transform rotates constant magnetisation, so prepared labels must be
    // rebuilt as well as the mesh and operator.
    invalidate(Dirty::Mesh | Dirty::PreparedCircuits | Dirty::Operator |
               Dirty::RightHandSide);
}

void AnalysisSession::setInstanceRegionOverrides(
    std::size_t instance, std::vector<mesh::InstanceRegionOverride> overrides)
{
    if (!m_instanced || instance >= m_instanced->instances.size())
        throw std::out_of_range("invalid instance index");
    for (const auto &override : overrides) {
        if (override.magnetisationRotationDegrees &&
            !std::isfinite(*override.magnetisationRotationDegrees))
            throw std::invalid_argument("magnetisation rotation must be finite");
        if (override.currentScale && !std::isfinite(*override.currentScale))
            throw std::invalid_argument("current scale must be finite");
        if (override.circuit)
            requireCircuit(CircuitId{*override.circuit});
    }
    m_instanced->instances[instance].regionOverrides = std::move(overrides);
    // Overrides are physics only: the topology and solver mesh are unchanged.
    invalidate(Dirty::PreparedCircuits | Dirty::Operator | Dirty::RightHandSide);
}

void AnalysisSession::addInstanceRegionOverride(
    std::size_t instance, const mesh::InstanceRegionOverride &override)
{
    if (!m_instanced || instance >= m_instanced->instances.size())
        throw std::out_of_range("invalid instance index");
    if (override.magnetisationRotationDegrees &&
        !std::isfinite(*override.magnetisationRotationDegrees))
        throw std::invalid_argument("magnetisation rotation must be finite");
    if (override.currentScale && !std::isfinite(*override.currentScale))
        throw std::invalid_argument("current scale must be finite");
    if (override.circuit)
        requireCircuit(CircuitId{*override.circuit});
    m_instanced->instances[instance].regionOverrides.push_back(override);
    invalidate(Dirty::PreparedCircuits | Dirty::Operator | Dirty::RightHandSide);
}

std::optional<mesh::ElementProvenance>
AnalysisSession::elementProvenance(std::size_t globalElement) const
{
    if (globalElement >= m_instancingProvenance.elementProvenance.size())
        return std::nullopt;
    return m_instancingProvenance.elementProvenance[globalElement];
}

std::optional<mesh::NodeProvenance>
AnalysisSession::nodeProvenance(std::size_t globalNode) const
{
    if (globalNode >= m_instancingProvenance.nodeProvenance.size())
        return std::nullopt;
    return m_instancingProvenance.nodeProvenance[globalNode];
}

std::vector<mesh::MeshIndex>
AnalysisSession::elementsForInstance(std::size_t instanceIndex) const
{
    std::vector<mesh::MeshIndex> elements;
    const auto &provenance = m_instancingProvenance.elementProvenance;
    for (mesh::MeshIndex i = 0; i < provenance.size(); ++i)
        if (provenance[i].instanceIndex == instanceIndex)
            elements.push_back(i);
    return elements;
}

std::vector<mesh::MeshIndex>
AnalysisSession::nodesForInstance(std::size_t instanceIndex) const
{
    std::vector<mesh::MeshIndex> nodes;
    const auto &provenance = m_instancingProvenance.nodeProvenance;
    for (mesh::MeshIndex i = 0; i < provenance.size(); ++i)
        if (provenance[i].instanceIndex == instanceIndex)
            nodes.push_back(i);
    return nodes;
}

std::size_t AnalysisSession::templateLabelCount() const
{
    std::size_t count = 0;
    for (const auto &property : m_model.problem().labellist) {
        const auto *label = dynamic_cast<const CMBlockLabel *>(property.get());
        if (label && !label->isHole())
            ++count;
    }
    return count;
}

mesh::SolverMesh
AnalysisSession::remapInstancedRegions(const mesh::SolverMesh &canonical) const
{
    mesh::SolverMesh mesh = canonical;
    const std::size_t labelCount = templateLabelCount();
    if (labelCount == 0)
        return mesh;
    for (std::size_t e = 0; e < mesh.elements.size(); ++e) {
        const std::size_t instance = m_instancingProvenance.elementProvenance[e].instanceIndex;
        const std::int32_t attribute = mesh.elements[e].regionAttribute;
        if (attribute > 0 && static_cast<std::size_t>(attribute) <= labelCount)
            mesh.elements[e].regionAttribute =
                static_cast<std::int32_t>(instance * labelCount + static_cast<std::size_t>(attribute));
    }
    return mesh;
}

std::shared_ptr<const mesh::SolverMesh> AnalysisSession::ensureInstancedMesh()
{
    if (m_mesh && !has(m_dirty, Dirty::Mesh))
        return m_mesh;
    const std::uint64_t templateId = mesh::templateTopologyIdentity(*m_instanced);
    const std::uint64_t layoutId = mesh::instanceLayoutIdentity(*m_instanced);
    if (!m_canonicalMesh || templateId != m_templateTopologyIdentity ||
        layoutId != m_instanceLayoutIdentity) {
        auto result = m_instanced->materialize();
        m_meshDiagnostics.clear();
        if (!result.succeeded()) {
            for (const auto &diagnostic : result.diagnostics)
                m_meshDiagnostics.push_back(
                    {mesh::MeshDiagnosticSeverity::Error, diagnostic.message, "InstancedMesh",
                     static_cast<int>(diagnostic.category)});
            throw std::runtime_error("instanced mesh materialisation failed");
        }
        m_canonicalMesh =
            std::make_shared<const mesh::SolverMesh>(std::move(result.mesh));
        m_instancingProvenance = std::move(result.provenance);
        m_templateTopologyIdentity = templateId;
        m_instanceLayoutIdentity = layoutId;
        ++m_materializationCount;
        ++m_meshGenerations;
    }
    m_mesh = std::make_shared<const mesh::SolverMesh>(remapInstancedRegions(*m_canonicalMesh));
    m_meshTopologyIdentity = mesh::materializedTopologyIdentity(*m_mesh);
    m_dirty = static_cast<Dirty>(bits(m_dirty) & ~bits(Dirty::Mesh));
    return m_mesh;
}

void AnalysisSession::rebuildInstancedPrepared(PreparedAnalysis &candidate) const
{
    candidate.labels.clear();
    candidate.circuits.clear();
    if (!m_instanced)
        return;

    std::vector<std::size_t> rawIndices;
    std::vector<const CMBlockLabel *> templateLabels;
    for (std::size_t raw = 0; raw < m_model.problem().labellist.size(); ++raw) {
        const auto *label =
            dynamic_cast<const CMBlockLabel *>(m_model.problem().labellist[raw].get());
        if (label && !label->isHole()) {
            templateLabels.push_back(label);
            rawIndices.push_back(raw);
        }
    }
    const std::size_t labelCount = templateLabels.size();
    const std::size_t instanceCount = m_instanced->instances.size();
    candidate.labels.reserve(labelCount * instanceCount);
    for (std::size_t k = 0; k < instanceCount; ++k) {
        const auto &instance = m_instanced->instances[k];
        for (std::size_t j = 0; j < labelCount; ++j) {
            CMBlockLabel label = *templateLabels[j];
            // A rigid transform rotates directional quantities. A constant
            // magnetisation direction therefore follows the instance rotation;
            // a positional MagDirFctn is evaluated at the transformed
            // coordinates and must not be rotated again. Overrides are applied
            // afterwards as an explicit delta.
            if (label.MagDirFctn.empty())
                label.MagDir += instance.transform.rotationDegrees;
            for (const auto &override : instance.regionOverrides) {
                if (override.sourceBlockLabel != rawIndices[j])
                    continue;
                if (override.circuit)
                    label.InCircuit = static_cast<int>(*override.circuit);
                if (override.magnetisationRotationDegrees)
                    label.MagDir += *override.magnetisationRotationDegrees;
                if (override.currentScale)
                    label.Turns = static_cast<int>(
                        std::lround(label.Turns * *override.currentScale));
            }
            candidate.labels.push_back(std::move(label));
        }
    }

    for (std::size_t i = 0; i < candidate.labels.size(); ++i) {
        const auto &label = candidate.labels[i];
        if (label.InCircuit < 0)
            continue;
        CircuitId id{static_cast<std::size_t>(label.InCircuit)};
        requireCircuit(id);
        const auto constraint = m_parameters.circuitConstraints.at(id);
        candidate.circuits.push_back({id, i, label.Turns, constraint});
    }
}

void AnalysisSession::requireCircuit(CircuitId id) const
{
    if (id.value >= m_model.problem().circproplist.size())
        throw std::out_of_range("invalid circuit id");
}

void AnalysisSession::requireAirGap(AirGapId id) const
{
    if (id.value >= m_model.problem().lineproplist.size())
        throw std::out_of_range("invalid air-gap id");
    const auto *boundary = dynamic_cast<const CMBoundaryProp *>(m_model.problem().lineproplist[id.value].get());
    if (!boundary || (boundary->BdryFormat != 6 && boundary->BdryFormat != 7))
        throw std::invalid_argument("boundary id does not identify an air-gap element");
}

void AnalysisSession::invalidate(Dirty dirty)
{
    m_dirty = m_dirty | dirty | Dirty::SolveState;
}

void AnalysisSession::setCircuitCurrent(CircuitId id, CComplex value)
{
    requireCircuit(id);
    m_parameters.circuitConstraints[id] = {CircuitConstraintKind::PrescribedCurrent, value};
    ++m_parameterRevision;
    invalidate(Dirty::RightHandSide);
}

void AnalysisSession::setCircuitVoltage(CircuitId id, CComplex value)
{
    requireCircuit(id);
    m_parameters.circuitConstraints[id] = {CircuitConstraintKind::PrescribedVoltage, value};
    ++m_parameterRevision;
    invalidate(Dirty::RightHandSide | Dirty::Operator);
}

void AnalysisSession::setCircuitOpen(CircuitId id)
{
    requireCircuit(id);
    m_parameters.circuitConstraints[id] = {CircuitConstraintKind::OpenCircuit, {}};
    ++m_parameterRevision;
    invalidate(Dirty::RightHandSide | Dirty::Operator);
}

void AnalysisSession::setCircuitCoupled(CircuitId id)
{
    requireCircuit(id);
    m_parameters.circuitConstraints[id] = {CircuitConstraintKind::CoupledUnknown, {}};
    ++m_parameterRevision;
    invalidate(Dirty::RightHandSide | Dirty::Operator);
}

void AnalysisSession::setAirGapAngle(AirGapId id, double inner, double outer)
{
    requireAirGap(id);
    if (!std::isfinite(inner) || !std::isfinite(outer))
        throw std::invalid_argument("air-gap angles must be finite");
    m_parameters.airGapPositions[id] = {inner, outer};
    ++m_parameterRevision;
    invalidate(Dirty::AirGapCoupling | Dirty::Operator);
}

void AnalysisSession::setFrequency(double frequency)
{
    if (!std::isfinite(frequency) || frequency < 0)
        throw std::invalid_argument("frequency must be finite and nonnegative");
    m_parameters.frequency = frequency;
    ++m_parameterRevision;
    invalidate(Dirty::PreparedMaterials | Dirty::Operator | Dirty::RightHandSide);
}

void AnalysisSession::setTime(double time)
{
    if (!std::isfinite(time))
        throw std::invalid_argument("time must be finite");
    m_parameters.time = time;
    ++m_parameterRevision;
    invalidate(Dirty::SolveState);
}

void AnalysisSession::setInitialState(std::shared_ptr<const AcceptedState> state)
{
    if (state && state->modelRevision != m_modelRevision)
        throw std::invalid_argument("initial solution belongs to another model revision");
    m_parameters.initialState = std::move(state);
    ++m_parameterRevision;
    invalidate(Dirty::InitialState);
}

void AnalysisSession::setMaterialProperty(MaterialId id, MaterialProperty property,
                                          const MaterialPropertyValue &value)
{
    if (id.value >= m_model.m_problem->blockproplist.size())
        throw std::out_of_range("invalid material id");
    auto *material = dynamic_cast<CMMaterialProp *>(m_model.m_problem->blockproplist[id.value].get());
    if (!material)
        throw std::invalid_argument("material id does not identify a magnetic material");

    switch (property) {
    case MaterialProperty::RelativePermeabilityX: material->mu_x = as<double>(value, "mu_x"); break;
    case MaterialProperty::RelativePermeabilityY: material->mu_y = as<double>(value, "mu_y"); break;
    case MaterialProperty::Coercivity: material->H_c = as<double>(value, "coercivity"); break;
    case MaterialProperty::Conductivity: material->Cduct = as<double>(value, "conductivity"); break;
    case MaterialProperty::LaminationType: material->LamType = as<int>(value, "lamination type"); break;
    case MaterialProperty::LaminationFill: material->LamFill = as<double>(value, "lamination fill"); break;
    case MaterialProperty::LaminationThickness: material->Lam_d = as<double>(value, "lamination thickness"); break;
    case MaterialProperty::AppliedCurrentDensity: material->J = as<CComplex>(value, "current density"); break;
    case MaterialProperty::BHCurve: {
        const auto &curve = as<BHCurve>(value, "B-H curve");
        if (curve.fluxDensity.size() != curve.fieldStrength.size())
            throw std::invalid_argument("B-H curve arrays must have equal length");
        material->Bdata = curve.fluxDensity;
        material->Hdata = curve.fieldStrength;
        material->BHpoints = static_cast<int>(curve.fluxDensity.size());
        material->clearSlopes();
        break;
    }
    }
    ++m_modelRevision;
    invalidate(Dirty::PreparedMaterials | Dirty::Operator);
}

void AnalysisSession::updateSolveParameters(const std::function<void(SolveParameters &)> &update)
{
    SolveParameters candidate = m_parameters;
    update(candidate);
    if (!std::isfinite(candidate.frequency) || candidate.frequency < 0 || !std::isfinite(candidate.time))
        throw std::invalid_argument("invalid frequency or time in solve parameter update");
    for (const auto &entry : candidate.circuitConstraints)
        requireCircuit(entry.first);
    for (const auto &entry : candidate.airGapPositions) {
        requireAirGap(entry.first);
        if (!std::isfinite(entry.second.innerAngle) || !std::isfinite(entry.second.outerAngle))
            throw std::invalid_argument("air-gap angles must be finite");
    }
    if (candidate.initialState && candidate.initialState->modelRevision != m_modelRevision)
        throw std::invalid_argument("initial solution belongs to another model revision");

    Dirty changed = Dirty::SolveState;
    if (candidate.frequency != m_parameters.frequency)
        changed = changed | Dirty::PreparedMaterials | Dirty::Operator | Dirty::RightHandSide;
    if (!sameConstraints(candidate.circuitConstraints, m_parameters.circuitConstraints))
        changed = changed | Dirty::RightHandSide | Dirty::Operator;
    if (!sameAirGaps(candidate.airGapPositions, m_parameters.airGapPositions))
        changed = changed | Dirty::AirGapCoupling | Dirty::Operator;
    if (candidate.initialState != m_parameters.initialState)
        changed = changed | Dirty::InitialState;

    m_parameters = std::move(candidate);
    ++m_parameterRevision;
    invalidate(changed);
}

void AnalysisSession::rebuildMaterials(PreparedAnalysis &candidate) const
{
    candidate.materials.clear();
    candidate.materials.reserve(m_model.problem().blockproplist.size());
    for (const auto &property : m_model.problem().blockproplist) {
        const auto *material = dynamic_cast<const CMMaterialProp *>(property.get());
        if (!material)
            throw std::invalid_argument("magnetic model contains a non-magnetic material");
        candidate.materials.emplace_back(*material);
        auto &prepared = candidate.materials.back();
        prepared.clearSlopes();
        if (prepared.BHpoints > 0)
            prepared.GetSlopes(m_parameters.frequency * 2.0 * 3.14159265358979323846);
    }
    candidate.frequency = m_parameters.frequency;
}

void AnalysisSession::rebuildCircuits(PreparedAnalysis &candidate) const
{
    candidate.circuits.clear();
    for (std::size_t labelIndex = 0; labelIndex < m_model.problem().labellist.size(); ++labelIndex) {
        const auto *label = dynamic_cast<const CMBlockLabel *>(m_model.problem().labellist[labelIndex].get());
        if (!label || label->InCircuit < 0)
            continue;
        CircuitId id{static_cast<std::size_t>(label->InCircuit)};
        requireCircuit(id);
        const auto constraint = m_parameters.circuitConstraints.at(id);
        candidate.circuits.push_back({id, labelIndex, label->Turns, constraint});
    }
}

void AnalysisSession::synchronize()
{
    if (m_dirty == Dirty::None || m_dirty == Dirty::SolveState)
        return;

    const Dirty requested = m_dirty;
    const auto immutableMesh = ensureMesh();
    PreparedAnalysis candidate = m_prepared;
    if (has(m_dirty, Dirty::PreparedMaterials))
        rebuildMaterials(candidate);
    if (m_instanced) {
        rebuildInstancedPrepared(candidate);
    } else if (has(m_dirty, Dirty::PreparedCircuits) || has(m_dirty, Dirty::RightHandSide)) {
        rebuildCircuits(candidate);
    }
    if (has(m_dirty, Dirty::AirGapCoupling))
        candidate.airGapPositions = m_parameters.airGapPositions;

    const Dirty rebuilt = static_cast<Dirty>(bits(requested) & ~bits(Dirty::SolveState));
    m_backend->synchronize(m_model, m_parameters, candidate, immutableMesh,
                           m_meshTopologyIdentity, rebuilt);
    m_prepared = std::move(candidate);
    m_dirty = Dirty::SolveState;
}

TrialSolution AnalysisSession::solve()
{
    synchronize();
    TrialSolution result = m_backend->solve(m_model, m_parameters, m_prepared);
    result.id = m_nextSolutionId++;
    result.sessionId = m_sessionId;
    result.modelRevision = m_modelRevision;
    result.parameterRevision = m_parameterRevision;
    result.time = m_parameters.time;
    m_dirty = Dirty::None;
    return result;
}

std::shared_ptr<const AcceptedState> AnalysisSession::acceptSolution(const TrialSolution &solution)
{
    if (solution.sessionId != m_sessionId || solution.modelRevision != m_modelRevision ||
        solution.parameterRevision != m_parameterRevision || solution.id == 0)
        throw std::invalid_argument("trial solution is stale or does not belong to this session");
    auto accepted = std::make_shared<AcceptedState>(AcceptedState{solution.modelRevision,
                                                                  solution.parameterRevision,
                                                                  solution.time,
                                                                  solution.backendState});
    m_parameters.initialState = accepted;
    ++m_parameterRevision;
    invalidate(Dirty::InitialState);
    return accepted;
}

} // namespace femm
