#include "InstancedMesh.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <map>
#include <set>
#include <tuple>
#include <utility>

namespace femm {
namespace mesh {
namespace {

constexpr double Pi = 3.141592653589793238462643383279502884;

void addDiagnostic(std::vector<MaterializationDiagnostic> &diagnostics,
                   MaterializationDiagnosticCategory category,
                   std::size_t instanceIndex, std::size_t objectIndex,
                   std::size_t localIndex, std::string message)
{
    diagnostics.push_back(
        {category, instanceIndex, objectIndex, localIndex, std::move(message)});
}

std::pair<MeshIndex, MeshIndex> edgeKey(MeshIndex first, MeshIndex second)
{
    return first < second ? std::make_pair(first, second)
                          : std::make_pair(second, first);
}

double wrapDegrees(double angle)
{
    double wrapped = std::fmod(angle, 360.0);
    if (wrapped < 0.0)
        wrapped += 360.0;
    return wrapped;
}

/** Deterministic union-find used to weld declared seam nodes. */
class UnionFind {
public:
    explicit UnionFind(std::size_t count) : m_parent(count), m_rank(count, 0)
    {
        for (std::size_t i = 0; i < count; ++i)
            m_parent[i] = i;
    }

    std::size_t find(std::size_t value)
    {
        while (m_parent[value] != value) {
            m_parent[value] = m_parent[m_parent[value]];
            value = m_parent[value];
        }
        return value;
    }

    void unite(std::size_t first, std::size_t second)
    {
        first = find(first);
        second = find(second);
        if (first == second)
            return;
        if (m_rank[first] < m_rank[second])
            std::swap(first, second);
        m_parent[second] = first;
        if (m_rank[first] == m_rank[second])
            ++m_rank[first];
    }

private:
    std::vector<std::size_t> m_parent;
    std::vector<std::size_t> m_rank;
};

const TemplateSeam *seamAt(const InstancedMesh &instanced, std::size_t instance,
                           std::size_t seam)
{
    const auto &placed = instanced.instances[instance];
    if (placed.templateIndex >= instanced.templates.size())
        return nullptr;
    const auto &meshTemplate = instanced.templates[placed.templateIndex];
    if (seam >= meshTemplate.seams.size())
        return nullptr;
    return &meshTemplate.seams[seam];
}

/** FNV-1a over a canonical byte encoding of value data only. */
class StableHash {
public:
    void addByte(unsigned char value)
    {
        m_value ^= value;
        m_value *= 1099511628211ull;
    }

    void addUint64(std::uint64_t value)
    {
        for (int i = 0; i < 8; ++i)
            addByte(static_cast<unsigned char>((value >> (8 * i)) & 0xffu));
    }

    void addSize(std::size_t value)
    {
        addUint64(static_cast<std::uint64_t>(value));
    }

    void addDouble(double value)
    {
        std::uint64_t bits = 0;
        std::memcpy(&bits, &value, sizeof(bits));
        addUint64(bits);
    }

    void addString(const std::string &value)
    {
        addSize(value.size());
        for (unsigned char byte : value)
            addByte(byte);
    }

    void addPeriodicity(SolverMesh::Periodicity value)
    {
        addByte(value == SolverMesh::Periodicity::Antiperiodic ? 1u : 0u);
    }

    std::uint64_t value() const { return m_value; }

private:
    std::uint64_t m_value = 1469598103934665603ull;
};

void hashSolverMesh(StableHash &hash, const SolverMesh &mesh)
{
    hash.addSize(mesh.nodes.size());
    for (const auto &node : mesh.nodes) {
        hash.addDouble(node.x);
        hash.addDouble(node.y);
        hash.addUint64(static_cast<std::uint64_t>(node.boundaryMarker));
    }
    hash.addSize(mesh.elements.size());
    for (const auto &element : mesh.elements) {
        for (MeshIndex node : element.nodes)
            hash.addSize(node);
        hash.addUint64(static_cast<std::uint64_t>(element.regionAttribute));
    }
    hash.addSize(mesh.edges.size());
    for (const auto &edge : mesh.edges) {
        hash.addSize(edge.first);
        hash.addSize(edge.second);
        hash.addUint64(static_cast<std::uint64_t>(edge.boundaryMarker));
    }
    hash.addSize(mesh.periodicConstraints.size());
    for (const auto &constraint : mesh.periodicConstraints) {
        hash.addSize(constraint.first);
        hash.addSize(constraint.second);
        hash.addPeriodicity(constraint.periodicity);
    }
    hash.addSize(mesh.airGaps.size());
    for (const auto &gap : mesh.airGaps) {
        hash.addString(gap.boundaryName);
        hash.addPeriodicity(gap.periodicity);
        hash.addSize(gap.totalArcElements);
        hash.addDouble(gap.totalArcLengthDegrees);
        hash.addDouble(gap.innerRadius);
        hash.addDouble(gap.outerRadius);
        hash.addDouble(gap.innerAngleDegrees);
        hash.addDouble(gap.outerAngleDegrees);
        hash.addDouble(gap.innerShift);
        hash.addDouble(gap.outerShift);
        hash.addDouble(gap.centerX);
        hash.addDouble(gap.centerY);
        hash.addSize(gap.quadraturePoints.size());
        for (const auto &point : gap.quadraturePoints) {
            for (MeshIndex node : point.nodes)
                hash.addSize(node);
            for (double weight : point.weights)
                hash.addDouble(weight);
        }
        const std::vector<SolverMesh::AirGapRingPoint> *rings[] = {
            &gap.innerRing, &gap.outerRing};
        for (const auto *ring : rings) {
            hash.addSize(ring->size());
            for (const auto &point : *ring) {
                hash.addSize(point.node);
                hash.addDouble(point.elementPosition);
                hash.addDouble(point.weight);
            }
        }
        hash.addSize(gap.nodeIndices.size());
        for (MeshIndex node : gap.nodeIndices)
            hash.addSize(node);
    }
}

} // namespace

bool RigidTransform2D::isFinite() const
{
    return std::isfinite(rotationDegrees) && std::isfinite(translationXMetres) &&
           std::isfinite(translationYMetres);
}

void RigidTransform2D::applyPoint(double x, double y, double &outX, double &outY) const
{
    const double radians = rotationDegrees * Pi / 180.0;
    const double cosine = std::cos(radians);
    const double sine = std::sin(radians);
    outX = cosine * x - sine * y + translationXMetres;
    outY = sine * x + cosine * y + translationYMetres;
}

void RigidTransform2D::applyDirection(double x, double y, double &outX, double &outY) const
{
    const double radians = rotationDegrees * Pi / 180.0;
    const double cosine = std::cos(radians);
    const double sine = std::sin(radians);
    outX = cosine * x - sine * y;
    outY = sine * x + cosine * y;
}

RigidTransform2D RigidTransform2D::rotationAbout(double centerXMetres, double centerYMetres,
                                                 double angleDegrees)
{
    RigidTransform2D transform;
    transform.rotationDegrees = angleDegrees;
    const double radians = angleDegrees * Pi / 180.0;
    const double cosine = std::cos(radians);
    const double sine = std::sin(radians);
    transform.translationXMetres = centerXMetres - (cosine * centerXMetres - sine * centerYMetres);
    transform.translationYMetres = centerYMetres - (sine * centerXMetres + cosine * centerYMetres);
    return transform;
}

std::vector<MaterializationDiagnostic> InstancedMesh::validate() const
{
    std::vector<MaterializationDiagnostic> diagnostics;

    for (std::size_t t = 0; t < templates.size(); ++t) {
        const auto &meshTemplate = templates[t];
        const auto &localMesh = meshTemplate.localMesh;
        for (std::size_t e = 0; e < localMesh.elements.size(); ++e)
            for (std::size_t k = 0; k < localMesh.elements[e].nodes.size(); ++k)
                if (localMesh.elements[e].nodes[k] >= localMesh.nodes.size())
                    addDiagnostic(diagnostics,
                                  MaterializationDiagnosticCategory::InvalidElementNode, t, e,
                                  k, "element references an invalid local node");
        for (std::size_t e = 0; e < localMesh.edges.size(); ++e) {
            if (localMesh.edges[e].first >= localMesh.nodes.size())
                addDiagnostic(diagnostics, MaterializationDiagnosticCategory::InvalidEdgeNode,
                              t, e, 0, "edge references an invalid local node");
            if (localMesh.edges[e].second >= localMesh.nodes.size())
                addDiagnostic(diagnostics, MaterializationDiagnosticCategory::InvalidEdgeNode,
                              t, e, 1, "edge references an invalid local node");
        }
        for (std::size_t s = 0; s < meshTemplate.seams.size(); ++s) {
            const auto &seam = meshTemplate.seams[s];
            if (seam.orderedNodes.empty()) {
                addDiagnostic(diagnostics, MaterializationDiagnosticCategory::EmptySeam,
                              t, s, 0, "seam '" + seam.name + "' has no nodes");
                continue;
            }
            std::set<MeshIndex> unique;
            for (std::size_t i = 0; i < seam.orderedNodes.size(); ++i) {
                const MeshIndex node = seam.orderedNodes[i];
                if (node >= meshTemplate.localMesh.nodes.size()) {
                    addDiagnostic(diagnostics,
                                  MaterializationDiagnosticCategory::InvalidSeamNode,
                                  t, s, i, "seam '" + seam.name + "' references an invalid node");
                } else if (!unique.insert(node).second) {
                    addDiagnostic(diagnostics,
                                  MaterializationDiagnosticCategory::DuplicateSeamNode,
                                  t, s, i, "seam '" + seam.name + "' repeats a node");
                }
            }
            if (seam.orderedNodes.size() < 2)
                continue;
            std::set<std::pair<MeshIndex, MeshIndex>> localEdges;
            for (const auto &edge : meshTemplate.localMesh.edges)
                localEdges.insert(edgeKey(edge.first, edge.second));
            for (std::size_t i = 1; i < seam.orderedNodes.size(); ++i) {
                const auto key = edgeKey(seam.orderedNodes[i - 1], seam.orderedNodes[i]);
                if (localEdges.count(key) == 0) {
                    addDiagnostic(diagnostics,
                                  MaterializationDiagnosticCategory::DisconnectedSeamChain,
                                  t, s, i, "seam '" + seam.name + "' is not a connected chain");
                    break;
                }
            }
        }
    }

    for (std::size_t c = 0; c < airGapCouplings.size(); ++c) {
        const auto &coupling = airGapCouplings[c];
        if (coupling.innerTemplate >= templates.size() ||
            coupling.outerTemplate >= templates.size()) {
            addDiagnostic(diagnostics, MaterializationDiagnosticCategory::InvalidAirGapCoupling,
                          c, 0, 0, "air-gap coupling references an invalid template");
            continue;
        }
        if (coupling.innerSeam >= templates[coupling.innerTemplate].seams.size() ||
            coupling.outerSeam >= templates[coupling.outerTemplate].seams.size()) {
            addDiagnostic(diagnostics, MaterializationDiagnosticCategory::InvalidAirGapCoupling,
                          c, 0, 0, "air-gap coupling references an invalid seam");
            continue;
        }
        if (coupling.innerTemplate == coupling.outerTemplate) {
            addDiagnostic(diagnostics, MaterializationDiagnosticCategory::InvalidAirGapCoupling,
                          c, 0, 0, "air-gap coupling must join two distinct templates");
        }
    }

    for (std::size_t c = 0; c < periodicClosures.size(); ++c) {
        const auto &closure = periodicClosures[c];
        if (closure.firstInstance >= instances.size() ||
            closure.secondInstance >= instances.size()) {
            addDiagnostic(diagnostics, MaterializationDiagnosticCategory::InvalidPeriodicClosure,
                          c, 0, 0, "periodic closure references an invalid instance");
            continue;
        }
        if (closure.firstSeam >=
                templates[instances[closure.firstInstance].templateIndex].seams.size() ||
            closure.secondSeam >=
                templates[instances[closure.secondInstance].templateIndex].seams.size()) {
            addDiagnostic(diagnostics, MaterializationDiagnosticCategory::InvalidPeriodicClosure,
                          c, 0, 0, "periodic closure references an invalid seam");
            continue;
        }
        if (closure.firstInstance == closure.secondInstance &&
            closure.firstSeam == closure.secondSeam) {
            addDiagnostic(diagnostics, MaterializationDiagnosticCategory::InvalidPeriodicClosure,
                          c, 0, 0, "periodic closure links a seam to itself");
            continue;
        }
        const auto &first =
            templates[instances[closure.firstInstance].templateIndex].seams[closure.firstSeam];
        const auto &second =
            templates[instances[closure.secondInstance].templateIndex].seams[closure.secondSeam];
        if (first.orderedNodes.size() != second.orderedNodes.size()) {
            addDiagnostic(diagnostics, MaterializationDiagnosticCategory::InvalidPeriodicClosure,
                          c, 0, 0, "periodic closure seams have different node counts");
        }
    }

    using SeamKey = std::pair<std::pair<std::size_t, std::size_t>,
                              std::pair<std::size_t, std::size_t>>;
    std::map<SeamKey, int> orientations;
    for (std::size_t i = 0; i < instances.size(); ++i) {
        const auto &instance = instances[i];
        if (instance.templateIndex >= templates.size()) {
            addDiagnostic(diagnostics, MaterializationDiagnosticCategory::InvalidTemplateIndex,
                          i, instance.templateIndex, 0, "instance references an invalid template");
            continue;
        }
        if (!instance.transform.isFinite()) {
            addDiagnostic(diagnostics, MaterializationDiagnosticCategory::NonFiniteTransform,
                          i, 0, 0, "instance transform is not finite");
        }
        for (std::size_t c = 0; c < instance.seamConnections.size(); ++c) {
            const auto &connection = instance.seamConnections[c];
            if (connection.seam >= templates[instance.templateIndex].seams.size()) {
                addDiagnostic(diagnostics, MaterializationDiagnosticCategory::InvalidInstanceSeam,
                              i, c, connection.seam, "connection references an invalid seam");
                continue;
            }
            if (connection.otherInstance >= instances.size()) {
                addDiagnostic(diagnostics,
                              MaterializationDiagnosticCategory::InvalidConnectionInstance,
                              i, c, connection.otherInstance,
                              "connection references an invalid instance");
                continue;
            }
            if (connection.otherInstance == i && connection.seam == connection.otherSeam) {
                addDiagnostic(diagnostics, MaterializationDiagnosticCategory::SelfConnection,
                              i, c, 0, "connection welds a seam to itself");
                continue;
            }
            const TemplateSeam *other =
                seamAt(*this, connection.otherInstance, connection.otherSeam);
            if (!other) {
                addDiagnostic(diagnostics, MaterializationDiagnosticCategory::InvalidConnectionSeam,
                              i, c, connection.otherSeam,
                              "connection references an invalid other seam");
                continue;
            }
            const auto &seam = templates[instance.templateIndex].seams[connection.seam];
            if (seam.orderedNodes.size() != other->orderedNodes.size()) {
                addDiagnostic(diagnostics, MaterializationDiagnosticCategory::SeamCardinalityMismatch,
                              i, c, 0, "connected seams have different node counts");
                continue;
            }
            const auto first = std::make_pair(i, connection.seam);
            const auto second = std::make_pair(connection.otherInstance, connection.otherSeam);
            const SeamKey key = first <= second ? SeamKey(first, second)
                                                : SeamKey(second, first);
            const int orientation =
                connection.orientation == SeamOrientation::Forward ? 0 : 1;
            const auto existing = orientations.find(key);
            if (existing == orientations.end()) {
                orientations.emplace(key, orientation);
            } else if (existing->second != orientation) {
                addDiagnostic(diagnostics,
                              MaterializationDiagnosticCategory::ConflictingConnectionOrientation,
                              i, c, 0, "seam pair is declared with conflicting orientations");
            }
        }
    }
    return diagnostics;
}

MaterializationResult InstancedMesh::materialize() const
{
    MaterializationResult result;
    result.diagnostics = validate();
    if (!result.diagnostics.empty())
        return result;

    // Flatten (instance, local node) into one key space for welding.
    std::vector<std::size_t> instanceBase(instances.size(), 0);
    std::size_t totalNodes = 0;
    for (std::size_t i = 0; i < instances.size(); ++i) {
        instanceBase[i] = totalNodes;
        totalNodes += templates[instances[i].templateIndex].localMesh.nodes.size();
    }
    UnionFind welds(totalNodes);

    // Apply each unique seam connection once. Orientation is symmetric under
    // swapping the two endpoints, so the canonical key ignores declaration order.
    using SeamKey = std::pair<std::pair<std::size_t, std::size_t>,
                              std::pair<std::size_t, std::size_t>>;
    std::map<SeamKey, SeamOrientation> applied;
    for (std::size_t i = 0; i < instances.size(); ++i) {
        const auto &instance = instances[i];
        for (const auto &connection : instance.seamConnections) {
            const auto first = std::make_pair(i, connection.seam);
            const auto second =
                std::make_pair(connection.otherInstance, connection.otherSeam);
            const SeamKey key = first <= second ? SeamKey(first, second)
                                                : SeamKey(second, first);
            applied.emplace(key, connection.orientation);
        }
    }
    for (const auto &entry : applied) {
        const std::size_t instanceA = entry.first.first.first;
        const std::size_t seamA = entry.first.first.second;
        const std::size_t instanceB = entry.first.second.first;
        const std::size_t seamB = entry.first.second.second;
        const TemplateSeam &first = templates[instances[instanceA].templateIndex]
                                        .seams[seamA];
        const TemplateSeam &second = templates[instances[instanceB].templateIndex]
                                         .seams[seamB];
        const bool reverse = entry.second == SeamOrientation::Reverse;
        const std::size_t count = first.orderedNodes.size();
        for (std::size_t i = 0; i < count; ++i) {
            const std::size_t j = reverse ? count - 1 - i : i;
            welds.unite(instanceBase[instanceA] + first.orderedNodes[i],
                        instanceBase[instanceB] + second.orderedNodes[j]);
        }
    }

    // Validate that every declared weld actually places its nodes together.
    for (std::size_t i = 0; i < instances.size(); ++i) {
        const auto &instance = instances[i];
        const auto &transform = instance.transform;
        for (const auto &connection : instance.seamConnections) {
            const TemplateSeam &first =
                templates[instance.templateIndex].seams[connection.seam];
            const TemplateSeam &second =
                templates[instances[connection.otherInstance].templateIndex]
                    .seams[connection.otherSeam];
            const bool reverse = connection.orientation == SeamOrientation::Reverse;
            const std::size_t count = first.orderedNodes.size();
            const auto &firstMesh =
                templates[instance.templateIndex].localMesh;
            const auto &secondMesh =
                templates[instances[connection.otherInstance].templateIndex].localMesh;
            const auto &secondTransform = instances[connection.otherInstance].transform;
            for (std::size_t k = 0; k < count; ++k) {
                const std::size_t j = reverse ? count - 1 - k : k;
                const auto &a = firstMesh.nodes[first.orderedNodes[k]];
                const auto &b = secondMesh.nodes[second.orderedNodes[j]];
                double ax, ay, bx, by;
                transform.applyPoint(a.x, a.y, ax, ay);
                secondTransform.applyPoint(b.x, b.y, bx, by);
                if (std::abs(ax - bx) > InstancedMeshWeldToleranceMetres ||
                    std::abs(ay - by) > InstancedMeshWeldToleranceMetres) {
                    addDiagnostic(result.diagnostics,
                                  MaterializationDiagnosticCategory::SeamCoordinateMismatch,
                                  i, connection.seam, k,
                                  "declared weld nodes do not coincide");
                }
            }
        }
    }
    if (!result.diagnostics.empty())
        return result;

    // Assign deterministic global node indices: first encounter wins.
    std::vector<MeshIndex> rootToGlobal(totalNodes, InvalidMeshIndex);
    result.provenance.nodeMap.resize(instances.size());
    for (std::size_t i = 0; i < instances.size(); ++i) {
        const auto &instance = instances[i];
        const auto &localMesh = templates[instance.templateIndex].localMesh;
        result.provenance.nodeMap[i].resize(localMesh.nodes.size(), InvalidMeshIndex);
        for (MeshIndex n = 0; n < localMesh.nodes.size(); ++n) {
            const std::size_t root = welds.find(instanceBase[i] + n);
            if (rootToGlobal[root] == InvalidMeshIndex) {
                if (result.mesh.nodes.size() >= InvalidMeshIndex) {
                    addDiagnostic(result.diagnostics,
                                  MaterializationDiagnosticCategory::IndexOverflow, i, n, 0,
                                  "global node index overflow");
                    return result;
                }
                const auto &node = localMesh.nodes[n];
                double x = 0.0, y = 0.0;
                instance.transform.applyPoint(node.x, node.y, x, y);
                rootToGlobal[root] = result.mesh.nodes.size();
                result.mesh.nodes.push_back({x, y, node.boundaryMarker});
                result.provenance.nodeProvenance.push_back(
                    {instance.templateIndex, i, n});
            }
            result.provenance.nodeMap[i][n] = rootToGlobal[root];
        }
    }

    // Materialise elements, rejecting degenerate or reversed triangles.
    result.provenance.elementMap.resize(instances.size());
    for (std::size_t i = 0; i < instances.size(); ++i) {
        const auto &instance = instances[i];
        const auto &localMesh = templates[instance.templateIndex].localMesh;
        result.provenance.elementMap[i].resize(localMesh.elements.size(), InvalidMeshIndex);
        for (MeshIndex e = 0; e < localMesh.elements.size(); ++e) {
            const auto &element = localMesh.elements[e];
            std::array<MeshIndex, 3> global{};
            bool valid = true;
            for (std::size_t k = 0; k < element.nodes.size(); ++k) {
                if (element.nodes[k] >= localMesh.nodes.size()) {
                    addDiagnostic(result.diagnostics,
                                  MaterializationDiagnosticCategory::InvalidElementNode, i, e,
                                  k, "element references an invalid local node");
                    valid = false;
                    continue;
                }
                global[k] = result.provenance.nodeMap[i][element.nodes[k]];
            }
            if (!valid)
                continue;
            const auto &a = result.mesh.nodes[global[0]];
            const auto &b = result.mesh.nodes[global[1]];
            const auto &c = result.mesh.nodes[global[2]];
            const double twiceArea =
                (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
            if (!std::isfinite(twiceArea) || twiceArea == 0.0) {
                addDiagnostic(result.diagnostics,
                              MaterializationDiagnosticCategory::DegenerateElement, i, e, 0,
                              "materialised element is degenerate");
                continue;
            }
            if (twiceArea < 0.0) {
                addDiagnostic(result.diagnostics,
                              MaterializationDiagnosticCategory::ReversedElement, i, e, 0,
                              "materialised element has reversed orientation");
                continue;
            }
            result.provenance.elementMap[i][e] = result.mesh.elements.size();
            result.mesh.elements.push_back({global, element.regionAttribute});
            result.provenance.elementProvenance.push_back({instance.templateIndex, i, e});
        }
    }

    // Materialise edges, canonicalising duplicates and rejecting marker conflicts.
    std::map<std::pair<MeshIndex, MeshIndex>, std::int32_t> edgeMarkers;
    for (std::size_t i = 0; i < instances.size(); ++i) {
        const auto &localMesh = templates[instances[i].templateIndex].localMesh;
        for (std::size_t e = 0; e < localMesh.edges.size(); ++e) {
            const auto &edge = localMesh.edges[e];
            if (edge.first >= localMesh.nodes.size() ||
                edge.second >= localMesh.nodes.size()) {
                addDiagnostic(result.diagnostics,
                              MaterializationDiagnosticCategory::InvalidEdgeNode, i, e, 0,
                              "edge references an invalid local node");
                continue;
            }
            const auto key = edgeKey(result.provenance.nodeMap[i][edge.first],
                                     result.provenance.nodeMap[i][edge.second]);
            const auto existing = edgeMarkers.find(key);
            if (existing == edgeMarkers.end()) {
                edgeMarkers.emplace(key, edge.boundaryMarker);
            } else if (existing->second != edge.boundaryMarker) {
                addDiagnostic(result.diagnostics,
                              MaterializationDiagnosticCategory::EdgeMarkerConflict, i, e, 0,
                              "a duplicated edge has conflicting boundary markers");
            }
        }
    }
    for (const auto &entry : edgeMarkers)
        result.mesh.edges.push_back({entry.first.first, entry.first.second, entry.second});

    // Remap ordinary periodic constraints.
    std::set<std::tuple<MeshIndex, MeshIndex, int>> seenConstraints;
    for (std::size_t i = 0; i < instances.size(); ++i) {
        const auto &localMesh = templates[instances[i].templateIndex].localMesh;
        for (std::size_t c = 0; c < localMesh.periodicConstraints.size(); ++c) {
            const auto &constraint = localMesh.periodicConstraints[c];
            if (constraint.first >= localMesh.nodes.size() ||
                constraint.second >= localMesh.nodes.size()) {
                addDiagnostic(result.diagnostics,
                              MaterializationDiagnosticCategory::InvalidPeriodicNode, i, c, 0,
                              "periodic constraint references an invalid local node");
                continue;
            }
            const MeshIndex first = result.provenance.nodeMap[i][constraint.first];
            const MeshIndex second = result.provenance.nodeMap[i][constraint.second];
            const int type =
                constraint.periodicity == SolverMesh::Periodicity::Antiperiodic ? 1 : 0;
            if (seenConstraints.insert({std::min(first, second), std::max(first, second), type})
                    .second) {
                result.mesh.periodicConstraints.push_back(
                    {first, second, constraint.periodicity});
            }
        }
    }

    // Emit the periodic/antiperiodic links that close an open sector. The two
    // end seams are not welded, so their nodes stay distinct and the field is
    // identified only by the explicit constraint.
    for (const auto &closure : periodicClosures) {
        const auto &first =
            templates[instances[closure.firstInstance].templateIndex].seams[closure.firstSeam];
        const auto &second =
            templates[instances[closure.secondInstance].templateIndex].seams[closure.secondSeam];
        const bool reverse = closure.orientation == SeamOrientation::Reverse;
        const std::size_t count = first.orderedNodes.size();
        for (std::size_t i = 0; i < count; ++i) {
            const std::size_t j = reverse ? count - 1 - i : i;
            const MeshIndex a =
                result.provenance.nodeMap[closure.firstInstance][first.orderedNodes[i]];
            const MeshIndex b =
                result.provenance.nodeMap[closure.secondInstance][second.orderedNodes[j]];
            const int type =
                closure.periodicity == SolverMesh::Periodicity::Antiperiodic ? 1 : 0;
            if (seenConstraints.insert({std::min(a, b), std::max(a, b), type}).second) {
                result.mesh.periodicConstraints.push_back({a, b, closure.periodicity});
            }
        }
    }

    // Remap AGE rings, quadrature nodes, and node-index lists.
    for (std::size_t i = 0; i < instances.size(); ++i) {
        const auto &instance = instances[i];
        const auto &localMesh = templates[instance.templateIndex].localMesh;
        for (std::size_t g = 0; g < localMesh.airGaps.size(); ++g) {
            const auto &gap = localMesh.airGaps[g];
            SolverMesh::AirGap mapped;
            mapped.boundaryName = gap.boundaryName;
            mapped.periodicity = gap.periodicity;
            mapped.totalArcElements = gap.totalArcElements;
            mapped.totalArcLengthDegrees = gap.totalArcLengthDegrees;
            mapped.innerRadius = gap.innerRadius;
            mapped.outerRadius = gap.outerRadius;
            mapped.innerShift = gap.innerShift;
            mapped.outerShift = gap.outerShift;
            double centerX = 0.0, centerY = 0.0;
            instance.transform.applyPoint(gap.centerX, gap.centerY, centerX, centerY);
            mapped.centerX = centerX;
            mapped.centerY = centerY;
            mapped.innerAngleDegrees =
                wrapDegrees(gap.innerAngleDegrees + instance.transform.rotationDegrees);
            mapped.outerAngleDegrees =
                wrapDegrees(gap.outerAngleDegrees + instance.transform.rotationDegrees);

            bool structureValid = true;
            const auto remapNode = [&](MeshIndex node, std::size_t localIndex,
                                       const char *description) -> MeshIndex {
                if (node >= localMesh.nodes.size()) {
                    addDiagnostic(result.diagnostics,
                                  MaterializationDiagnosticCategory::InvalidAirGapNode, i, g,
                                  localIndex, std::string("AGE ") + description +
                                                  " references an invalid local node");
                    structureValid = false;
                    return InvalidMeshIndex;
                }
                return result.provenance.nodeMap[i][node];
            };

            mapped.nodeIndices.reserve(gap.nodeIndices.size());
            for (std::size_t n = 0; n < gap.nodeIndices.size(); ++n)
                mapped.nodeIndices.push_back(
                    remapNode(gap.nodeIndices[n], n, "nodeIndices"));

            mapped.quadraturePoints.reserve(gap.quadraturePoints.size());
            for (std::size_t q = 0; q < gap.quadraturePoints.size(); ++q) {
                SolverMesh::AirGapQuadraturePoint point;
                for (std::size_t k = 0; k < 4; ++k)
                    point.nodes[k] = remapNode(gap.quadraturePoints[q].nodes[k],
                                               q * 4 + k, "quadrature");
                point.weights = gap.quadraturePoints[q].weights;
                mapped.quadraturePoints.push_back(point);
            }

            const std::vector<SolverMesh::AirGapRingPoint> *sourceRings[] = {
                &gap.innerRing, &gap.outerRing};
            std::vector<SolverMesh::AirGapRingPoint> *targetRings[] = {
                &mapped.innerRing, &mapped.outerRing};
            for (std::size_t r = 0; r < 2; ++r) {
                targetRings[r]->reserve(sourceRings[r]->size());
                for (std::size_t p = 0; p < sourceRings[r]->size(); ++p) {
                    SolverMesh::AirGapRingPoint point = (*sourceRings[r])[p];
                    point.node = remapNode(point.node, p, "ring");
                    targetRings[r]->push_back(point);
                }
            }

            if (gap.quadraturePoints.empty() ||
                gap.quadraturePoints.size() - 1 != gap.totalArcElements ||
                gap.innerRing.empty() != gap.outerRing.empty() ||
                (!gap.innerRing.empty() &&
                 (gap.innerRing.size() != gap.outerRing.size() ||
                  gap.innerRing.size() < gap.totalArcElements))) {
                addDiagnostic(result.diagnostics,
                              MaterializationDiagnosticCategory::InvalidAirGapStructure, i, g,
                              0, "AGE structure is inconsistent");
                structureValid = false;
            }
            if (structureValid)
                result.mesh.airGaps.push_back(std::move(mapped));
        }
    }

    // Assemble cross-template air-gap couplings from the per-instance seam
    // nodes. The inner and outer rings keep independent unknowns; only the
    // solver's air-gap element links them.
    for (std::size_t c = 0; c < airGapCouplings.size(); ++c) {
        const auto &coupling = airGapCouplings[c];
        SolverMesh::AirGap gap;
        gap.boundaryName = coupling.boundaryName;
        gap.periodicity = coupling.periodicity;
        gap.totalArcLengthDegrees = coupling.totalArcLengthDegrees;
        gap.innerRadius = coupling.innerRadiusMetres;
        gap.outerRadius = coupling.outerRadiusMetres;
        gap.centerX = coupling.centerXMetres;
        gap.centerY = coupling.centerYMetres;
        gap.innerAngleDegrees = coupling.innerAngleDegrees;
        gap.outerAngleDegrees = coupling.outerAngleDegrees;
        gap.innerShift = coupling.innerShift;
        gap.outerShift = coupling.outerShift;

        const auto collectRing = [&](std::size_t templateIndex, std::size_t seamIndex,
                                     std::vector<SolverMesh::AirGapRingPoint> &ring) {
            const auto &seam = templates[templateIndex].seams[seamIndex];
            for (std::size_t i = 0; i < instances.size(); ++i) {
                if (instances[i].templateIndex != templateIndex)
                    continue;
                for (MeshIndex local : seam.orderedNodes) {
                    const MeshIndex global = result.provenance.nodeMap[i][local];
                    const auto &node = result.mesh.nodes[global];
                    double angle = std::atan2(node.y - coupling.centerYMetres,
                                              node.x - coupling.centerXMetres) *
                                   180.0 / Pi;
                    if (angle < 0.0)
                        angle += 360.0;
                    ring.push_back({global, angle, 1.0});
                }
            }
            std::stable_sort(ring.begin(), ring.end(),
                             [](const SolverMesh::AirGapRingPoint &left,
                                const SolverMesh::AirGapRingPoint &right) {
                                 return left.elementPosition < right.elementPosition;
                             });
            // Adjacent welded instances contribute the same boundary node
            // twice; keep one ring entry per global node.
            std::vector<SolverMesh::AirGapRingPoint> unique;
            unique.reserve(ring.size());
            for (const auto &point : ring)
                if (unique.empty() || unique.back().node != point.node)
                    unique.push_back(point);
            ring = std::move(unique);
        };
        collectRing(coupling.innerTemplate, coupling.innerSeam, gap.innerRing);
        collectRing(coupling.outerTemplate, coupling.outerSeam, gap.outerRing);

        if (gap.innerRing.empty() || gap.innerRing.size() != gap.outerRing.size()) {
            addDiagnostic(result.diagnostics,
                          MaterializationDiagnosticCategory::InvalidAirGapCoupling, c, 0, 0,
                          "air-gap coupling rings have incompatible cardinality");
            continue;
        }
        const std::size_t count = gap.innerRing.size();
        const double step = coupling.totalArcLengthDegrees / static_cast<double>(count);
        if (!(step > 0.0)) {
            addDiagnostic(result.diagnostics,
                          MaterializationDiagnosticCategory::InvalidAirGapCoupling, c, 0, 0,
                          "air-gap coupling has a non-positive angular step");
            continue;
        }
        for (auto &point : gap.innerRing)
            point.elementPosition /= step;
        for (auto &point : gap.outerRing)
            point.elementPosition /= step;
        gap.innerShift = gap.innerRing.front().elementPosition;
        gap.outerShift = gap.outerRing.front().elementPosition;
        gap.totalArcElements = count;
        gap.nodeIndices.reserve(count * 2);
        for (const auto &point : gap.innerRing)
            gap.nodeIndices.push_back(point.node);
        for (const auto &point : gap.outerRing)
            gap.nodeIndices.push_back(point.node);
        gap.quadraturePoints.reserve(count + 1);
        for (std::size_t i = 0; i <= count; ++i) {
            const std::size_t next = i % count;
            const std::size_t previous = next == 0 ? count - 1 : next - 1;
            SolverMesh::AirGapQuadraturePoint point;
            point.nodes = {{gap.innerRing[previous].node, gap.innerRing[next].node,
                            gap.outerRing[previous].node, gap.outerRing[next].node}};
            point.weights = {{gap.innerRing[previous].weight, gap.innerRing[next].weight,
                              gap.outerRing[previous].weight, gap.outerRing[next].weight}};
            gap.quadraturePoints.push_back(point);
        }
        result.mesh.airGaps.push_back(std::move(gap));
    }

    return result;
}

std::uint64_t templateTopologyIdentity(const InstancedMesh &instanced)
{
    StableHash hash;
    hash.addSize(instanced.templates.size());
    for (const auto &meshTemplate : instanced.templates) {
        hashSolverMesh(hash, meshTemplate.localMesh);
        hash.addSize(meshTemplate.seams.size());
        for (const auto &seam : meshTemplate.seams) {
            hash.addString(seam.name);
            hash.addSize(seam.orderedNodes.size());
            for (MeshIndex node : seam.orderedNodes)
                hash.addSize(node);
        }
    }
    return hash.value();
}

std::uint64_t instanceLayoutIdentity(const InstancedMesh &instanced)
{
    StableHash hash;
    hash.addSize(instanced.instances.size());
    for (const auto &instance : instanced.instances) {
        hash.addSize(instance.templateIndex);
        hash.addDouble(instance.transform.rotationDegrees);
        hash.addDouble(instance.transform.translationXMetres);
        hash.addDouble(instance.transform.translationYMetres);
        hash.addSize(instance.seamConnections.size());
        for (const auto &connection : instance.seamConnections) {
            hash.addSize(connection.seam);
            hash.addSize(connection.otherInstance);
            hash.addSize(connection.otherSeam);
            hash.addByte(connection.orientation == SeamOrientation::Reverse ? 1u : 0u);
        }
    }
    hash.addSize(instanced.periodicClosures.size());
    for (const auto &closure : instanced.periodicClosures) {
        hash.addSize(closure.firstInstance);
        hash.addSize(closure.firstSeam);
        hash.addSize(closure.secondInstance);
        hash.addSize(closure.secondSeam);
        hash.addByte(closure.orientation == SeamOrientation::Reverse ? 1u : 0u);
        hash.addByte(closure.periodicity == SolverMesh::Periodicity::Antiperiodic ? 1u : 0u);
    }
    return hash.value();
}

std::uint64_t materializedTopologyIdentity(const SolverMesh &mesh)
{
    StableHash hash;
    hashSolverMesh(hash, mesh);
    return hash.value();
}

} // namespace mesh
} // namespace femm
