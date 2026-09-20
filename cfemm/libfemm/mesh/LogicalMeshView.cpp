#include "LogicalMeshView.h"

#include "InstancedMeshDetail.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <set>

namespace femm {
namespace mesh {
namespace {

std::pair<MeshIndex, MeshIndex> edgeKey(MeshIndex first, MeshIndex second)
{
    return first < second ? std::make_pair(first, second)
                          : std::make_pair(second, first);
}

} // namespace

LogicalMeshView LogicalMeshView::build(const InstancedMesh &instanced,
                                       std::vector<MaterializationDiagnostic> &diagnostics)
{
    LogicalMeshView view;
    detail::WeldPlan plan;
    if (!detail::buildWeldPlan(instanced, plan, diagnostics)) {
        view.m_diagnostics = diagnostics;
        return view;
    }

    view.m_templates = instanced.templates;
    view.m_instances = instanced.instances;
    const std::size_t instanceCount = view.m_instances.size();

    // Per-template local edge markers let an element's boundary sides be read
    // without mapping its nodes to global indices first.
    view.m_templateEdgeMarkers.resize(view.m_templates.size());
    for (std::size_t t = 0; t < view.m_templates.size(); ++t) {
        auto &markers = view.m_templateEdgeMarkers[t];
        for (const auto &edge : view.m_templates[t].localMesh.edges)
            markers.emplace(edgeKey(edge.first, edge.second), edge.boundaryMarker);
    }

    // Assign global node indices exactly as InstancedMesh::materialize(): first
    // encounter wins over instances and then local nodes. Record every node that
    // welds onto an already-assigned representative.
    view.m_instanceNodeBase.assign(instanceCount + 1, 0);
    view.m_weldedNodes.resize(instanceCount);
    std::vector<MeshIndex> rootToGlobal(plan.totalNodes, InvalidMeshIndex);
    std::vector<NodeProvenance> representative;
    std::vector<char> weldedGlobal;
    MeshIndex nextNode = 0;
    for (std::size_t i = 0; i < instanceCount; ++i) {
        view.m_instanceNodeBase[i] = nextNode;
        const auto &localMesh = view.m_templates[view.m_instances[i].templateIndex].localMesh;
        for (MeshIndex local = 0; local < localMesh.nodes.size(); ++local) {
            const std::size_t root = plan.root[plan.instanceBase[i] + local];
            if (rootToGlobal[root] == InvalidMeshIndex) {
                rootToGlobal[root] = nextNode++;
                representative.push_back({view.m_instances[i].templateIndex, i, local});
                weldedGlobal.push_back(0);
            } else {
                const MeshIndex global = rootToGlobal[root];
                view.m_weldedNodes[i].push_back({local, global});
                weldedGlobal[global] = 1;
            }
        }
    }
    view.m_instanceNodeBase[instanceCount] = nextNode;
    view.m_nodeCount = nextNode;
    for (MeshIndex global = 0; global < nextNode; ++global) {
        if (!weldedGlobal[global])
            continue;
        const NodeProvenance &node = representative[global];
        view.m_weldedGlobalToLocal.push_back(
            {global, {node.instanceIndex, node.localNode}});
    }
    std::sort(view.m_weldedGlobalToLocal.begin(), view.m_weldedGlobalToLocal.end(),
              [](const std::pair<MeshIndex, std::pair<std::size_t, MeshIndex>> &left,
                 const std::pair<MeshIndex, std::pair<std::size_t, MeshIndex>> &right) {
                  return left.first < right.first;
              });

    view.m_instanceElementBase.assign(instanceCount + 1, 0);
    std::size_t elementTotal = 0;
    for (std::size_t i = 0; i < instanceCount; ++i) {
        view.m_instanceElementBase[i] = elementTotal;
        elementTotal +=
            view.m_templates[view.m_instances[i].templateIndex].localMesh.elements.size();
    }
    view.m_instanceElementBase[instanceCount] = elementTotal;
    view.m_elementCount = elementTotal;

    // Store only the physical boundary edges; unmarked interior edges are
    // recomputed from element connectivity when they are needed.
    std::map<std::pair<MeshIndex, MeshIndex>, std::int32_t> marked;
    for (std::size_t i = 0; i < instanceCount; ++i) {
        const auto &localMesh = view.m_templates[view.m_instances[i].templateIndex].localMesh;
        for (const auto &edge : localMesh.edges) {
            if (edge.first >= localMesh.nodes.size() || edge.second >= localMesh.nodes.size())
                continue; // already rejected by buildWeldPlan
            if (edge.boundaryMarker == 0)
                continue;
            const MeshIndex first = view.nodeFor(i, edge.first);
            const MeshIndex second = view.nodeFor(i, edge.second);
            const auto key = edgeKey(first, second);
            const auto inserted = marked.emplace(key, edge.boundaryMarker);
            if (!inserted.second && inserted.first->second != edge.boundaryMarker) {
                diagnostics.push_back({MaterializationDiagnosticCategory::EdgeMarkerConflict, i,
                                       0, 0,
                                       "a duplicated edge has conflicting boundary markers"});
                view.m_diagnostics = diagnostics;
                return view;
            }
        }
    }
    view.m_markedEdges.assign(marked.begin(), marked.end());

    detail::TopologyRemap topology;
    if (!detail::remapTopology(
            instanced,
            [&view](std::size_t instance, MeshIndex local) {
                return view.nodeFor(instance, local);
            },
            [&view](MeshIndex node, double &x, double &y) {
                view.nodeCoordinates(node, x, y);
            },
            topology, diagnostics)) {
        view.m_diagnostics = diagnostics;
        return view;
    }
    view.m_periodicConstraints = std::move(topology.periodicConstraints);
    view.m_airGaps = std::move(topology.airGaps);
    view.m_valid = true;
    return view;
}

MeshIndex LogicalMeshView::nodeFor(std::size_t instance, MeshIndex localNode) const
{
    const auto &welded = m_weldedNodes[instance];
    const auto found = std::lower_bound(
        welded.begin(), welded.end(), localNode,
        [](const std::pair<MeshIndex, MeshIndex> &entry, MeshIndex value) {
            return entry.first < value;
        });
    if (found != welded.end() && found->first == localNode)
        return found->second;
    const std::size_t weldedBefore = static_cast<std::size_t>(found - welded.begin());
    return static_cast<MeshIndex>(m_instanceNodeBase[instance] + localNode - weldedBefore);
}

MeshIndex LogicalMeshView::elementFor(std::size_t instance, MeshIndex localElement) const
{
    return static_cast<MeshIndex>(m_instanceElementBase[instance] + localElement);
}

NodeProvenance LogicalMeshView::nodeProvenance(MeshIndex node) const
{
    const auto welded = std::lower_bound(
        m_weldedGlobalToLocal.begin(), m_weldedGlobalToLocal.end(), node,
        [](const std::pair<MeshIndex, std::pair<std::size_t, MeshIndex>> &entry,
           MeshIndex value) { return entry.first < value; });
    if (welded != m_weldedGlobalToLocal.end() && welded->first == node) {
        return {m_instances[welded->second.first].templateIndex, welded->second.first,
                welded->second.second};
    }

    const auto base =
        std::upper_bound(m_instanceNodeBase.begin(), m_instanceNodeBase.end(), node);
    const std::size_t instance =
        static_cast<std::size_t>(base - m_instanceNodeBase.begin()) - 1;
    const std::size_t rank = node - m_instanceNodeBase[instance];

    const auto &weldedLocal = m_weldedNodes[instance];
    const std::size_t localCount =
        m_templates[m_instances[instance].templateIndex].localMesh.nodes.size();
    std::size_t low = 0;
    std::size_t high = localCount;
    while (low < high) {
        const std::size_t middle = low + (high - low) / 2;
        const auto found = std::lower_bound(
            weldedLocal.begin(), weldedLocal.end(), static_cast<MeshIndex>(middle),
            [](const std::pair<MeshIndex, MeshIndex> &entry, MeshIndex value) {
                return entry.first < value;
            });
        const std::size_t weldedBefore = static_cast<std::size_t>(found - weldedLocal.begin());
        if (middle - weldedBefore < rank)
            low = middle + 1;
        else
            high = middle;
    }
    // f(n) = n - (welded nodes before n) is flat across a welded node, so the
    // first candidate with f(n) >= rank can be welded; advance to the first
    // non-welded local node, which carries the requested rank.
    while (low < localCount) {
        const auto found = std::lower_bound(
            weldedLocal.begin(), weldedLocal.end(), static_cast<MeshIndex>(low),
            [](const std::pair<MeshIndex, MeshIndex> &entry, MeshIndex value) {
                return entry.first < value;
            });
        if (found == weldedLocal.end() || found->first != low)
            break;
        ++low;
    }
    return {m_instances[instance].templateIndex, instance, static_cast<MeshIndex>(low)};
}

void LogicalMeshView::nodeCoordinates(MeshIndex node, double &x, double &y) const
{
    const NodeProvenance provenance = nodeProvenance(node);
    const auto &local =
        m_templates[provenance.templateIndex].localMesh.nodes[provenance.localNode];
    m_instances[provenance.instanceIndex].transform.applyPoint(local.x, local.y, x, y);
}

void LogicalMeshView::allNodeCoordinates(std::vector<double> &x, std::vector<double> &y) const
{
    x.assign(m_nodeCount, 0.0);
    y.assign(m_nodeCount, 0.0);
    for (std::size_t i = 0; i < m_instances.size(); ++i) {
        const auto &localMesh = m_templates[m_instances[i].templateIndex].localMesh;
        const auto &transform = m_instances[i].transform;
        for (MeshIndex local = 0; local < localMesh.nodes.size(); ++local) {
            const MeshIndex global = nodeFor(i, local);
            transform.applyPoint(localMesh.nodes[local].x, localMesh.nodes[local].y, x[global],
                                 y[global]);
        }
    }
}

std::int32_t LogicalMeshView::nodeBoundaryMarker(MeshIndex node) const
{
    const NodeProvenance provenance = nodeProvenance(node);
    return m_templates[provenance.templateIndex]
        .localMesh.nodes[provenance.localNode]
        .boundaryMarker;
}

ElementProvenance LogicalMeshView::elementProvenance(MeshIndex element) const
{
    const auto base =
        std::upper_bound(m_instanceElementBase.begin(), m_instanceElementBase.end(), element);
    const std::size_t instance =
        static_cast<std::size_t>(base - m_instanceElementBase.begin()) - 1;
    return {m_instances[instance].templateIndex, instance,
            static_cast<MeshIndex>(element - m_instanceElementBase[instance])};
}

std::array<MeshIndex, 3> LogicalMeshView::elementNodes(MeshIndex element) const
{
    const ElementProvenance provenance = elementProvenance(element);
    const auto &local =
        m_templates[provenance.templateIndex].localMesh.elements[provenance.localElement];
    std::array<MeshIndex, 3> nodes{};
    for (std::size_t k = 0; k < nodes.size(); ++k)
        nodes[k] = nodeFor(provenance.instanceIndex, local.nodes[k]);
    return nodes;
}

void LogicalMeshView::elementCoordinates(MeshIndex element, double x[3], double y[3]) const
{
    const ElementProvenance provenance = elementProvenance(element);
    const auto &local =
        m_templates[provenance.templateIndex].localMesh.elements[provenance.localElement];
    const auto &transform = m_instances[provenance.instanceIndex].transform;
    for (std::size_t k = 0; k < 3; ++k) {
        const auto &node =
            m_templates[provenance.templateIndex].localMesh.nodes[local.nodes[k]];
        transform.applyPoint(node.x, node.y, x[k], y[k]);
    }
}

std::array<std::int32_t, 3> LogicalMeshView::elementEdgeMarkers(MeshIndex element) const
{
    const ElementProvenance provenance = elementProvenance(element);
    const auto &local =
        m_templates[provenance.templateIndex].localMesh.elements[provenance.localElement];
    const auto &markers = m_templateEdgeMarkers[provenance.templateIndex];
    std::array<std::int32_t, 3> result{};
    for (std::size_t k = 0; k < 3; ++k) {
        const auto found = markers.find(edgeKey(local.nodes[k], local.nodes[(k + 1) % 3]));
        result[k] = found == markers.end() ? 0 : found->second;
    }
    return result;
}

std::int32_t LogicalMeshView::elementRegionAttribute(MeshIndex element) const
{
    const ElementProvenance provenance = elementProvenance(element);
    return m_templates[provenance.templateIndex]
        .localMesh.elements[provenance.localElement]
        .regionAttribute;
}

std::size_t LogicalMeshView::elementInstance(MeshIndex element) const
{
    return elementProvenance(element).instanceIndex;
}

std::int32_t LogicalMeshView::boundaryMarkerForEdge(MeshIndex first, MeshIndex second) const
{
    const auto key = edgeKey(first, second);
    const auto found = std::lower_bound(
        m_markedEdges.begin(), m_markedEdges.end(), key,
        [](const std::pair<std::pair<MeshIndex, MeshIndex>, std::int32_t> &entry,
           const std::pair<MeshIndex, MeshIndex> &value) { return entry.first < value; });
    if (found != m_markedEdges.end() && found->first == key)
        return found->second;
    return 0;
}

std::size_t LogicalMeshView::edgeCount() const
{
    std::set<std::pair<MeshIndex, MeshIndex>> edges;
    for (MeshIndex element = 0; element < m_elementCount; ++element) {
        const auto nodes = elementNodes(element);
        for (std::size_t k = 0; k < 3; ++k)
            edges.insert(edgeKey(nodes[k], nodes[(k + 1) % 3]));
    }
    return edges.size();
}

LogicalMeshView::Adjacency LogicalMeshView::buildAdjacency() const
{
    // Two-pass compressed sparse row construction. Counting degrees first and
    // filling once avoids one heap allocation per node.
    Adjacency result;
    result.offsets.assign(m_nodeCount + 1, 0);
    for (MeshIndex element = 0; element < m_elementCount; ++element) {
        const auto nodes = elementNodes(element);
        for (std::size_t k = 0; k < 3; ++k) {
            ++result.offsets[nodes[k] + 1];
            ++result.offsets[nodes[(k + 1) % 3] + 1];
        }
    }
    for (std::size_t node = 0; node < m_nodeCount; ++node)
        result.offsets[node + 1] += result.offsets[node];
    result.neighbors.resize(result.offsets[m_nodeCount]);

    std::vector<std::size_t> cursor(result.offsets.begin(), result.offsets.end() - 1);
    for (MeshIndex element = 0; element < m_elementCount; ++element) {
        const auto nodes = elementNodes(element);
        for (std::size_t k = 0; k < 3; ++k) {
            const MeshIndex first = nodes[k];
            const MeshIndex second = nodes[(k + 1) % 3];
            result.neighbors[cursor[first]++] = second;
            result.neighbors[cursor[second]++] = first;
        }
    }

    // Sort and deduplicate each row, then compact the rows in place.
    std::vector<std::size_t> compacted(m_nodeCount + 1, 0);
    std::size_t write = 0;
    for (std::size_t node = 0; node < m_nodeCount; ++node) {
        const std::size_t begin = result.offsets[node];
        const std::size_t end = result.offsets[node + 1];
        std::sort(result.neighbors.begin() + begin, result.neighbors.begin() + end);
        std::size_t uniqueEnd = begin;
        for (std::size_t i = begin; i < end; ++i)
            if (i == begin || result.neighbors[i] != result.neighbors[i - 1])
                result.neighbors[uniqueEnd++] = result.neighbors[i];
        for (std::size_t i = begin; i < uniqueEnd; ++i)
            result.neighbors[write++] = result.neighbors[i];
        compacted[node + 1] = write;
    }
    result.neighbors.resize(write);
    result.offsets = std::move(compacted);
    return result;
}

std::vector<MeshIndex> LogicalMeshView::cuthillMcKeeOrdering() const
{
    return cuthillMcKeeOrdering(buildAdjacency());
}

std::vector<MeshIndex> LogicalMeshView::cuthillMcKeeOrdering(const Adjacency &adjacency) const
{
    std::vector<MeshIndex> order;
    order.reserve(m_nodeCount);
    std::vector<char> visited(m_nodeCount, 0);
    std::vector<MeshIndex> frontier;

    while (order.size() < m_nodeCount) {
        // Start each component at the minimum-degree unvisited node, breaking
        // ties by the smallest global index for determinism.
        std::size_t start = m_nodeCount;
        std::size_t startDegree = 0;
        for (std::size_t node = 0; node < m_nodeCount; ++node) {
            if (visited[node])
                continue;
            const std::size_t degree = adjacency.degree(node);
            if (start == m_nodeCount || degree < startDegree) {
                start = node;
                startDegree = degree;
            }
        }
        visited[start] = 1;
        order.push_back(static_cast<MeshIndex>(start));
        frontier.clear();
        frontier.push_back(static_cast<MeshIndex>(start));

        // Index-based queue: erasing from the front would make the traversal
        // quadratic in the number of nodes.
        for (std::size_t head = 0; head < frontier.size(); ++head) {
            const MeshIndex current = frontier[head];
            std::vector<MeshIndex> discovered;
            for (std::size_t k = adjacency.offsets[current];
                 k < adjacency.offsets[current + 1]; ++k) {
                const MeshIndex neighbor = adjacency.neighbors[k];
                if (visited[neighbor])
                    continue;
                visited[neighbor] = 1;
                discovered.push_back(neighbor);
            }
            std::sort(discovered.begin(), discovered.end(),
                      [&adjacency](MeshIndex left, MeshIndex right) {
                          const std::size_t leftDegree = adjacency.degree(left);
                          const std::size_t rightDegree = adjacency.degree(right);
                          return leftDegree != rightDegree ? leftDegree < rightDegree
                                                           : left < right;
                      });
            for (MeshIndex node : discovered) {
                order.push_back(node);
                frontier.push_back(node);
            }
        }
    }

    std::vector<MeshIndex> oldToNew(m_nodeCount, InvalidMeshIndex);
    for (std::size_t position = 0; position < order.size(); ++position)
        oldToNew[order[position]] = static_cast<MeshIndex>(position);
    return oldToNew;
}

std::size_t LogicalMeshView::bandwidth(const Adjacency &adjacency,
                                       const std::vector<MeshIndex> &oldToNew)
{
    std::size_t widest = 0;
    for (std::size_t node = 0; node < oldToNew.size(); ++node) {
        const std::size_t position = oldToNew[node];
        for (std::size_t k = adjacency.offsets[node]; k < adjacency.offsets[node + 1]; ++k) {
            const std::size_t neighbor = oldToNew[adjacency.neighbors[k]];
            const std::size_t distance =
                position > neighbor ? position - neighbor : neighbor - position;
            widest = std::max(widest, distance);
        }
    }
    return widest + 1;
}

std::size_t LogicalMeshView::profile(const Adjacency &adjacency,
                                     const std::vector<MeshIndex> &oldToNew)
{
    std::size_t total = 0;
    for (std::size_t node = 0; node < oldToNew.size(); ++node) {
        const std::size_t position = oldToNew[node];
        std::size_t lowest = position;
        for (std::size_t k = adjacency.offsets[node]; k < adjacency.offsets[node + 1]; ++k)
            lowest = std::min(lowest, static_cast<std::size_t>(oldToNew[adjacency.neighbors[k]]));
        total += position - lowest;
    }
    return total;
}

std::size_t LogicalMeshView::storedByteCount() const
{
    std::size_t bytes = sizeof(LogicalMeshView);
    for (const auto &meshTemplate : m_templates) {
        const auto &local = meshTemplate.localMesh;
        bytes += local.nodes.size() * sizeof(SolverMesh::Node);
        bytes += local.elements.size() * sizeof(SolverMesh::Element);
        bytes += local.edges.size() * sizeof(SolverMesh::Edge);
        bytes += local.periodicConstraints.size() * sizeof(SolverMesh::PeriodicConstraint);
        for (const auto &gap : local.airGaps)
            bytes += sizeof(SolverMesh::AirGap) + gap.quadraturePoints.size() *
                                                      sizeof(SolverMesh::AirGapQuadraturePoint);
        for (const auto &seam : meshTemplate.seams) {
            bytes += seam.name.size();
            bytes += seam.orderedNodes.size() * sizeof(MeshIndex);
        }
    }
    bytes += m_instances.size() * sizeof(MeshInstance);
    for (const auto &instance : m_instances)
        bytes += instance.seamConnections.size() * sizeof(SeamConnection);
    bytes += m_instanceNodeBase.size() * sizeof(std::size_t);
    for (const auto &welds : m_weldedNodes)
        bytes += welds.size() * sizeof(std::pair<MeshIndex, MeshIndex>);
    bytes += m_weldedGlobalToLocal.size() *
             sizeof(std::pair<MeshIndex, std::pair<std::size_t, MeshIndex>>);
    bytes += m_instanceElementBase.size() * sizeof(std::size_t);
    bytes += m_markedEdges.size() *
             sizeof(std::pair<std::pair<MeshIndex, MeshIndex>, std::int32_t>);
    bytes += m_periodicConstraints.size() * sizeof(SolverMesh::PeriodicConstraint);
    for (const auto &gap : m_airGaps)
        bytes += sizeof(SolverMesh::AirGap) +
                 (gap.quadraturePoints.size() + gap.innerRing.size() + gap.outerRing.size() +
                  gap.nodeIndices.size()) *
                     sizeof(SolverMesh::AirGapQuadraturePoint);
    return bytes;
}

std::size_t LogicalMeshView::expandedByteCount() const
{
    std::size_t bytes = sizeof(SolverMesh);
    bytes += m_nodeCount * sizeof(SolverMesh::Node);
    bytes += m_elementCount * sizeof(SolverMesh::Element);
    bytes += edgeCount() * sizeof(SolverMesh::Edge);
    bytes += m_periodicConstraints.size() * sizeof(SolverMesh::PeriodicConstraint);
    for (const auto &gap : m_airGaps) {
        bytes += sizeof(SolverMesh::AirGap);
        bytes += gap.quadraturePoints.size() * sizeof(SolverMesh::AirGapQuadraturePoint);
        bytes += gap.innerRing.size() * sizeof(SolverMesh::AirGapRingPoint);
        bytes += gap.outerRing.size() * sizeof(SolverMesh::AirGapRingPoint);
        bytes += gap.nodeIndices.size() * sizeof(MeshIndex);
    }
    return bytes;
}

} // namespace mesh
} // namespace femm
