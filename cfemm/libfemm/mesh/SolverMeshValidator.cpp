#include "SolverMeshValidator.h"

#include <algorithm>
#include <cmath>
#include <set>

namespace femm {
namespace mesh {
namespace {

void add(SolverMeshValidationResult &result, SolverMeshValidationCategory category,
         std::size_t objectIndex, std::size_t localIndex, const std::string &context)
{
    result.diagnostics.push_back({category, objectIndex, localIndex, context});
}

bool validNode(MeshIndex node, std::size_t count)
{
    return node < count;
}

bool validPeriodicity(SolverMesh::Periodicity value)
{
    return value == SolverMesh::Periodicity::Periodic ||
           value == SolverMesh::Periodicity::Antiperiodic;
}

} // namespace

SolverMeshValidationResult validateSolverMesh(const SolverMesh &mesh)
{
    SolverMeshValidationResult result;
    std::set<std::pair<MeshIndex, MeshIndex>> elementEdges;
    for (std::size_t i = 0; i < mesh.nodes.size(); ++i) {
        if (!std::isfinite(mesh.nodes[i].x))
            add(result, SolverMeshValidationCategory::NonFiniteNodeCoordinate, i, 0, "x");
        if (!std::isfinite(mesh.nodes[i].y))
            add(result, SolverMeshValidationCategory::NonFiniteNodeCoordinate, i, 1, "y");
    }

    for (std::size_t i = 0; i < mesh.elements.size(); ++i) {
        const auto &element = mesh.elements[i];
        bool connectivityValid = true;
        for (std::size_t j = 0; j < element.nodes.size(); ++j) {
            if (!validNode(element.nodes[j], mesh.nodes.size())) {
                add(result, SolverMeshValidationCategory::InvalidElementNode, i, j, "element node");
                connectivityValid = false;
            }
        }
        if (!connectivityValid)
            continue;
        if (element.nodes[0] == element.nodes[1] || element.nodes[1] == element.nodes[2] ||
            element.nodes[2] == element.nodes[0]) {
            add(result, SolverMeshValidationCategory::RepeatedElementNode, i, 0, "element");
            continue;
        }
        const auto &a = mesh.nodes[element.nodes[0]];
        const auto &b = mesh.nodes[element.nodes[1]];
        const auto &c = mesh.nodes[element.nodes[2]];
        const double twiceArea = (b.x - a.x) * (c.y - a.y) -
                                 (b.y - a.y) * (c.x - a.x);
        if (!std::isfinite(twiceArea))
            add(result, SolverMeshValidationCategory::NonFiniteElementArea, i, 0, "signed area");
        else if (twiceArea <= 0.0)
            add(result, SolverMeshValidationCategory::NonPositiveElementArea, i, 0, "signed area");
        for (std::size_t j = 0; j < 3; ++j) {
            const MeshIndex first = element.nodes[j];
            const MeshIndex second = element.nodes[(j + 1) % 3];
            elementEdges.insert(std::minmax(first, second));
        }
    }

    for (std::size_t i = 0; i < mesh.edges.size(); ++i) {
        const auto &edge = mesh.edges[i];
        bool connectivityValid = true;
        if (!validNode(edge.first, mesh.nodes.size())) {
            add(result, SolverMeshValidationCategory::InvalidEdgeNode, i, 0, "first");
            connectivityValid = false;
        }
        if (!validNode(edge.second, mesh.nodes.size())) {
            add(result, SolverMeshValidationCategory::InvalidEdgeNode, i, 1, "second");
            connectivityValid = false;
        }
        if (!connectivityValid)
            continue;
        if (edge.first == edge.second) {
            add(result, SolverMeshValidationCategory::SelfEdge, i, 0, "edge");
            continue;
        }
        if (elementEdges.count(std::minmax(edge.first, edge.second)) == 0)
            add(result, SolverMeshValidationCategory::EdgeNotOwnedByElement, i, 0, "edge");
    }

    for (std::size_t i = 0; i < mesh.periodicConstraints.size(); ++i) {
        const auto &constraint = mesh.periodicConstraints[i];
        if (!validNode(constraint.first, mesh.nodes.size()))
            add(result, SolverMeshValidationCategory::InvalidPeriodicNode, i, 0, "first");
        if (!validNode(constraint.second, mesh.nodes.size()))
            add(result, SolverMeshValidationCategory::InvalidPeriodicNode, i, 1, "second");
        if (validNode(constraint.first, mesh.nodes.size()) && constraint.first == constraint.second)
            add(result, SolverMeshValidationCategory::SelfPeriodicConstraint, i, 0, "constraint");
        if (!validPeriodicity(constraint.periodicity))
            add(result, SolverMeshValidationCategory::InvalidPeriodicity, i, 0, "constraint");
    }

    for (std::size_t i = 0; i < mesh.airGaps.size(); ++i) {
        const auto &gap = mesh.airGaps[i];
        const std::string context = gap.boundaryName.empty() ? "<unnamed>" : gap.boundaryName;
        const double scalars[] = {gap.totalArcLengthDegrees, gap.innerRadius, gap.outerRadius,
                                  gap.innerAngleDegrees, gap.outerAngleDegrees, gap.innerShift,
                                  gap.outerShift, gap.centerX, gap.centerY};
        for (std::size_t j = 0; j < sizeof(scalars) / sizeof(scalars[0]); ++j)
            if (!std::isfinite(scalars[j]))
                add(result, SolverMeshValidationCategory::NonFiniteAirGapValue, i, j, context);
        if (!validPeriodicity(gap.periodicity))
            add(result, SolverMeshValidationCategory::InvalidPeriodicity, i, 0, context);
        if (gap.totalArcElements == 0 || gap.totalArcLengthDegrees <= 0.0 ||
            gap.innerRadius <= 0.0 || gap.outerRadius <= gap.innerRadius)
            add(result, SolverMeshValidationCategory::InvalidAirGapGeometry, i, 0, context);
        if (gap.innerRing.empty() != gap.outerRing.empty() ||
            (!gap.innerRing.empty() && gap.innerRing.size() != gap.outerRing.size()))
            add(result, SolverMeshValidationCategory::InvalidAirGapStructure, i, 0, context);

        for (std::size_t q = 0; q < gap.quadraturePoints.size(); ++q) {
            for (std::size_t j = 0; j < 4; ++j) {
                if (!validNode(gap.quadraturePoints[q].nodes[j], mesh.nodes.size()))
                    add(result, SolverMeshValidationCategory::InvalidAirGapQuadratureNode, i, q * 4 + j, context);
                if (!std::isfinite(gap.quadraturePoints[q].weights[j]))
                    add(result, SolverMeshValidationCategory::NonFiniteAirGapQuadratureWeight, i, q * 4 + j, context);
            }
        }
        const std::vector<SolverMesh::AirGapRingPoint> *rings[] = {&gap.innerRing, &gap.outerRing};
        for (std::size_t r = 0; r < 2; ++r) {
            const std::string ringContext = context + (r == 0 ? "/innerRing" : "/outerRing");
            for (std::size_t j = 0; j < rings[r]->size(); ++j) {
                const auto &point = (*rings[r])[j];
                if (!validNode(point.node, mesh.nodes.size()))
                    add(result, SolverMeshValidationCategory::InvalidAirGapRingNode, i, j, ringContext);
                if (!std::isfinite(point.elementPosition) || !std::isfinite(point.weight))
                    add(result, SolverMeshValidationCategory::NonFiniteAirGapRingValue, i, j, ringContext);
            }
        }
        for (std::size_t j = 0; j < gap.nodeIndices.size(); ++j)
            if (!validNode(gap.nodeIndices[j], mesh.nodes.size()))
                add(result, SolverMeshValidationCategory::InvalidAirGapNodeIndex, i, j, context);
    }
    return result;
}

} // namespace mesh
} // namespace femm
