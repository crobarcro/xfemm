#include "mesh/SolverMeshValidator.h"

#include <cmath>
#include <iostream>
#include <limits>

using femm::mesh::SolverMesh;
using femm::mesh::SolverMeshValidationCategory;

namespace {

SolverMesh ordinaryMesh()
{
    SolverMesh mesh;
    mesh.nodes = {{0.0, 0.0, 0}, {1.0, 0.0, 0}, {0.0, 1.0, 0}};
    mesh.elements.push_back({{{0, 1, 2}}, 1});
    mesh.edges = {{0, 1, 0}, {1, 2, 0}, {2, 0, 0}};
    return mesh;
}

SolverMesh ageMesh()
{
    SolverMesh mesh = ordinaryMesh();
    SolverMesh::AirGap gap;
    gap.boundaryName = "synthetic gap";
    gap.totalArcElements = 1;
    gap.totalArcLengthDegrees = 90.0;
    gap.innerRadius = 0.5;
    gap.outerRadius = 1.0;
    gap.quadraturePoints.push_back({{{0, 1, 1, 2}}, {{0.5, 0.5, 0.5, 0.5}}});
    gap.quadraturePoints.push_back({{{1, 2, 2, 0}}, {{0.5, 0.5, 0.5, 0.5}}});
    gap.innerRing = {{0, 0.0, 1.0}, {1, 1.0, 1.0}};
    gap.outerRing = {{1, 0.0, 1.0}, {2, 1.0, 1.0}};
    gap.nodeIndices = {0, 1, 1, 2};
    mesh.airGaps.push_back(gap);
    return mesh;
}

bool has(const SolverMesh &mesh, SolverMeshValidationCategory category)
{
    const auto result = femm::mesh::validateSolverMesh(mesh);
    for (const auto &diagnostic : result.diagnostics)
        if (diagnostic.category == category)
            return true;
    return false;
}

bool hasDiagnostic(const SolverMesh &mesh, SolverMeshValidationCategory category,
                   std::size_t objectIndex, std::size_t localIndex, const std::string &context)
{
    const auto result = femm::mesh::validateSolverMesh(mesh);
    for (const auto &diagnostic : result.diagnostics)
        if (diagnostic.category == category && diagnostic.objectIndex == objectIndex &&
            diagnostic.localIndex == localIndex && diagnostic.context == context)
            return true;
    return false;
}

int fail(const char *message)
{
    std::cerr << message << '\n';
    return 1;
}

#define EXPECT_CATEGORY(mesh, category) \
    if (!has((mesh), SolverMeshValidationCategory::category)) return fail("missing " #category)

} // namespace

int main()
{
    if (!femm::mesh::validateSolverMesh(ordinaryMesh()).valid())
        return fail("ordinary mesh was rejected");

    auto mesh = ordinaryMesh();
    mesh.periodicConstraints.push_back({0, 1, SolverMesh::Periodicity::Antiperiodic});
    if (!femm::mesh::validateSolverMesh(mesh).valid() || !mesh.airGaps.empty())
        return fail("ordinary periodic mesh without AGE was rejected");
    mesh.periodicConstraints = {
        {0, 0, SolverMesh::Periodicity::Periodic},
        {1, 1, SolverMesh::Periodicity::Antiperiodic}
    };
    if (!femm::mesh::validateSolverMesh(mesh).valid())
        return fail("valid periodic and antiperiodic self-pairs were rejected");
    if (!femm::mesh::validateSolverMesh(ageMesh()).valid())
        return fail("synthetic AGE mesh was rejected");
    mesh = ageMesh();
    mesh.airGaps[0].innerRing.clear();
    mesh.airGaps[0].outerRing.clear();
    if (!femm::mesh::validateSolverMesh(mesh).valid())
        return fail("valid legacy-compatible AGE without reusable rings was rejected");

    mesh = ordinaryMesh(); mesh.nodes[0].x = std::numeric_limits<double>::quiet_NaN();
    EXPECT_CATEGORY(mesh, NonFiniteNodeCoordinate);
    if (!hasDiagnostic(mesh, SolverMeshValidationCategory::NonFiniteNodeCoordinate, 0, 0, "x"))
        return fail("node diagnostic did not identify its coordinate");
    mesh = ordinaryMesh(); mesh.nodes[0].y = std::numeric_limits<double>::infinity();
    EXPECT_CATEGORY(mesh, NonFiniteNodeCoordinate);
    mesh = ordinaryMesh(); mesh.elements[0].nodes[2] = 99;
    EXPECT_CATEGORY(mesh, InvalidElementNode);
    if (!hasDiagnostic(mesh, SolverMeshValidationCategory::InvalidElementNode, 0, 2, "element node"))
        return fail("element diagnostic did not identify its local node");
    mesh = ordinaryMesh(); mesh.elements[0].nodes[2] = 1;
    EXPECT_CATEGORY(mesh, RepeatedElementNode);
    mesh = ordinaryMesh(); mesh.nodes[2] = {2.0, 0.0, 0};
    EXPECT_CATEGORY(mesh, NonPositiveElementArea);
    mesh = ordinaryMesh(); std::swap(mesh.elements[0].nodes[1], mesh.elements[0].nodes[2]);
    EXPECT_CATEGORY(mesh, NonPositiveElementArea);
    mesh = ordinaryMesh(); mesh.nodes[1].x = std::numeric_limits<double>::max();
    mesh.nodes[2].y = std::numeric_limits<double>::max();
    EXPECT_CATEGORY(mesh, NonFiniteElementArea);

    mesh = ordinaryMesh(); mesh.edges[0].first = 99;
    EXPECT_CATEGORY(mesh, InvalidEdgeNode);
    mesh = ordinaryMesh(); mesh.edges[0].second = mesh.edges[0].first;
    EXPECT_CATEGORY(mesh, SelfEdge);
    mesh = ordinaryMesh(); mesh.nodes.push_back({2.0, 2.0, 0}); mesh.edges[0] = {0, 3, 0};
    EXPECT_CATEGORY(mesh, EdgeNotOwnedByElement);
    mesh = ordinaryMesh(); mesh.periodicConstraints.push_back({0, 99, SolverMesh::Periodicity::Periodic});
    EXPECT_CATEGORY(mesh, InvalidPeriodicNode);
    mesh = ordinaryMesh(); mesh.periodicConstraints.push_back({0, 1, static_cast<SolverMesh::Periodicity>(99)});
    EXPECT_CATEGORY(mesh, InvalidPeriodicity);

    mesh = ageMesh(); mesh.airGaps[0].centerX = std::numeric_limits<double>::quiet_NaN();
    EXPECT_CATEGORY(mesh, NonFiniteAirGapValue);
    mesh = ageMesh(); mesh.airGaps[0].outerRadius = mesh.airGaps[0].innerRadius;
    EXPECT_CATEGORY(mesh, InvalidAirGapGeometry);
    mesh = ageMesh(); mesh.airGaps[0].outerRing.pop_back();
    EXPECT_CATEGORY(mesh, InvalidAirGapStructure);
    mesh = ageMesh(); mesh.airGaps[0].quadraturePoints.pop_back();
    EXPECT_CATEGORY(mesh, InvalidAirGapStructure);
    mesh = ageMesh(); mesh.airGaps[0].quadraturePoints.push_back(mesh.airGaps[0].quadraturePoints.back());
    EXPECT_CATEGORY(mesh, InvalidAirGapStructure);
    mesh = ageMesh(); mesh.airGaps[0].totalArcElements = 3;
    mesh.airGaps[0].quadraturePoints.resize(4, mesh.airGaps[0].quadraturePoints.back());
    EXPECT_CATEGORY(mesh, InvalidAirGapStructure);
    mesh = ageMesh(); mesh.airGaps[0].quadraturePoints[0].nodes[2] = 99;
    EXPECT_CATEGORY(mesh, InvalidAirGapQuadratureNode);
    if (!hasDiagnostic(mesh, SolverMeshValidationCategory::InvalidAirGapQuadratureNode,
                       0, 2, "synthetic gap"))
        return fail("AGE diagnostic did not identify its gap and quadrature node");
    mesh = ageMesh(); mesh.airGaps[0].quadraturePoints[0].weights[1] = std::numeric_limits<double>::infinity();
    EXPECT_CATEGORY(mesh, NonFiniteAirGapQuadratureWeight);
    mesh = ageMesh(); mesh.airGaps[0].innerRing[0].node = 99;
    EXPECT_CATEGORY(mesh, InvalidAirGapRingNode);
    if (!hasDiagnostic(mesh, SolverMeshValidationCategory::InvalidAirGapRingNode,
                       0, 0, "synthetic gap/innerRing"))
        return fail("AGE ring diagnostic did not identify its ring");
    mesh = ageMesh(); mesh.airGaps[0].outerRing[0].elementPosition = std::numeric_limits<double>::quiet_NaN();
    EXPECT_CATEGORY(mesh, NonFiniteAirGapRingValue);
    mesh = ageMesh(); mesh.airGaps[0].innerRing[0].weight = std::numeric_limits<double>::infinity();
    EXPECT_CATEGORY(mesh, NonFiniteAirGapRingValue);
    mesh = ageMesh(); mesh.airGaps[0].nodeIndices[2] = 99;
    EXPECT_CATEGORY(mesh, InvalidAirGapNodeIndex);

    return 0;
}
