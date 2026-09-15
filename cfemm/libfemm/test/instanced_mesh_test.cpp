#include "mesh/InstancedMesh.h"
#include "mesh/SolverMeshValidator.h"

#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

namespace {

using namespace femm::mesh;

constexpr double Pi = 3.141592653589793238462643383279502884;

#define REQUIRE(condition)                                                            \
    do {                                                                              \
        if (!(condition)) {                                                           \
            std::cerr << __FILE__ << ':' << __LINE__ << ": " << #condition << '\n';   \
            return false;                                                             \
        }                                                                             \
    } while (false)

bool hasCategory(const std::vector<MaterializationDiagnostic> &diagnostics,
                 MaterializationDiagnosticCategory category)
{
    for (const auto &diagnostic : diagnostics)
        if (diagnostic.category == category)
            return true;
    return false;
}

SolverMesh makeTriangleMesh(double ax, double ay, double bx, double by, double cx, double cy,
                            std::int32_t boundaryMarker = 0)
{
    SolverMesh mesh;
    mesh.nodes = {{ax, ay, boundaryMarker}, {bx, by, 0}, {cx, cy, 0}};
    mesh.elements.push_back({{{0, 1, 2}}, 1});
    mesh.edges = {{0, 1, 0}, {1, 2, 0}, {2, 0, 0}};
    return mesh;
}

MeshTemplate makeTriangleTemplate(double ax, double ay, double bx, double by, double cx,
                                  double cy)
{
    MeshTemplate meshTemplate;
    meshTemplate.localMesh = makeTriangleMesh(ax, ay, bx, by, cx, cy);
    meshTemplate.seams = {{"edge", {0, 1}}};
    return meshTemplate;
}

MeshTemplate makeWedgeTemplate(double innerRadius, double outerRadius, double angleDegrees)
{
    const double radians = angleDegrees * Pi / 180.0;
    const double cosine = std::cos(radians);
    const double sine = std::sin(radians);
    SolverMesh mesh;
    mesh.nodes = {
        {innerRadius, 0.0, 0},
        {outerRadius, 0.0, 0},
        {innerRadius * cosine, innerRadius * sine, 0},
        {outerRadius * cosine, outerRadius * sine, 0},
    };
    mesh.elements.push_back({{{0, 1, 3}}, 1});
    mesh.elements.push_back({{{0, 3, 2}}, 1});
    mesh.edges = {{0, 1, 0}, {1, 3, 0}, {2, 3, 0}, {0, 2, 0}, {0, 3, 0}};
    MeshTemplate meshTemplate;
    meshTemplate.localMesh = mesh;
    meshTemplate.seams = {{"left", {0, 1}}, {"right", {2, 3}}};
    return meshTemplate;
}

InstancedMesh makeRing(double angleStep, std::size_t count)
{
    InstancedMesh instanced;
    instanced.templates.push_back(makeWedgeTemplate(1.0, 2.0, angleStep));
    for (std::size_t k = 0; k < count; ++k) {
        MeshInstance instance;
        instance.templateIndex = 0;
        instance.transform.rotationDegrees = static_cast<double>(k) * angleStep;
        SeamConnection connection;
        connection.seam = 1;
        connection.otherInstance = (k + 1) % count;
        connection.otherSeam = 0;
        instance.seamConnections.push_back(connection);
        instanced.instances.push_back(instance);
    }
    return instanced;
}

// C1: transforms and non-finite rejection.
bool testTransform()
{
    RigidTransform2D transform;
    transform.rotationDegrees = 90.0;
    transform.translationXMetres = 1.0;
    transform.translationYMetres = 2.0;
    double x = 0.0, y = 0.0;
    transform.applyPoint(1.0, 0.0, x, y);
    REQUIRE(std::abs(x - 1.0) < 1e-12);
    REQUIRE(std::abs(y - 3.0) < 1e-12);
    transform.applyDirection(1.0, 0.0, x, y);
    REQUIRE(std::abs(x) < 1e-12);
    REQUIRE(std::abs(y - 1.0) < 1e-12);
    REQUIRE(transform.isFinite());

    RigidTransform2D nonFinite;
    nonFinite.rotationDegrees = std::nan("");
    REQUIRE(!nonFinite.isFinite());

    InstancedMesh instanced;
    instanced.templates.push_back(makeTriangleTemplate(0, 0, 1, 0, 0, 1));
    MeshInstance instance;
    instance.transform = nonFinite;
    instanced.instances.push_back(instance);
    REQUIRE(hasCategory(instanced.validate(),
                        MaterializationDiagnosticCategory::NonFiniteTransform));
    return true;
}

// C2: one failure per template/instance/seam invariant.
bool testValidationCategories()
{
    InstancedMesh base;
    base.templates.push_back(makeTriangleTemplate(0, 0, 1, 0, 0, 1));
    MeshInstance first;
    first.seamConnections.push_back({0, 1, 0, SeamOrientation::Forward});
    MeshInstance second;
    second.transform.translationXMetres = 1.0;
    base.instances.push_back(first);
    base.instances.push_back(second);

    InstancedMesh invalidTemplate = base;
    invalidTemplate.instances[0].templateIndex = 99;
    REQUIRE(hasCategory(invalidTemplate.validate(),
                        MaterializationDiagnosticCategory::InvalidTemplateIndex));

    InstancedMesh emptySeam = base;
    emptySeam.templates[0].seams[0].orderedNodes.clear();
    REQUIRE(hasCategory(emptySeam.validate(),
                        MaterializationDiagnosticCategory::EmptySeam));

    InstancedMesh invalidSeamNode = base;
    invalidSeamNode.templates[0].seams[0].orderedNodes = {0, 99};
    REQUIRE(hasCategory(invalidSeamNode.validate(),
                        MaterializationDiagnosticCategory::InvalidSeamNode));

    InstancedMesh duplicateSeamNode = base;
    duplicateSeamNode.templates[0].seams[0].orderedNodes = {0, 0};
    REQUIRE(hasCategory(duplicateSeamNode.validate(),
                        MaterializationDiagnosticCategory::DuplicateSeamNode));

    InstancedMesh disconnected;
    MeshTemplate disconnectedTemplate;
    disconnectedTemplate.localMesh.nodes = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {2, 2, 0}};
    disconnectedTemplate.localMesh.elements.push_back({{{0, 1, 2}}, 1});
    disconnectedTemplate.localMesh.edges = {{0, 1, 0}, {1, 2, 0}, {2, 0, 0}};
    disconnectedTemplate.seams = {{"edge", {0, 3}}};
    disconnected.templates.push_back(disconnectedTemplate);
    disconnected.instances.push_back({});
    REQUIRE(hasCategory(disconnected.validate(),
                        MaterializationDiagnosticCategory::DisconnectedSeamChain));

    InstancedMesh invalidElementNode = base;
    invalidElementNode.templates[0].localMesh.elements[0].nodes[1] = 99;
    REQUIRE(hasCategory(invalidElementNode.validate(),
                        MaterializationDiagnosticCategory::InvalidElementNode));

    InstancedMesh invalidEdgeNode = base;
    invalidEdgeNode.templates[0].localMesh.edges[0].second = 99;
    REQUIRE(hasCategory(invalidEdgeNode.validate(),
                        MaterializationDiagnosticCategory::InvalidEdgeNode));

    InstancedMesh invalidSeam = base;
    invalidSeam.instances[0].seamConnections[0].seam = 99;
    REQUIRE(hasCategory(invalidSeam.validate(),
                        MaterializationDiagnosticCategory::InvalidInstanceSeam));

    InstancedMesh invalidInstance = base;
    invalidInstance.instances[0].seamConnections[0].otherInstance = 99;
    REQUIRE(hasCategory(invalidInstance.validate(),
                        MaterializationDiagnosticCategory::InvalidConnectionInstance));

    InstancedMesh invalidOtherSeam = base;
    invalidOtherSeam.instances[0].seamConnections[0].otherSeam = 99;
    REQUIRE(hasCategory(invalidOtherSeam.validate(),
                        MaterializationDiagnosticCategory::InvalidConnectionSeam));

    InstancedMesh selfConnection = base;
    selfConnection.instances[0].seamConnections[0].otherInstance = 0;
    REQUIRE(hasCategory(selfConnection.validate(),
                        MaterializationDiagnosticCategory::SelfConnection));

    // Two templates whose seams differ in length cannot be welded.
    InstancedMesh cardinality;
    MeshTemplate longSeam = makeTriangleTemplate(0, 0, 1, 0, 0, 1);
    MeshTemplate shortSeam;
    shortSeam.localMesh = makeTriangleMesh(0, 0, 1, 0, 0, 1);
    shortSeam.seams = {{"point", {0}}};
    cardinality.templates.push_back(longSeam);
    cardinality.templates.push_back(shortSeam);
    MeshInstance longInstance;
    longInstance.seamConnections.push_back({0, 1, 0, SeamOrientation::Forward});
    MeshInstance shortInstance;
    shortInstance.templateIndex = 1;
    cardinality.instances.push_back(longInstance);
    cardinality.instances.push_back(shortInstance);
    REQUIRE(hasCategory(cardinality.validate(),
                        MaterializationDiagnosticCategory::SeamCardinalityMismatch));

    InstancedMesh conflicting = base;
    conflicting.instances[1].seamConnections.push_back(
        {0, 0, 0, SeamOrientation::Reverse});
    REQUIRE(hasCategory(conflicting.validate(),
                        MaterializationDiagnosticCategory::ConflictingConnectionOrientation));
    return true;
}

// C3: translated pair keeps distinct unknowns; reversed seams weld.
bool testWelding()
{
    {
        InstancedMesh instanced;
        instanced.templates.push_back(makeTriangleTemplate(0, 0, 1, 0, 0, 1));
        MeshInstance first;
        MeshInstance second;
        second.transform.translationXMetres = 2.0;
        instanced.instances.push_back(first);
        instanced.instances.push_back(second);
        auto result = instanced.materialize();
        REQUIRE(result.succeeded());
        REQUIRE(result.mesh.nodes.size() == 6);
        REQUIRE(result.mesh.elements.size() == 2);
        REQUIRE(result.nodeMap[0][0] != result.nodeMap[1][0]);
        REQUIRE(result.nodeProvenance.size() == 6);
        REQUIRE(result.elementProvenance.size() == 2);
        REQUIRE(result.elementProvenance[1].instanceIndex == 1);
    }

    {
        // Reversed seam ordering must weld node i to node n-1-i.
        InstancedMesh instanced;
        instanced.templates.push_back(makeTriangleTemplate(1, 0, 1, 1, 0, 0));
        instanced.templates.push_back(makeTriangleTemplate(1, 1, 1, 0, 2, 0));
        MeshInstance first;
        first.seamConnections.push_back({0, 1, 0, SeamOrientation::Reverse});
        MeshInstance second;
        second.templateIndex = 1;
        instanced.instances.push_back(first);
        instanced.instances.push_back(second);
        auto result = instanced.materialize();
        REQUIRE(result.succeeded());
        REQUIRE(result.mesh.nodes.size() == 4);
        REQUIRE(result.mesh.elements.size() == 2);
        REQUIRE(result.nodeMap[0][0] == result.nodeMap[1][1]);
        REQUIRE(result.nodeMap[0][1] == result.nodeMap[1][0]);
    }

    {
        // Close but undeclared nodes must remain distinct.
        InstancedMesh instanced;
        instanced.templates.push_back(makeTriangleTemplate(0, 0, 1, 0, 0, 1));
        MeshInstance first;
        MeshInstance second;
        second.transform.translationXMetres = 1e-12;
        instanced.instances.push_back(first);
        instanced.instances.push_back(second);
        auto result = instanced.materialize();
        REQUIRE(result.succeeded());
        REQUIRE(result.mesh.nodes.size() == 6);
    }
    return true;
}

// C3/C4: a closed annular ring has the predicted counts and positive areas.
bool testClosedRing()
{
    const std::size_t count = 8;
    InstancedMesh instanced = makeRing(360.0 / static_cast<double>(count), count);
    auto result = instanced.materialize();
    REQUIRE(result.succeeded());
    REQUIRE(result.mesh.nodes.size() == 2 * count);
    REQUIRE(result.mesh.elements.size() == 2 * count);
    REQUIRE(result.mesh.edges.size() == 4 * count);
    for (const auto &element : result.mesh.elements) {
        const auto &a = result.mesh.nodes[element.nodes[0]];
        const auto &b = result.mesh.nodes[element.nodes[1]];
        const auto &c = result.mesh.nodes[element.nodes[2]];
        const double twiceArea = (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
        REQUIRE(twiceArea > 0.0);
    }
    REQUIRE(validateSolverMesh(result.mesh).valid());
    REQUIRE(result.nodeProvenance.size() == result.mesh.nodes.size());
    REQUIRE(result.elementProvenance.size() == result.mesh.elements.size());
    return true;
}

// C4: degenerate and reversed elements and marker conflicts are rejected.
bool testElementAndEdgeRejection()
{
    {
        InstancedMesh instanced;
        MeshTemplate meshTemplate;
        meshTemplate.localMesh.nodes = {{0, 0, 0}, {1, 0, 0}, {2, 0, 0}};
        meshTemplate.localMesh.elements.push_back({{{0, 1, 2}}, 1});
        meshTemplate.localMesh.edges = {{0, 1, 0}, {1, 2, 0}, {2, 0, 0}};
        instanced.templates.push_back(meshTemplate);
        instanced.instances.push_back({});
        REQUIRE(hasCategory(instanced.materialize().diagnostics,
                            MaterializationDiagnosticCategory::DegenerateElement));
    }
    {
        InstancedMesh instanced;
        MeshTemplate meshTemplate;
        meshTemplate.localMesh.nodes = {{0, 0, 0}, {0, 1, 0}, {1, 0, 0}};
        meshTemplate.localMesh.elements.push_back({{{0, 1, 2}}, 1});
        meshTemplate.localMesh.edges = {{0, 1, 0}, {1, 2, 0}, {2, 0, 0}};
        instanced.templates.push_back(meshTemplate);
        instanced.instances.push_back({});
        REQUIRE(hasCategory(instanced.materialize().diagnostics,
                            MaterializationDiagnosticCategory::ReversedElement));
    }
    {
        // A shared seam edge with disagreeing boundary markers is rejected.
        InstancedMesh instanced;
        MeshTemplate first = makeTriangleTemplate(1, 0, 1, 1, 0, 0);
        first.localMesh.edges[0].boundaryMarker = 1;
        MeshTemplate second = makeTriangleTemplate(1, 1, 1, 0, 2, 0);
        second.localMesh.edges[0].boundaryMarker = 2;
        instanced.templates.push_back(first);
        instanced.templates.push_back(second);
        MeshInstance firstInstance;
        firstInstance.seamConnections.push_back({0, 1, 0, SeamOrientation::Reverse});
        MeshInstance secondInstance;
        secondInstance.templateIndex = 1;
        instanced.instances.push_back(firstInstance);
        instanced.instances.push_back(secondInstance);
        REQUIRE(hasCategory(instanced.materialize().diagnostics,
                            MaterializationDiagnosticCategory::EdgeMarkerConflict));
    }
    return true;
}

// C5: periodic and AGE topology are remapped through the provenance table.
bool testPeriodicAndAgeRemap()
{
    MeshTemplate meshTemplate;
    SolverMesh &mesh = meshTemplate.localMesh;
    mesh.nodes = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}};
    mesh.elements.push_back({{{0, 1, 2}}, 1});
    mesh.elements.push_back({{{0, 2, 3}}, 1});
    mesh.edges = {{0, 1, 0}, {1, 2, 0}, {2, 3, 0}, {3, 0, 0}, {0, 2, 0}};
    mesh.periodicConstraints.push_back({0, 2, SolverMesh::Periodicity::Antiperiodic});

    SolverMesh::AirGap gap;
    gap.boundaryName = "gap";
    gap.periodicity = SolverMesh::Periodicity::Periodic;
    gap.totalArcElements = 1;
    gap.totalArcLengthDegrees = 360.0;
    gap.innerRadius = 1.0;
    gap.outerRadius = 2.0;
    gap.centerX = 0.5;
    gap.centerY = 0.5;
    gap.innerAngleDegrees = 10.0;
    gap.outerAngleDegrees = 20.0;
    gap.innerRing.push_back({0, 0.0, 1.0});
    gap.outerRing.push_back({1, 0.0, 1.0});
    gap.quadraturePoints.resize(2);
    gap.quadraturePoints[0].nodes = {{0, 0, 1, 1}};
    gap.quadraturePoints[0].weights = {{1, 1, 1, 1}};
    gap.quadraturePoints[1].nodes = {{0, 0, 1, 1}};
    gap.quadraturePoints[1].weights = {{1, 1, 1, 1}};
    gap.nodeIndices = {0, 1};
    mesh.airGaps.push_back(gap);

    InstancedMesh instanced;
    instanced.templates.push_back(meshTemplate);
    MeshInstance first;
    MeshInstance second;
    second.transform.translationXMetres = 10.0;
    instanced.instances.push_back(first);
    instanced.instances.push_back(second);

    auto result = instanced.materialize();
    REQUIRE(result.succeeded());
    REQUIRE(result.mesh.periodicConstraints.size() == 2);
    REQUIRE(result.mesh.airGaps.size() == 2);
    for (std::size_t i = 0; i < 2; ++i) {
        REQUIRE(result.mesh.periodicConstraints[i].first == result.nodeMap[i][0]);
        REQUIRE(result.mesh.periodicConstraints[i].second == result.nodeMap[i][2]);
        REQUIRE(result.mesh.periodicConstraints[i].periodicity ==
                SolverMesh::Periodicity::Antiperiodic);
        const auto &mapped = result.mesh.airGaps[i];
        REQUIRE(mapped.innerRing[0].node == result.nodeMap[i][0]);
        REQUIRE(mapped.outerRing[0].node == result.nodeMap[i][1]);
        REQUIRE(mapped.nodeIndices[0] == result.nodeMap[i][0]);
        REQUIRE(mapped.nodeIndices[1] == result.nodeMap[i][1]);
        REQUIRE(mapped.quadraturePoints[0].nodes[0] == result.nodeMap[i][0]);
        REQUIRE(mapped.quadraturePoints[0].nodes[2] == result.nodeMap[i][1]);
        REQUIRE(std::abs(mapped.centerX - (0.5 + 10.0 * static_cast<double>(i))) < 1e-12);
        REQUIRE(mapped.periodicity == SolverMesh::Periodicity::Periodic);
    }

    // Invalid local references are rejected.
    InstancedMesh badPeriodic;
    MeshTemplate badPeriodicTemplate = meshTemplate;
    badPeriodicTemplate.localMesh.periodicConstraints[0].first = 99;
    badPeriodic.instances.push_back({});
    badPeriodic.templates.push_back(badPeriodicTemplate);
    REQUIRE(hasCategory(badPeriodic.materialize().diagnostics,
                        MaterializationDiagnosticCategory::InvalidPeriodicNode));

    InstancedMesh badAge;
    MeshTemplate badAgeTemplate = meshTemplate;
    badAgeTemplate.localMesh.airGaps[0].innerRing[0].node = 99;
    badAge.instances.push_back({});
    badAge.templates.push_back(badAgeTemplate);
    REQUIRE(hasCategory(badAge.materialize().diagnostics,
                        MaterializationDiagnosticCategory::InvalidAirGapNode));
    return true;
}

// A rotated instance must rotate its AGE centre and starting angles.
bool testAgeRotation()
{
    MeshTemplate meshTemplate;
    SolverMesh &mesh = meshTemplate.localMesh;
    mesh.nodes = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}};
    mesh.elements.push_back({{{0, 1, 2}}, 1});
    mesh.edges = {{0, 1, 0}, {1, 2, 0}, {2, 0, 0}};
    SolverMesh::AirGap gap;
    gap.boundaryName = "gap";
    gap.totalArcElements = 1;
    gap.totalArcLengthDegrees = 360.0;
    gap.innerRadius = 1.0;
    gap.outerRadius = 2.0;
    gap.centerX = 1.0;
    gap.centerY = 0.0;
    gap.innerAngleDegrees = 10.0;
    gap.outerAngleDegrees = 20.0;
    gap.innerRing.push_back({0, 0.0, 1.0});
    gap.outerRing.push_back({1, 0.0, 1.0});
    gap.quadraturePoints.resize(2);
    gap.quadraturePoints[0].nodes = {{0, 0, 1, 1}};
    gap.quadraturePoints[0].weights = {{1, 1, 1, 1}};
    gap.quadraturePoints[1].nodes = {{0, 0, 1, 1}};
    gap.quadraturePoints[1].weights = {{1, 1, 1, 1}};
    mesh.airGaps.push_back(gap);

    InstancedMesh instanced;
    instanced.templates.push_back(meshTemplate);
    MeshInstance instance;
    instance.transform.rotationDegrees = 90.0;
    instanced.instances.push_back(instance);
    auto result = instanced.materialize();
    REQUIRE(result.succeeded());
    REQUIRE(result.mesh.airGaps.size() == 1);
    REQUIRE(std::abs(result.mesh.airGaps[0].centerX) < 1e-12);
    REQUIRE(std::abs(result.mesh.airGaps[0].centerY - 1.0) < 1e-12);
    REQUIRE(std::abs(result.mesh.airGaps[0].innerAngleDegrees - 100.0) < 1e-12);
    REQUIRE(std::abs(result.mesh.airGaps[0].outerAngleDegrees - 110.0) < 1e-12);
    return true;
}

} // namespace

int main()
{
    const std::vector<std::pair<std::string, bool (*)()>> tests = {
        {"transform", testTransform},
        {"validation-categories", testValidationCategories},
        {"welding", testWelding},
        {"closed-ring", testClosedRing},
        {"element-edge-rejection", testElementAndEdgeRejection},
        {"periodic-age-remap", testPeriodicAndAgeRemap},
        {"age-rotation", testAgeRotation},
    };
    bool success = true;
    for (const auto &test : tests) {
        if (!test.second()) {
            std::cerr << "FAILED: " << test.first << '\n';
            success = false;
        }
    }
    if (success)
        std::cout << "instanced mesh materializer tests passed\n";
    return success ? 0 : 1;
}
