#include "mesh/InstancedMesh.h"
#include "mesh/LogicalMeshView.h"
#include "mesh/SolverMeshValidator.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <set>
#include <string>
#include <vector>

namespace {

using namespace femm::mesh;

constexpr double Pi = 3.141592653589793238462643383279502884;

#define REQUIRE(condition)                                                          \
    do {                                                                            \
        if (!(condition)) {                                                         \
            std::cerr << __FILE__ << ':' << __LINE__ << ": " << #condition << '\n'; \
            return false;                                                           \
        }                                                                           \
    } while (false)

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

MeshTemplate makeBlockTemplate(int nx, int ny)
{
    SolverMesh mesh;
    for (int j = 0; j <= ny; ++j)
        for (int i = 0; i <= nx; ++i)
            mesh.nodes.push_back({static_cast<double>(i), static_cast<double>(j), 0});
    const auto index = [nx](int i, int j) {
        return static_cast<MeshIndex>(j * (nx + 1) + i);
    };
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            mesh.elements.push_back({{{index(i, j), index(i + 1, j), index(i + 1, j + 1)}}, 1});
            mesh.elements.push_back({{{index(i, j), index(i + 1, j + 1), index(i, j + 1)}}, 1});
        }
    }
    for (int i = 0; i < nx; ++i) {
        mesh.edges.push_back({index(i, 0), index(i + 1, 0), 0});
        mesh.edges.push_back({index(i + 1, ny), index(i, ny), 0});
    }
    for (int j = 0; j < ny; ++j) {
        mesh.edges.push_back({index(nx, j), index(nx, j + 1), 0});
        mesh.edges.push_back({index(0, j + 1), index(0, j), 0});
    }
    MeshTemplate meshTemplate;
    meshTemplate.localMesh = mesh;
    return meshTemplate;
}

MeshTemplate makeSquareTemplate()
{
    MeshTemplate meshTemplate;
    SolverMesh &mesh = meshTemplate.localMesh;
    mesh.nodes = {{0, 0, 2}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}};
    mesh.elements.push_back({{{0, 1, 2}}, 1});
    mesh.elements.push_back({{{0, 2, 3}}, 1});
    mesh.edges = {{0, 1, 0}, {1, 2, 0}, {2, 3, 0}, {3, 0, 0}, {0, 2, 0}};
    return meshTemplate;
}

MeshTemplate makeTriangleTemplate(double ax, double ay, double bx, double by, double cx,
                                  double cy)
{
    MeshTemplate meshTemplate;
    SolverMesh &mesh = meshTemplate.localMesh;
    mesh.nodes = {{ax, ay, 0}, {bx, by, 0}, {cx, cy, 0}};
    mesh.elements.push_back({{{0, 1, 2}}, 1});
    mesh.edges = {{0, 1, 0}, {1, 2, 0}, {2, 0, 0}};
    meshTemplate.seams = {{"edge", {0, 1}}};
    return meshTemplate;
}

std::vector<std::set<MeshIndex>> meshAdjacency(const SolverMesh &mesh)
{
    std::vector<std::set<MeshIndex>> adjacency(mesh.nodes.size());
    for (const auto &element : mesh.elements) {
        for (std::size_t k = 0; k < 3; ++k) {
            const MeshIndex first = element.nodes[k];
            const MeshIndex second = element.nodes[(k + 1) % 3];
            adjacency[first].insert(second);
            adjacency[second].insert(first);
        }
    }
    return adjacency;
}

bool compareAgainstMaterialization(const InstancedMesh &instanced)
{
    std::vector<MaterializationDiagnostic> diagnostics;
    const LogicalMeshView view = LogicalMeshView::build(instanced, diagnostics);
    REQUIRE(diagnostics.empty());
    REQUIRE(view.valid());

    const auto result = instanced.materialize();
    REQUIRE(result.succeeded());
    const SolverMesh &mesh = result.mesh;

    REQUIRE(view.nodeCount() == mesh.nodes.size());
    REQUIRE(view.elementCount() == mesh.elements.size());
    REQUIRE(view.instanceCount() == instanced.instances.size());
    REQUIRE(view.edgeCount() == mesh.edges.size());
    REQUIRE(view.periodicConstraints().size() == mesh.periodicConstraints.size());
    REQUIRE(view.airGaps().size() == mesh.airGaps.size());

    for (std::size_t i = 0; i < instanced.instances.size(); ++i) {
        for (MeshIndex n = 0;
             n < instanced.templates[instanced.instances[i].templateIndex].localMesh.nodes.size();
             ++n)
            REQUIRE(view.nodeFor(i, n) == result.provenance.nodeMap[i][n]);
        for (MeshIndex e = 0; e < instanced.templates[instanced.instances[i].templateIndex]
                                      .localMesh.elements.size();
             ++e)
            REQUIRE(view.elementFor(i, e) == result.provenance.elementMap[i][e]);
    }

    for (MeshIndex n = 0; n < mesh.nodes.size(); ++n) {
        const NodeProvenance provenance = view.nodeProvenance(n);
        REQUIRE(provenance.templateIndex == result.provenance.nodeProvenance[n].templateIndex);
        REQUIRE(provenance.instanceIndex == result.provenance.nodeProvenance[n].instanceIndex);
        REQUIRE(provenance.localNode == result.provenance.nodeProvenance[n].localNode);
        double x = 0.0, y = 0.0;
        view.nodeCoordinates(n, x, y);
        REQUIRE(std::abs(x - mesh.nodes[n].x) < 1e-12);
        REQUIRE(std::abs(y - mesh.nodes[n].y) < 1e-12);
    }

    for (MeshIndex e = 0; e < mesh.elements.size(); ++e) {
        const ElementProvenance provenance = view.elementProvenance(e);
        REQUIRE(provenance.templateIndex == result.provenance.elementProvenance[e].templateIndex);
        REQUIRE(provenance.instanceIndex == result.provenance.elementProvenance[e].instanceIndex);
        REQUIRE(provenance.localElement == result.provenance.elementProvenance[e].localElement);
        const auto nodes = view.elementNodes(e);
        REQUIRE(nodes == mesh.elements[e].nodes);
        REQUIRE(view.elementRegionAttribute(e) == mesh.elements[e].regionAttribute);
        REQUIRE(view.elementInstance(e) == result.provenance.elementProvenance[e].instanceIndex);
    }

    for (const auto &edge : mesh.edges) {
        REQUIRE(view.boundaryMarkerForEdge(edge.first, edge.second) == edge.boundaryMarker);
        REQUIRE(view.boundaryMarkerForEdge(edge.second, edge.first) == edge.boundaryMarker);
    }

    for (std::size_t c = 0; c < mesh.periodicConstraints.size(); ++c) {
        REQUIRE(view.periodicConstraints()[c].first == mesh.periodicConstraints[c].first);
        REQUIRE(view.periodicConstraints()[c].second == mesh.periodicConstraints[c].second);
        REQUIRE(view.periodicConstraints()[c].periodicity ==
                mesh.periodicConstraints[c].periodicity);
    }
    for (std::size_t g = 0; g < mesh.airGaps.size(); ++g) {
        REQUIRE(view.airGaps()[g].boundaryName == mesh.airGaps[g].boundaryName);
        REQUIRE(view.airGaps()[g].totalArcElements == mesh.airGaps[g].totalArcElements);
        REQUIRE(view.airGaps()[g].nodeIndices == mesh.airGaps[g].nodeIndices);
        REQUIRE(view.airGaps()[g].innerRing.size() == mesh.airGaps[g].innerRing.size());
    }
    return true;
}

// G1: global DOFs and welded seams agree with the materialised oracle.
bool testViewMatchesMaterialization()
{
    REQUIRE(compareAgainstMaterialization(makeRing(45.0, 8)));

    // Independent translated instances keep distinct unknowns.
    InstancedMesh translated;
    translated.templates.push_back(makeSquareTemplate());
    MeshInstance first;
    MeshInstance second;
    second.transform.translationXMetres = 2.0;
    translated.instances.push_back(first);
    translated.instances.push_back(second);
    REQUIRE(compareAgainstMaterialization(translated));

    // An open chain leaves a node welded to a previous instance followed by an
    // un-welded node; the inverse provenance must skip the welded node.
    InstancedMesh chain;
    chain.templates.push_back(makeTriangleTemplate(0, 0, 1, 0, 0, 1));
    MeshInstance chainFirst;
    MeshInstance chainSecond;
    chainFirst.seamConnections.push_back({0, 1, 0, SeamOrientation::Forward});
    chain.instances.push_back(chainFirst);
    chain.instances.push_back(chainSecond);
    REQUIRE(compareAgainstMaterialization(chain));
    return true;
}

// G1/C5: periodic constraints and AGE topology survive the view.
bool testPeriodicAndAgeView()
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
    REQUIRE(compareAgainstMaterialization(instanced));
    return true;
}

// G2: adjacency matches the materialised graph; Cuthill reduces profile.
bool testAdjacencyAndOrdering()
{
    const std::size_t count = 12;
    InstancedMesh instanced = makeRing(360.0 / static_cast<double>(count), count);
    std::vector<MaterializationDiagnostic> diagnostics;
    const LogicalMeshView view = LogicalMeshView::build(instanced, diagnostics);
    REQUIRE(view.valid());

    const auto result = instanced.materialize();
    REQUIRE(result.succeeded());
    const auto expected = meshAdjacency(result.mesh);
    const auto adjacency = view.buildAdjacency();
    REQUIRE(adjacency.offsets.size() == view.nodeCount() + 1);
    for (std::size_t node = 0; node < view.nodeCount(); ++node) {
        std::set<MeshIndex> neighbors;
        for (std::size_t k = adjacency.offsets[node]; k < adjacency.offsets[node + 1]; ++k)
            neighbors.insert(adjacency.neighbors[k]);
        REQUIRE(neighbors == expected[node]);
    }

    // Natural ordering then Cuthill-McKee: bandwidth and profile must not grow.
    std::vector<MeshIndex> natural(view.nodeCount());
    for (std::size_t i = 0; i < natural.size(); ++i)
        natural[i] = static_cast<MeshIndex>(i);
    const std::size_t naturalBandwidth = LogicalMeshView::bandwidth(adjacency, natural);
    const std::size_t naturalProfile = LogicalMeshView::profile(adjacency, natural);
    const std::vector<MeshIndex> ordering = view.cuthillMcKeeOrdering();
    const std::size_t orderedBandwidth = LogicalMeshView::bandwidth(adjacency, ordering);
    const std::size_t orderedProfile = LogicalMeshView::profile(adjacency, ordering);
    REQUIRE(orderedBandwidth <= naturalBandwidth);
    REQUIRE(orderedProfile <= naturalProfile);
    return true;
}

// G6: the compact representation does not scale with logical element count.
bool testMemoryScaling()
{
    const MeshTemplate block = makeBlockTemplate(6, 6);
    const auto makeInstances = [&block](std::size_t count) {
        InstancedMesh instanced;
        instanced.templates.push_back(block);
        for (std::size_t k = 0; k < count; ++k) {
            MeshInstance instance;
            instance.transform.translationXMetres = 2.0 * static_cast<double>(k);
            instanced.instances.push_back(instance);
        }
        return instanced;
    };

    std::vector<MaterializationDiagnostic> diagnostics;
    const LogicalMeshView smallView = LogicalMeshView::build(makeInstances(8), diagnostics);
    REQUIRE(smallView.valid());
    const LogicalMeshView largeView = LogicalMeshView::build(makeInstances(128), diagnostics);
    REQUIRE(largeView.valid());

    REQUIRE(largeView.elementCount() == 128 * block.localMesh.elements.size());
    REQUIRE(largeView.storedByteCount() < largeView.expandedByteCount());

    // Adding instances expands the materialised mesh much faster than it grows
    // the compact view: the view adds instance metadata, the mesh adds full
    // element connectivity.
    const std::size_t storedDelta = largeView.storedByteCount() - smallView.storedByteCount();
    const std::size_t expandedDelta =
        largeView.expandedByteCount() - smallView.expandedByteCount();
    REQUIRE(storedDelta < expandedDelta);
    return true;
}

} // namespace

int main()
{
    const std::vector<std::pair<std::string, bool (*)()>> tests = {
        {"view-matches-materialization", testViewMatchesMaterialization},
        {"periodic-and-age-view", testPeriodicAndAgeView},
        {"adjacency-and-ordering", testAdjacencyAndOrdering},
        {"memory-scaling", testMemoryScaling},
    };
    bool success = true;
    for (const auto &test : tests) {
        if (!test.second()) {
            std::cerr << "FAILED: " << test.first << '\n';
            success = false;
        }
    }
    if (success)
        std::cout << "logical mesh view tests passed\n";
    return success ? 0 : 1;
}
