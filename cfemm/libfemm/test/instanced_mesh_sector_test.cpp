#include "mesh/InstancedMesh.h"
#include "mesh/SolverMeshValidator.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <set>
#include <string>
#include <utility>
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

/** One annular sector tile; the inner or outer arc is the air-gap seam. */
MeshTemplate annularTemplate(double innerRadius, double outerRadius, double sectorDegrees,
                             int angularDivisions, bool seamOnInner)
{
    MeshTemplate meshTemplate;
    SolverMesh &mesh = meshTemplate.localMesh;
    std::vector<MeshIndex> innerNodes;
    std::vector<MeshIndex> outerNodes;
    for (int j = 0; j <= angularDivisions; ++j) {
        const double angle = sectorDegrees * Pi / 180.0 * j / angularDivisions;
        innerNodes.push_back(mesh.nodes.size());
        mesh.nodes.push_back({innerRadius * std::cos(angle), innerRadius * std::sin(angle), 0});
        outerNodes.push_back(mesh.nodes.size());
        mesh.nodes.push_back({outerRadius * std::cos(angle), outerRadius * std::sin(angle), 0});
    }
    for (int j = 0; j < angularDivisions; ++j) {
        mesh.elements.push_back({{{innerNodes[j], outerNodes[j], outerNodes[j + 1]}}, 1});
        mesh.elements.push_back(
            {{{innerNodes[j], outerNodes[j + 1], innerNodes[j + 1]}}, 1});
    }
    std::set<std::pair<MeshIndex, MeshIndex>> edges;
    for (const auto &element : mesh.elements)
        for (std::size_t j = 0; j < 3; ++j)
            edges.insert(std::minmax(element.nodes[j], element.nodes[(j + 1) % 3]));
    for (const auto &edge : edges)
        mesh.edges.push_back({edge.first, edge.second, 0});

    meshTemplate.seams.push_back({"airgap", seamOnInner ? innerNodes : outerNodes});
    meshTemplate.seams.push_back({"left", {innerNodes.front(), outerNodes.front()}});
    meshTemplate.seams.push_back({"right", {innerNodes.back(), outerNodes.back()}});
    return meshTemplate;
}

/** A sector ring: internal joins are welded, the ends are periodic. */
void addSectorRing(InstancedMesh &instanced, std::size_t templateIndex, std::size_t count,
                   double sectorDegrees)
{
    const std::size_t base = instanced.instances.size();
    const double step = sectorDegrees / static_cast<double>(count);
    for (std::size_t k = 0; k < count; ++k) {
        MeshInstance instance;
        instance.templateIndex = templateIndex;
        instance.transform =
            RigidTransform2D::rotationAbout(0.0, 0.0, static_cast<double>(k) * step);
        if (k + 1 < count)
            instance.seamConnections.push_back(
                {2, base + k + 1, 1, SeamOrientation::Forward});
        instanced.instances.push_back(std::move(instance));
    }
    PeriodicClosure closure;
    closure.firstInstance = base + count - 1;
    closure.firstSeam = 2;
    closure.secondInstance = base;
    closure.secondSeam = 1;
    closure.periodicity = SolverMesh::Periodicity::Periodic;
    instanced.periodicClosures.push_back(closure);
}

InstancedMesh makeSector()
{
    InstancedMesh instanced;
    // Rotor pole tile (30 degrees) and stator slot tile (10 degrees) both
    // discretise the air-gap arc into 7 nodes over the 60-degree sector.
    instanced.templates.push_back(annularTemplate(0.04, 0.05, 30.0, 3, false));
    instanced.templates.push_back(annularTemplate(0.06, 0.08, 10.0, 1, true));
    addSectorRing(instanced, 0, 2, 60.0);
    addSectorRing(instanced, 1, 6, 60.0);

    AirGapCoupling coupling;
    coupling.innerTemplate = 0;
    coupling.innerSeam = 0;
    coupling.outerTemplate = 1;
    coupling.outerSeam = 0;
    coupling.boundaryName = "sliding-gap";
    coupling.totalArcLengthDegrees = 60.0;
    coupling.innerRadiusMetres = 0.05;
    coupling.outerRadiusMetres = 0.06;
    instanced.airGapCouplings.push_back(coupling);
    return instanced;
}

bool testSectorClosure()
{
    InstancedMesh instanced = makeSector();
    auto result = instanced.materialize();
    if (!result.succeeded())
        for (const auto &diagnostic : result.diagnostics)
            std::cerr << "diag cat=" << static_cast<int>(diagnostic.category)
                      << " msg=" << diagnostic.message << '\n';
    REQUIRE(result.succeeded());
    REQUIRE(validateSolverMesh(result.mesh).valid());

    // Internal joins are welded; the two end seams stay distinct and are
    // linked by four periodic constraints (two nodes per template).
    REQUIRE(result.mesh.periodicConstraints.size() == 4);
    for (const auto &constraint : result.mesh.periodicConstraints) {
        REQUIRE(constraint.first != constraint.second);
        REQUIRE(constraint.periodicity == SolverMesh::Periodicity::Periodic);
        REQUIRE(constraint.first < result.mesh.nodes.size());
        REQUIRE(constraint.second < result.mesh.nodes.size());
    }

    REQUIRE(result.mesh.airGaps.size() == 1);
    const auto &gap = result.mesh.airGaps.front();
    REQUIRE(gap.innerRing.size() == 7);
    REQUIRE(gap.outerRing.size() == 7);
    REQUIRE(gap.totalArcElements == 7);

    // The end ring nodes are not welded: every inner ring node is distinct
    // from every other, and likewise for the outer ring.
    std::set<MeshIndex> inner;
    for (const auto &point : gap.innerRing)
        REQUIRE(inner.insert(point.node).second);
    std::set<MeshIndex> outer;
    for (const auto &point : gap.outerRing)
        REQUIRE(outer.insert(point.node).second);
    return true;
}

bool testInvalidClosureRejected()
{
    InstancedMesh instanced = makeSector();
    instanced.periodicClosures.front().secondSeam = 99;
    const auto result = instanced.materialize();
    REQUIRE(!result.succeeded());
    bool found = false;
    for (const auto &diagnostic : result.diagnostics)
        found = found || diagnostic.category ==
                             MaterializationDiagnosticCategory::InvalidPeriodicClosure;
    REQUIRE(found);
    return true;
}

} // namespace

int main()
{
    if (!testSectorClosure() || !testInvalidClosureRejected()) {
        std::cerr << "FAILED: instanced sector closure\n";
        return 1;
    }
    std::cout << "instanced sector closure tests passed\n";
    return 0;
}
