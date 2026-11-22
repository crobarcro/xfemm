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
        for (std::size_t j = 0; j < 3; ++j) {
            const MeshIndex first = element.nodes[j];
            const MeshIndex second = element.nodes[(j + 1) % 3];
            edges.insert(std::minmax(first, second));
        }
    for (const auto &edge : edges)
        mesh.edges.push_back({edge.first, edge.second, 0});

    meshTemplate.seams.push_back({"airgap", seamOnInner ? innerNodes : outerNodes});
    meshTemplate.seams.push_back({"left", {innerNodes.front(), outerNodes.front()}});
    meshTemplate.seams.push_back({"right", {innerNodes.back(), outerNodes.back()}});
    return meshTemplate;
}

InstancedMesh makeTwoDomainRing()
{
    InstancedMesh instanced;
    // Template 0: rotor pole (90 degrees). Template 1: stator slot (45 degrees).
    // Both air-gap arcs are discretised into 16 elements per full ring.
    instanced.templates.push_back(annularTemplate(0.05, 0.06, 90.0, 4, false));
    instanced.templates.push_back(annularTemplate(0.07, 0.09, 45.0, 2, true));

    const auto addRing = [&](std::size_t templateIndex, std::size_t count,
                             double angleStep) {
        const std::size_t base = instanced.instances.size();
        for (std::size_t k = 0; k < count; ++k) {
            MeshInstance instance;
            instance.templateIndex = templateIndex;
            instance.transform =
                RigidTransform2D::rotationAbout(0.0, 0.0, static_cast<double>(k) * angleStep);
            instance.seamConnections.push_back(
                {2, base + (k + 1) % count, 1, SeamOrientation::Forward});
            instanced.instances.push_back(std::move(instance));
        }
    };
    addRing(0, 4, 90.0);
    addRing(1, 8, 45.0);

    AirGapCoupling coupling;
    coupling.innerTemplate = 0;
    coupling.innerSeam = 0;
    coupling.outerTemplate = 1;
    coupling.outerSeam = 0;
    coupling.boundaryName = "sliding-gap";
    coupling.centerXMetres = 0.0;
    coupling.centerYMetres = 0.0;
    coupling.innerRadiusMetres = 0.06;
    coupling.outerRadiusMetres = 0.07;
    coupling.totalArcLengthDegrees = 360.0;
    instanced.airGapCouplings.push_back(coupling);
    return instanced;
}

bool testIndependentCountsAndAirGap()
{
    InstancedMesh instanced = makeTwoDomainRing();
    auto result = instanced.materialize();
    if (!result.succeeded())
        for (const auto &diagnostic : result.diagnostics)
            std::cerr << "diag cat=" << static_cast<int>(diagnostic.category) << " inst="
                      << diagnostic.instanceIndex << " obj=" << diagnostic.objectIndex
                      << " msg=" << diagnostic.message << '\n';
    REQUIRE(result.succeeded());

    // Independent instance counts: 4 rotor poles and 8 stator slots.
    REQUIRE(instanced.instances.size() == 12);
    REQUIRE(result.mesh.airGaps.size() == 1);
    const auto &gap = result.mesh.airGaps.front();
    REQUIRE(gap.boundaryName == "sliding-gap");
    REQUIRE(gap.innerRing.size() == 16);
    REQUIRE(gap.outerRing.size() == 16);
    REQUIRE(gap.totalArcElements == 16);
    REQUIRE(gap.quadraturePoints.size() == 17);
    REQUIRE(gap.nodeIndices.size() == 32);
    for (const auto &point : gap.innerRing)
        REQUIRE(point.node < result.mesh.nodes.size());
    for (const auto &point : gap.outerRing)
        REQUIRE(point.node < result.mesh.nodes.size());

    // The two domains are not welded: every rotor air-gap node is distinct
    // from every stator air-gap node.
    for (const auto &inner : gap.innerRing)
        for (const auto &outer : gap.outerRing)
            REQUIRE(inner.node != outer.node);

    REQUIRE(validateSolverMesh(result.mesh).valid());

    // An invalid coupling is rejected before materialisation.
    InstancedMesh invalid = instanced;
    invalid.airGapCouplings.front().innerSeam = 99;
    REQUIRE(!invalid.materialize().succeeded());
    return true;
}

} // namespace

int main()
{
    if (!testIndependentCountsAndAirGap()) {
        std::cerr << "FAILED: independent-counts-airgap\n";
        return 1;
    }
    std::cout << "instanced multi-domain air-gap tests passed\n";
    return 0;
}
