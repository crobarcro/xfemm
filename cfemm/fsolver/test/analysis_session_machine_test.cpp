#include "AnalysisSession.h"
#include "FSolverAnalysisBackend.h"

#include "CBlockLabel.h"
#include "CBoundaryProp.h"
#include "CCircuit.h"
#include "CMaterialProp.h"
#include "CPointProp.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using namespace femm;
using namespace femm::mesh;

constexpr double Pi = 3.141592653589793238462643383279502884;

std::unique_ptr<FemmProblem> makeModel()
{
    auto problem = std::make_unique<FemmProblem>(FileType::MagneticsFile);

    auto point = std::make_unique<CMPointProp>();
    point->PointName = "fixed";
    point->A = 0;
    problem->nodeproplist.push_back(std::move(point));

    auto coil = std::make_unique<CMMaterialProp>();
    coil->BlockName = "coil";
    coil->mu_x = coil->mu_y = 1;
    coil->J = CComplex(1, 0);
    problem->blockproplist.push_back(std::move(coil));

    auto magnet = std::make_unique<CMMaterialProp>();
    magnet->BlockName = "magnet";
    magnet->mu_x = magnet->mu_y = 1;
    magnet->H_c = 1;
    problem->blockproplist.push_back(std::move(magnet));

    auto circuit = std::make_unique<CMCircuit>();
    circuit->CircName = "phase";
    circuit->CircType = 1;
    circuit->Amps = 1;
    problem->circproplist.push_back(std::move(circuit));

    auto coilLabel = std::make_unique<CMBlockLabel>();
    coilLabel->BlockType = 0;
    coilLabel->BlockTypeName = "coil";
    coilLabel->InCircuit = 0;
    coilLabel->Turns = 1;
    problem->labellist.push_back(std::move(coilLabel));

    auto magnetLabel = std::make_unique<CMBlockLabel>();
    magnetLabel->BlockType = 1;
    magnetLabel->BlockTypeName = "magnet";
    problem->labellist.push_back(std::move(magnetLabel));

    auto age = std::make_unique<CMBoundaryProp>();
    age->BdryName = "sliding-gap";
    age->BdryFormat = 6;
    problem->lineproplist.push_back(std::move(age));
    return problem;
}

/** One annular sector tile with a radial seam on each edge. */
MeshTemplate annularTemplate(double innerRadius, double outerRadius, double sectorDegrees,
                             int angularDivisions, bool airGapOnInner,
                             std::int32_t regionAttribute)
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
        mesh.elements.push_back({{{innerNodes[j], outerNodes[j], outerNodes[j + 1]}}, regionAttribute});
        mesh.elements.push_back(
            {{{innerNodes[j], outerNodes[j + 1], innerNodes[j + 1]}}, regionAttribute});
    }
    std::set<std::pair<MeshIndex, MeshIndex>> edges;
    for (const auto &element : mesh.elements)
        for (std::size_t j = 0; j < 3; ++j)
            edges.insert(std::minmax(element.nodes[j], element.nodes[(j + 1) % 3]));
    for (const auto &edge : edges)
        mesh.edges.push_back({edge.first, edge.second, 0});

    // Pin A = 0 at one node (marker 2 selects the first point property) so the
    // pure-Neumann problem is well posed.
    mesh.nodes.front().boundaryMarker = 2;

    meshTemplate.seams.push_back({"airgap", airGapOnInner ? innerNodes : outerNodes});
    meshTemplate.seams.push_back({"left", {innerNodes.front(), outerNodes.front()}});
    meshTemplate.seams.push_back({"right", {innerNodes.back(), outerNodes.back()}});
    return meshTemplate;
}

InstancedMesh makeMachine()
{
    InstancedMesh instanced;
    // Rotor pole template (180 degrees, magnet) and stator slot template
    // (90 degrees, coil). Independent counts: 2 poles and 4 slots.
    instanced.templates.push_back(annularTemplate(0.04, 0.05, 180.0, 4, false, 2));
    instanced.templates.push_back(annularTemplate(0.06, 0.08, 90.0, 2, true, 1));

    const auto addRing = [&](std::size_t templateIndex, std::size_t count, double angleStep) {
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
    addRing(0, 2, 180.0);
    addRing(1, 4, 90.0);

    AirGapCoupling coupling;
    coupling.innerTemplate = 0;
    coupling.innerSeam = 0;
    coupling.outerTemplate = 1;
    coupling.outerSeam = 0;
    coupling.boundaryName = "sliding-gap";
    coupling.innerRadiusMetres = 0.05;
    coupling.outerRadiusMetres = 0.06;
    coupling.totalArcLengthDegrees = 360.0;
    instanced.airGapCouplings.push_back(coupling);
    return instanced;
}

bool testIndependentDomainsAndSweep()
{
    auto solver = std::make_shared<FSolverAnalysisBackend>();
    AnalysisSession session(ModelDefinition(makeModel()), solver);
    session.setInstancedMesh(makeMachine());

    const auto mesh = session.ensureMesh();
    assert(mesh->airGaps.size() == 1);
    const auto &gap = mesh->airGaps.front();
    assert(gap.innerRing.size() == 8);
    assert(gap.outerRing.size() == 8);
    assert(gap.quadraturePoints.size() == 9);
    assert(session.materializationCount() == 1);

    session.solve();
    assert(solver->topologyImportCount() == 1);
    const std::size_t importedTopology = session.meshTopologyIdentity();
    const std::size_t materializations = session.materializationCount();

    // Rotor-position sweep: changing the AGE angle reuses the topology and
    // does not re-materialise or re-import the solver mesh.
    for (double angle : {15.0, 30.0, 45.0}) {
        session.setAirGapAngle(session.model().airGap("sliding-gap"), angle, 0.0);
        session.solve();
    }
    assert(solver->solveCount() == 4);
    assert(session.materializationCount() == materializations);
    assert(solver->topologyImportCount() == 1);
    assert(session.meshTopologyIdentity() == importedTopology);

    // A transform change re-materialises and re-imports. Rotate the whole ring
    // so the declared seams stay coincident.
    const std::size_t instanceCount = session.instancedMesh()->instances.size();
    for (std::size_t i = 0; i < instanceCount; ++i) {
        RigidTransform2D moved = session.instancedMesh()->instances[i].transform;
        moved.rotationDegrees += 5.0;
        session.setInstanceTransform(i, moved);
    }
    session.solve();
    assert(session.materializationCount() == materializations + 1);
    assert(solver->topologyImportCount() == 2);
    return true;
}

} // namespace

int main()
{
    try {
        if (!testIndependentDomainsAndSweep())
            return 1;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    std::cout << "instanced machine session smoke test passed\n";
    return 0;
}
