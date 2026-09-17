#include "AnalysisSession.h"
#include "FSolverAnalysisBackend.h"

#include "CBlockLabel.h"
#include "CCircuit.h"
#include "CMaterialProp.h"
#include "CPointProp.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

namespace {

using namespace femm;
using namespace femm::mesh;

std::unique_ptr<FemmProblem> makeProblem()
{
    auto problem = std::make_unique<FemmProblem>(FileType::MagneticsFile);

    auto point = std::make_unique<CMPointProp>();
    point->PointName = "fixed";
    point->A = 0;
    problem->nodeproplist.push_back(std::move(point));

    auto material = std::make_unique<CMMaterialProp>();
    material->BlockName = "copper";
    material->mu_x = material->mu_y = 1;
    material->H_c = 1;
    problem->blockproplist.push_back(std::move(material));

    auto circuit = std::make_unique<CMCircuit>();
    circuit->CircName = "phase";
    circuit->CircType = 1;
    circuit->Amps = 1;
    problem->circproplist.push_back(std::move(circuit));

    auto label = std::make_unique<CMBlockLabel>();
    label->BlockType = 0;
    label->BlockTypeName = "copper";
    label->InCircuit = 0;
    label->Turns = 1;
    problem->labellist.push_back(std::move(label));
    return problem;
}

MeshTemplate makeSquareTemplate()
{
    MeshTemplate meshTemplate;
    SolverMesh &mesh = meshTemplate.localMesh;
    // Node 0 carries point-property marker 2 (SolverMesh marker m maps to point
    // property m - 2), pinning A = 0 so the patch is well posed.
    mesh.nodes = {{0, 0, 2}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}};
    mesh.elements.push_back({{{0, 1, 2}}, 1});
    mesh.elements.push_back({{{0, 2, 3}}, 1});
    mesh.edges = {{0, 1, 0}, {1, 2, 0}, {2, 3, 0}, {3, 0, 0}, {0, 2, 0}};
    return meshTemplate;
}

InstancedMesh makeInstanced()
{
    InstancedMesh instanced;
    instanced.templates.push_back(makeSquareTemplate());
    MeshInstance first;
    first.templateIndex = 0;
    MeshInstance second;
    second.templateIndex = 0;
    second.transform.translationXMetres = 2.0;
    instanced.instances.push_back(first);
    instanced.instances.push_back(second);
    return instanced;
}

bool testIndependentPhysicsAndCaches()
{
    auto solver = std::make_shared<FSolverAnalysisBackend>();
    AnalysisSession session(ModelDefinition(makeProblem()), solver);
    session.setInstancedMesh(makeInstanced());

    session.synchronize();
    const auto mesh = session.mesh();
    // Two instances of a one-label template produce two independent labels.
    for (const auto &element : mesh->elements)
        assert(element.regionAttribute >= 1 && element.regionAttribute <= 2);
    assert(session.prepared().labels.size() == 2);
    assert(session.prepared().circuits.size() == 2);
    assert(session.materializationCount() == 1);
    assert(session.templateTopologyIdentity() != 0);
    assert(session.instanceLayoutIdentity() != 0);
    assert(session.instancingProvenance().elementProvenance.size() == mesh->elements.size());
    assert(session.elementsForInstance(0).size() + session.elementsForInstance(1).size() ==
           mesh->elements.size());
    assert(session.elementsForInstance(0).size() == 2);
    assert(session.nodesForInstance(0).size() == 4);

    // Instancing alone introduces no periodic constraint.
    assert(mesh->periodicConstraints.empty());

    // Alternating physics on the second instance: reversed current and a
    // rotated magnetisation, resolved without touching the template.
    std::vector<InstanceRegionOverride> overrides;
    InstanceRegionOverride override0;
    override0.sourceBlockLabel = 0;
    override0.currentScale = -1.0;
    override0.magnetisationRotationDegrees = 180.0;
    overrides.push_back(override0);
    session.setInstanceRegionOverrides(1, overrides);

    session.synchronize();
    assert(session.prepared().labels[0].Turns == 1);
    assert(session.prepared().labels[1].Turns == -1);
    assert(session.prepared().labels[0].MagDir == 0.0);
    assert(session.prepared().labels[1].MagDir == 180.0);
    assert(session.prepared().labels[0].InCircuit == 0);
    assert(session.prepared().labels[1].InCircuit == 0);
    // Physics changes reuse topology: no re-materialisation and no re-import.
    assert(session.materializationCount() == 1);
    assert(solver->topologyImportCount() == 1);

    const auto trial = session.solve();
    assert(solver->solveCount() == 1);
    assert(trial.real);
    bool nonTrivial = false;
    for (double value : trial.real->nodal.magneticVectorPotential) {
        assert(std::isfinite(value));
        nonTrivial = nonTrivial || std::abs(value) > 1e-12;
    }
    assert(nonTrivial);

    // A transform change re-materialises and re-imports topology. A rigid
    // rotation also carries the constant magnetisation with it: the instance
    // rotation (90) is added before the explicit override delta (180).
    RigidTransform2D moved;
    moved.translationXMetres = 4.0;
    moved.rotationDegrees = 90.0;
    session.setInstanceTransform(1, moved);
    session.synchronize();
    assert(session.materializationCount() == 2);
    assert(solver->topologyImportCount() == 2);
    assert(session.prepared().labels[1].MagDir == 270.0);

    // Circuit, material, and solve-parameter changes never touch topology.
    session.setCircuitCurrent(session.model().circuit("phase"), CComplex(2, 0));
    session.setMaterialProperty(session.model().material("copper"),
                                MaterialProperty::Conductivity, 5.8);
    session.synchronize();
    assert(session.materializationCount() == 2);
    assert(solver->topologyImportCount() == 2);

    // Clearing instancing returns to the mesher path.
    session.clearInstancedMesh();
    assert(!session.instancedMesh());
    return true;
}

} // namespace

int main()
{
    try {
        testIndependentPhysicsAndCaches();
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
