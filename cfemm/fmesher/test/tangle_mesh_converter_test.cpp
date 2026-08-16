#include "TangleMeshConverter.h"

#include <tangle_mesh.h>

#include <cmath>
#include <iostream>
#include <string>

namespace {

int fail(const std::string &message)
{
    std::cerr << message << '\n';
    return 1;
}

Mesh ordinaryMesh()
{
    Mesh mesh;
    mesh.vertices.resize(3);
    mesh.vertices[0].x = 0.0; mesh.vertices[0].y = 0.0; mesh.vertices[0].id = 0; mesh.vertices[0].marker = 2;
    mesh.vertices[1].x = 2.0; mesh.vertices[1].y = 0.0; mesh.vertices[1].id = 1; mesh.vertices[1].marker = 3;
    mesh.vertices[2].x = 0.0; mesh.vertices[2].y = 2.0; mesh.vertices[2].id = 2; mesh.vertices[2].marker = 0;
    Triangle triangle;
    triangle.v = {{0, 1, 2}};
    triangle.region_attrib = 4.0;
    mesh.triangles.push_back(triangle);
    mesh.segments = {{0, 1, 9}};
    mesh.edges = {{0, 1}, {1, 2}, {2, 0}};
    return mesh;
}

} // namespace

int main()
{
    auto source = ordinaryMesh();
    auto result = fmesher::convertTangleMesh(source, 0.001);
    if (!result.succeeded())
        return fail("ordinary Tangle mesh conversion failed");
    if (result.mesh.nodes.size() != 3 || result.mesh.elements.size() != 1 ||
        result.mesh.edges.size() != 3)
        return fail("ordinary Tangle mesh conversion lost topology");
    if (std::abs(result.mesh.nodes[1].x - 0.002) > 1e-15 ||
        result.mesh.nodes[0].boundaryMarker != 2 ||
        result.mesh.elements[0].regionAttribute != 4 ||
        result.mesh.edges[0].boundaryMarker != 9 ||
        result.mesh.edges[1].boundaryMarker != 0)
        return fail("ordinary Tangle mesh conversion lost values");

    source.pbc_pairs.push_back({0, 1, 1});
    result = fmesher::convertTangleMesh(source, 1.0);
    if (!result.succeeded() || result.mesh.periodicConstraints.size() != 1 ||
        result.mesh.periodicConstraints[0].periodicity !=
            femm::mesh::SolverMesh::Periodicity::Antiperiodic)
        return fail("Tangle periodic pair conversion failed");

    Mesh annulus;
    annulus.vertices.resize(4);
    annulus.vertices[0].x = 1.0; annulus.vertices[0].y = 0.0; annulus.vertices[0].id = 0;
    annulus.vertices[1].x = 0.0; annulus.vertices[1].y = 1.0; annulus.vertices[1].id = 1;
    annulus.vertices[2].x = 2.0; annulus.vertices[2].y = 0.0; annulus.vertices[2].id = 2;
    annulus.vertices[3].x = 0.0; annulus.vertices[3].y = 2.0; annulus.vertices[3].id = 3;
    AGEDef age;
    age.name = "gap";
    age.format = 0;
    age.innerAngle = 0.0;
    age.outerAngle = 0.0;
    age.ri = 1.0;
    age.ro = 2.0;
    age.totalArcLength = 180.0;
    age.cx = 0.0;
    age.cy = 0.0;
    age.n = 2;
    age.innerNodes = {0, 1};
    age.outerNodes = {2, 3};
    annulus.age_defs.push_back(age);
    result = fmesher::convertTangleMesh(annulus, 0.01);
    if (!result.succeeded() || result.mesh.airGaps.size() != 1)
        return fail("Tangle AGE conversion failed");
    const auto &gap = result.mesh.airGaps[0];
    if (gap.boundaryName != "gap" || gap.totalArcElements != 2 ||
        gap.innerRing.size() != 4 || gap.outerRing.size() != 4 ||
        gap.quadraturePoints.size() != 3 || gap.nodeIndices.size() != 8 ||
        std::abs(gap.innerRadius - 0.01) > 1e-15 ||
        std::abs(gap.outerRadius - 0.02) > 1e-15)
        return fail("Tangle AGE conversion produced incomplete topology");

    annulus.age_defs[0].format = 1;
    result = fmesher::convertTangleMesh(annulus, 1.0);
    if (!result.succeeded() || result.mesh.airGaps[0].innerRing[2].weight != -1.0 ||
        result.mesh.airGaps[0].outerRing[2].weight != -1.0)
        return fail("Tangle antiperiodic AGE signs were not preserved");

    annulus.age_defs[0].innerNodes[1] = 99;
    if (fmesher::convertTangleMesh(annulus, 1.0).succeeded())
        return fail("invalid Tangle AGE node reference was not diagnosed");

    source = ordinaryMesh();
    source.triangles[0].v[2] = 99;
    result = fmesher::convertTangleMesh(source, 1.0);
    if (result.succeeded() || result.status != femm::mesh::MeshStatus::ConversionFailure ||
        result.diagnostics.empty() || result.diagnostics[0].backendName != "Tangle")
        return fail("invalid Tangle node reference was not diagnosed");

    source = ordinaryMesh();
    source.pbc_pairs.push_back({0, 1, 7});
    if (fmesher::convertTangleMesh(source, 1.0).succeeded())
        return fail("invalid Tangle periodicity was not diagnosed");

    if (fmesher::convertTangleMesh(ordinaryMesh(), 0.0).succeeded())
        return fail("invalid Tangle length scale was not diagnosed");

    return 0;
}
