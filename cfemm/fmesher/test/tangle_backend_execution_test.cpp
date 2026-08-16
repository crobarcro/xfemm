#include "TangleMesherBackend.h"

#include "FemmProblem.h"

#include <tangle_mesh.h>

#include <iostream>
#include <string>

namespace {

int fail(const std::string &message)
{
    std::cerr << message << '\n';
    return 1;
}

Mesh engineMesh()
{
    Mesh mesh;
    mesh.vertices.resize(3);
    mesh.vertices[0].x = 0.0; mesh.vertices[0].y = 0.0; mesh.vertices[0].id = 0;
    mesh.vertices[1].x = 1.0; mesh.vertices[1].y = 0.0; mesh.vertices[1].id = 1;
    mesh.vertices[2].x = 0.0; mesh.vertices[2].y = 1.0; mesh.vertices[2].id = 2;
    Triangle triangle;
    triangle.v = {{0, 1, 2}};
    triangle.region_attrib = 0.0;
    mesh.triangles.push_back(triangle);
    mesh.segments = {{0, 1, 3}};
    mesh.edges = {{0, 1}, {1, 2}, {2, 0}};
    mesh.pbc_pairs.push_back({0, 1, 0});
    return mesh;
}

} // namespace

int main()
{
    femm::FemmProblem problem(femm::FileType::MagneticsFile);
    problem.pathName = "in-memory-engine-probe.fem";
    int calls = 0;
    std::string receivedPath;
    fmesher::TangleMesherBackend backend(
        [&](const std::string &path, Mesh &mesh) {
            ++calls;
            receivedPath = path;
            mesh = engineMesh();
            return TANGLE_OK;
        });

    auto result = backend.mesh(problem, true);
    if (!result.succeeded() || calls != 1 || receivedPath != problem.pathName)
        return fail("TangleMesherBackend did not execute the injected engine exactly once");
    if (result.mesh.nodes.size() != 3 || result.mesh.elements.size() != 1 ||
        result.mesh.periodicConstraints.size() != 1)
        return fail("TangleMesherBackend did not convert the engine result");

    result = backend.mesh(problem, false);
    if (!result.succeeded() || calls != 2 || !result.mesh.periodicConstraints.empty())
        return fail("non-periodic Tangle request retained periodic constraints");

    fmesher::TangleMesherBackend failingBackend(
        [](const std::string &, Mesh &) { return TANGLE_ERR_MESH; });
    result = failingBackend.mesh(problem, true);
    if (result.succeeded() || result.status != femm::mesh::MeshStatus::BackendFailure ||
        result.diagnostics.size() != 1 ||
        result.diagnostics[0].backendName != "Tangle" ||
        result.diagnostics[0].backendErrorCode != TANGLE_ERR_MESH)
        return fail("Tangle engine failure was not mapped to an actionable diagnostic");

    problem.pathName.clear();
    result = backend.mesh(problem, true);
    if (result.succeeded() || result.status != femm::mesh::MeshStatus::InvalidInput ||
        calls != 2)
        return fail("pathless problem was not rejected before engine execution");

    return 0;
}
