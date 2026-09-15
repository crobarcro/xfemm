#include "TangleMesherBackend.h"

#include "FemmProblem.h"

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
    ::MeshOptions receivedOptions;
    fmesher::TangleMesherBackend backend(
        [&](const std::string &path, const ::MeshOptions &options, Mesh &mesh) {
            ++calls;
            receivedPath = path;
            receivedOptions = options;
            mesh = engineMesh();
            return TANGLE_OK;
        });

    femm::mesh::MeshingOptions options;
    options.minimumAngleDegrees = 28.0;
    options.defaultElementSize = 0.25;
    options.forceMaximumElementArea = true;
    options.suppressExteriorSteinerPoints = true;
    options.suppressUnusedVertices = true;
    options.verbose = true;
    auto result = backend.mesh(problem, true, options);
    if (!result.succeeded() || calls != 1 || receivedPath != problem.pathName)
        return fail("TangleMesherBackend did not execute the injected engine exactly once");
    if (result.mesh.nodes.size() != 3 || result.mesh.elements.size() != 1 ||
        result.mesh.periodicConstraints.size() != 1)
        return fail("TangleMesherBackend did not convert the engine result");
    if (receivedOptions.minimumAngleDegrees != 28.0 ||
        receivedOptions.maximumElementArea != 0.25 * 0.25 ||
        !receivedOptions.forceMaximumElementArea ||
        !receivedOptions.suppressExteriorSteinerPoints ||
        !receivedOptions.suppressUnusedVertices || !receivedOptions.verbose)
        return fail("TangleMesherBackend did not forward the meshing options");

    result = backend.mesh(problem, false, options);
    if (!result.succeeded() || calls != 2 || !result.mesh.periodicConstraints.empty())
        return fail("non-periodic Tangle request retained periodic constraints");

    // Invalid options are rejected before the engine runs.
    femm::mesh::MeshingOptions invalid;
    invalid.minimumAngleDegrees = -1.0;
    result = backend.mesh(problem, true, invalid);
    if (result.succeeded() || result.status != femm::mesh::MeshStatus::InvalidInput ||
        calls != 2 || result.diagnostics.empty())
        return fail("invalid meshing options were not rejected before engine execution");

    invalid = femm::mesh::MeshingOptions{};
    invalid.defaultElementSize = std::nan("");
    result = backend.mesh(problem, true, invalid);
    if (result.succeeded() || result.status != femm::mesh::MeshStatus::InvalidInput ||
        calls != 2)
        return fail("non-finite element size was not rejected before engine execution");

    fmesher::TangleMesherBackend failingBackend(
        [](const std::string &, const ::MeshOptions &, Mesh &) { return TANGLE_ERR_MESH; });
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

    if (std::string(backend.name()) != "Tangle")
        return fail("Tangle backend did not report its name");

    return 0;
}
