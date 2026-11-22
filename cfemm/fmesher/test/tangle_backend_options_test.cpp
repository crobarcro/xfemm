#include "TangleMesherBackend.h"

#include "FemmProblem.h"
#include "FemmReader.h"
#include "mesh/SolverMeshValidator.h"

#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

namespace {

using femm::mesh::MeshingOptions;

int fail(const std::string &message)
{
    std::cerr << message << '\n';
    return 1;
}

std::unique_ptr<femm::FemmProblem> readProblem(const std::string &path)
{
    auto problem = std::make_unique<femm::FemmProblem>(femm::FileType::MagneticsFile);
    auto view = std::shared_ptr<femm::FemmProblem>(problem.get(), [](femm::FemmProblem *) {});
    femm::MagneticsReader reader(view, std::cerr);
    if (reader.parse(path) != femm::F_FILE_OK)
        throw std::runtime_error("could not parse " + path);
    return problem;
}

} // namespace

int main(int argc, char **argv)
{
    if (argc != 2)
        return fail("expected a non-periodic FEMM fixture path");
    auto problem = readProblem(argv[1]);
    fmesher::TangleMesherBackend backend;

    const auto base = backend.mesh(*problem, false);
    if (!base.succeeded() || !femm::mesh::validateSolverMesh(base.mesh).valid())
        return fail("default Tangle meshing failed");

    MeshingOptions finerAngle;
    finerAngle.minimumAngleDegrees = 33.0;
    const auto angle = backend.mesh(*problem, false, finerAngle);
    if (!angle.succeeded() || !femm::mesh::validateSolverMesh(angle.mesh).valid())
        return fail("minimum-angle Tangle meshing failed");
    if (angle.mesh.elements.size() <= base.mesh.elements.size())
        return fail("minimum angle did not refine the mesh");

    MeshingOptions finerSize;
    finerSize.defaultElementSize = 0.02;
    const auto size = backend.mesh(*problem, false, finerSize);
    if (!size.succeeded() || !femm::mesh::validateSolverMesh(size.mesh).valid())
        return fail("element-size Tangle meshing failed");
    if (size.mesh.elements.size() <= base.mesh.elements.size())
        return fail("default element size did not refine the mesh");

    MeshingOptions forcedArea;
    forcedArea.defaultElementSize = 0.5;
    forcedArea.forceMaximumElementArea = true;
    const auto forced = backend.mesh(*problem, false, forcedArea);
    if (!forced.succeeded() || !femm::mesh::validateSolverMesh(forced.mesh).valid())
        return fail("forced-area Tangle meshing failed");
    if (forced.mesh.elements.size() >= base.mesh.elements.size())
        return fail("forced maximum element area did not coarsen the mesh");

    MeshingOptions noSteiner;
    noSteiner.suppressExteriorSteinerPoints = true;
    const auto suppressed = backend.mesh(*problem, false, noSteiner);
    if (!suppressed.succeeded() || !femm::mesh::validateSolverMesh(suppressed.mesh).valid())
        return fail("boundary Steiner suppression produced an invalid mesh");
    if (suppressed.mesh.elements.size() >= base.mesh.elements.size())
        return fail("boundary Steiner suppression did not reduce refinement");

    MeshingOptions verbose;
    verbose.verbose = true;
    const auto spoken = backend.mesh(*problem, false, verbose);
    if (!spoken.succeeded() || spoken.mesh.elements.size() != base.mesh.elements.size())
        return fail("verbose option changed the mesh");

    return 0;
}
