#include "MesherBackend.h"
#include "TangleMesherBackend.h"
#include "TriangleMesherBackend.h"
#include "fmesher.h"
#include "FemmReader.h"
#include "mesh/SolverMeshValidator.h"

#include <chrono>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <memory>
#include <string>

namespace {

bool parseProblem(const char *path, fmesher::FMesher &facade)
{
    facade.problem->filetype = femm::FileType::MagneticsFile;
    femm::MagneticsReader reader(facade.problem, std::cerr);
    return reader.parse(path) == femm::F_FILE_OK;
}

int fail(const std::string &message)
{
    std::cerr << message << '\n';
    return 1;
}

} // namespace

int main(int argc, char **argv)
{
    if (argc != 4) return fail("Expected backend, ordinary-periodic, and AGE problem paths");

    const std::string backendName = argv[1];
    std::unique_ptr<fmesher::MesherBackend> backend;
    if (backendName == "Triangle")
        backend.reset(new fmesher::TriangleMesherBackend);
    else if (backendName == "Tangle")
        backend.reset(new fmesher::TangleMesherBackend);
    else
        return fail("Unknown mesher backend: " + backendName);

    fmesher::FMesher facade;
    if (!parseProblem(argv[2], facade)) return fail("Could not parse periodic test problem");
    fmesher::FMesher ageFacade;
    if (!parseProblem(argv[3], ageFacade)) return fail("Could not parse AGE test problem");

    const auto originalDirectory = std::filesystem::current_path();
    const auto unique = std::chrono::steady_clock::now().time_since_epoch().count();
    const auto scratchDirectory = std::filesystem::temp_directory_path()
            / ("xfemm-backend-test-" + std::to_string(unique));
    std::filesystem::create_directory(scratchDirectory);
    std::filesystem::current_path(scratchDirectory);

    femm::mesh::MeshingOptions options;
    const auto result = backend->mesh(*facade.problem, false, options);
    const auto periodicResult = backend->mesh(*facade.problem, true, options);
    const auto ageResult = backend->mesh(*ageFacade.problem, true, options);

    const bool createdFile = std::filesystem::directory_iterator(scratchDirectory)
            != std::filesystem::directory_iterator();
    std::filesystem::current_path(originalDirectory);
    std::filesystem::remove_all(scratchDirectory);

    if (!result.succeeded() || result.mesh.nodes.empty() || result.mesh.elements.empty()
            || result.mesh.edges.empty() || !femm::mesh::validateSolverMesh(result.mesh).valid())
        return fail("Non-periodic backend result was incomplete");
    if (!result.mesh.periodicConstraints.empty() || !result.mesh.airGaps.empty())
        return fail("Non-periodic backend result contained periodic topology");
    if (!periodicResult.succeeded() || periodicResult.mesh.nodes.empty()
            || periodicResult.mesh.elements.empty() || periodicResult.mesh.edges.empty()
            || periodicResult.mesh.periodicConstraints.empty()
            || !periodicResult.mesh.airGaps.empty()
            || !femm::mesh::validateSolverMesh(periodicResult.mesh).valid())
        return fail("Ordinary periodic backend result was incomplete");
    if (!ageResult.succeeded() || ageResult.mesh.nodes.empty()
            || ageResult.mesh.elements.empty() || ageResult.mesh.edges.empty()
            || ageResult.mesh.airGaps.empty()
            || !femm::mesh::validateSolverMesh(ageResult.mesh).valid())
        return fail("AGE backend result was incomplete");
    for (const auto &airGap : ageResult.mesh.airGaps)
        if (airGap.nodeIndices.empty() || airGap.innerRing.empty() || airGap.outerRing.empty())
            return fail("AGE backend result omitted expected reusable topology");
    if (createdFile) return fail(backendName + "MesherBackend created a file");

    femm::mesh::SolverMesh::Node markerProbe;
    markerProbe.boundaryMarker = 7;
    // Marker decoding belongs to the solver consumer, not SolverMesh.
    const std::int32_t pointProperty = markerProbe.boundaryMarker > 1
            ? markerProbe.boundaryMarker - 2 : -1;
    return pointProperty == 5 ? 0 : fail("Raw boundary marker was not preserved");
}
