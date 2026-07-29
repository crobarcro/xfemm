#include "MesherBackend.h"
#include "TangleMesherBackend.h"
#include "TriangleMesherBackend.h"
#include "fmesher.h"
#include "FemmReader.h"

#include <chrono>
#include <cmath>
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

bool validNodeIndex(femm::mesh::MeshIndex index, const femm::mesh::SolverMesh &mesh)
{
    return index < mesh.nodes.size();
}

bool validGeometry(const femm::mesh::SolverMesh &mesh)
{
    if (mesh.nodes.empty() || mesh.elements.empty() || mesh.edges.empty()) return false;
    for (const auto &node : mesh.nodes)
        if (!std::isfinite(node.x) || !std::isfinite(node.y)) return false;
    for (const auto &element : mesh.elements)
        for (const auto node : element.nodes)
            if (!validNodeIndex(node, mesh)) return false;
    for (const auto &edge : mesh.edges)
        if (!validNodeIndex(edge.first, mesh) || !validNodeIndex(edge.second, mesh)) return false;
    return true;
}

bool validPeriodicTopology(const femm::mesh::SolverMesh &mesh)
{
    for (const auto &constraint : mesh.periodicConstraints)
        if (!validNodeIndex(constraint.first, mesh) || !validNodeIndex(constraint.second, mesh))
            return false;
    for (const auto &airGap : mesh.airGaps) {
        if (airGap.totalArcElements == 0 || airGap.quadraturePoints.empty()
                || airGap.nodeIndices.empty() || airGap.innerRing.empty()
                || airGap.outerRing.empty()
                || airGap.innerRing.size() != airGap.outerRing.size()) return false;
        for (const auto node : airGap.nodeIndices)
            if (!validNodeIndex(node, mesh)) return false;
        for (const auto &point : airGap.quadraturePoints)
            for (const auto node : point.nodes)
                if (!validNodeIndex(node, mesh)) return false;
        for (const auto &point : airGap.innerRing)
            if (!validNodeIndex(point.node, mesh)) return false;
        for (const auto &point : airGap.outerRing)
            if (!validNodeIndex(point.node, mesh)) return false;
    }
    return true;
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

    if (!result.succeeded() || !validGeometry(result.mesh))
        return fail("Non-periodic backend result was incomplete");
    if (!result.mesh.periodicConstraints.empty() || !result.mesh.airGaps.empty())
        return fail("Non-periodic backend result contained periodic topology");
    if (!periodicResult.succeeded() || !validGeometry(periodicResult.mesh)
            || periodicResult.mesh.periodicConstraints.empty()
            || !periodicResult.mesh.airGaps.empty()
            || !validPeriodicTopology(periodicResult.mesh))
        return fail("Ordinary periodic backend result was incomplete");
    if (!ageResult.succeeded() || !validGeometry(ageResult.mesh)
            || ageResult.mesh.airGaps.empty() || !validPeriodicTopology(ageResult.mesh))
        return fail("AGE backend result was incomplete");
    if (createdFile) return fail(backendName + "MesherBackend created a file");

    femm::mesh::SolverMesh::Node markerProbe;
    markerProbe.boundaryMarker = 7;
    // Marker decoding belongs to the solver consumer, not SolverMesh.
    const std::int32_t pointProperty = markerProbe.boundaryMarker > 1
            ? markerProbe.boundaryMarker - 2 : -1;
    return pointProperty == 5 ? 0 : fail("Raw boundary marker was not preserved");
}
