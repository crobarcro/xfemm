#include "TangleMesherBackend.h"
#include "TriangleMesherBackend.h"

#include "CMaterialProp.h"
#include "FemmProblem.h"
#include "FemmReader.h"
#include "mesh/SolverMeshValidator.h"

#include <cmath>
#include <iostream>
#include <memory>
#include <string>

namespace {

using namespace femm::mesh;

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

bool hasError(const std::vector<MeshDiagnostic> &diagnostics)
{
    for (const auto &diagnostic : diagnostics)
        if (diagnostic.severity == MeshDiagnosticSeverity::Error)
            return true;
    return false;
}

TemplateRequest validRequest()
{
    TemplateRequest request;
    request.instanceCount = 8;
    request.totalAngleDegrees = 360.0;
    return request;
}

bool expectRejected(fmesher::TangleMesherBackend &backend, femm::FemmProblem &problem,
                    const MeshingRequest &request, const std::string &label)
{
    auto result = backend.mesh(problem, request);
    if (result.succeeded() || !hasError(result.diagnostics)) {
        std::cerr << label << ": expected an error diagnostic\n";
        return false;
    }
    return true;
}

} // namespace

int main(int argc, char **argv)
{
    if (argc != 3)
        return fail("expected tile and AGE fixture paths");
    try {
        auto tile = readProblem(argv[1]);
        fmesher::TangleMesherBackend tangle;

        // A valid request materialises the whole ring and exposes provenance.
        MeshingRequest request;
        request.templates.push_back(validRequest());
        auto result = tangle.mesh(*tile, request);
        if (!result.succeeded() || !result.instancing)
            return fail("valid template request did not produce instancing provenance");
        if (result.instancing->nodeMap.size() != 8)
            return fail("valid template request produced the wrong instance count");
        if (result.mesh.nodes.empty() || !validateSolverMesh(result.mesh).valid())
            return fail("valid template request produced an invalid materialised mesh");

        // Incomplete coverage, overlap, non-finite centre, and bad seam selection.
        MeshingRequest tooFew;
        tooFew.templates.push_back(validRequest());
        tooFew.templates.front().instanceCount = 1;
        if (!expectRejected(tangle, *tile, tooFew, "single instance"))
            return 1;

        MeshingRequest partial;
        partial.templates.push_back(validRequest());
        partial.templates.front().totalAngleDegrees = 180.0;
        if (!expectRejected(tangle, *tile, partial, "incomplete coverage"))
            return 1;

        MeshingRequest nonFinite;
        nonFinite.templates.push_back(validRequest());
        nonFinite.templates.front().centerXMetres = std::nan("");
        if (!expectRejected(tangle, *tile, nonFinite, "non-finite centre"))
            return 1;

        MeshingRequest badSeam;
        badSeam.templates.push_back(validRequest());
        badSeam.templates.front().seamBoundaryProperties = {99};
        if (!expectRejected(tangle, *tile, badSeam, "unknown seam boundary"))
            return 1;

        MeshingRequest twoTemplates;
        twoTemplates.templates.push_back(validRequest());
        twoTemplates.templates.push_back(validRequest());
        auto multiple = tangle.mesh(*tile, twoTemplates);
        if (multiple.succeeded() || multiple.status != MeshStatus::Unsupported)
            return fail("multiple templates were not reported as unsupported");

        // Anisotropic materials are rejected.
        auto anisotropic = readProblem(argv[1]);
        for (auto &property : anisotropic->blockproplist) {
            auto *material = dynamic_cast<femm::CMMaterialProp *>(property.get());
            if (material)
                material->mu_y = 2.0;
        }
        if (!expectRejected(tangle, *anisotropic, request, "anisotropic material"))
            return 1;

        // An AGE inside the tile is rejected.
        auto age = readProblem(argv[2]);
        if (!expectRejected(tangle, *age, request, "AGE inside template"))
            return 1;

        // Triangle reports the capability gap rather than emulating it.
        fmesher::TriangleMesherBackend triangle;
        auto unsupported = triangle.mesh(*tile, request);
        if (unsupported.succeeded() || unsupported.status != MeshStatus::Unsupported ||
            unsupported.diagnostics.empty() ||
            unsupported.diagnostics.front().backendName != "Triangle")
            return fail("Triangle did not report rotational templates as unsupported");
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
