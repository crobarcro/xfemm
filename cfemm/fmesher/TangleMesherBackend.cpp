#include "TangleMesherBackend.h"

#include "TangleMeshConverter.h"

#include "FemmProblem.h"
#include "femmconstants.h"

#include <tangle_mesh.h>

#include <utility>

namespace fmesher {
namespace {

femm::mesh::MeshResult failure(femm::mesh::MeshStatus status,
                               const std::string &message, int engineStatus)
{
    femm::mesh::MeshResult result;
    result.status = status;
    result.diagnostics.push_back(
        {femm::mesh::MeshDiagnosticSeverity::Error, message, "Tangle", engineStatus});
    return result;
}

femm::mesh::MeshStatus statusFor(int engineStatus)
{
    switch (engineStatus) {
    case TANGLE_ERR_USAGE:
    case TANGLE_ERR_NO_FILE:
    case TANGLE_ERR_PARSE:
    case TANGLE_ERR_OPTION:
        return femm::mesh::MeshStatus::InvalidInput;
    default:
        return femm::mesh::MeshStatus::BackendFailure;
    }
}

const char *statusMessage(int engineStatus)
{
    switch (engineStatus) {
    case TANGLE_ERR_USAGE: return "Tangle rejected the meshing request";
    case TANGLE_ERR_NO_FILE: return "Tangle could not open the source FEMM file";
    case TANGLE_ERR_PARSE: return "Tangle could not parse the source FEMM file";
    case TANGLE_ERR_MESH: return "Tangle could not generate a conforming mesh";
    case TANGLE_ERR_OPTION: return "Tangle rejected an option for this FEMM problem";
    default: return "Tangle returned an unknown meshing error";
    }
}

} // namespace

TangleMesherBackend::TangleMesherBackend(Engine engine)
    : engine_(engine ? std::move(engine) : Engine(tangle_mesh_fem))
{
}

femm::mesh::MeshResult TangleMesherBackend::mesh(
        femm::FemmProblem &problem, bool periodic,
        const femm::mesh::MeshingOptions &options)
{
    (void)options; // Option mapping is A2.3; Tangle defaults are used for now.
    if (problem.getTitle().empty())
        return failure(femm::mesh::MeshStatus::InvalidInput,
                       "Tangle requires a FemmProblem loaded from a FEMM file",
                       TANGLE_ERR_NO_FILE);

    ::Mesh tangleMesh;
    const int engineStatus = engine_(problem.getTitle(), tangleMesh);
    if (engineStatus != TANGLE_OK)
        return failure(statusFor(engineStatus), statusMessage(engineStatus), engineStatus);

    auto result = convertTangleMesh(
        tangleMesh, femm::LengthConvMeters[problem.LengthUnits]);
    if (!periodic && result.succeeded()) {
        result.mesh.periodicConstraints.clear();
        result.mesh.airGaps.clear();
    }
    return result;
}

} // namespace fmesher
