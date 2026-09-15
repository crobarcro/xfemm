#include "TangleMesherBackend.h"

#include "BoundaryMatchValidator.h"
#include "TangleMeshConverter.h"

#include "CArcSegment.h"
#include "CBoundaryProp.h"
#include "CSegment.h"
#include "FemmProblem.h"
#include "femmconstants.h"

#include <tangle_mesh.h>

#include <cmath>
#include <limits>
#include <set>
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

bool validateOptions(const femm::mesh::MeshingOptions &options,
                     std::vector<femm::mesh::MeshDiagnostic> &diagnostics)
{
    bool valid = true;
    const auto error = [&](const std::string &message) {
        valid = false;
        diagnostics.push_back(
            {femm::mesh::MeshDiagnosticSeverity::Error, message, "Tangle", 0});
    };
    if (!std::isfinite(options.minimumAngleDegrees) || options.minimumAngleDegrees < 0.0)
        error("Tangle requires a finite, nonnegative minimum angle");
    if (!std::isfinite(options.defaultElementSize) || options.defaultElementSize < 0.0)
        error("Tangle requires a finite, nonnegative default element size");
    return valid;
}

::MeshOptions translateOptions(const femm::mesh::MeshingOptions &options)
{
    ::MeshOptions translated;
    translated.minimumAngleDegrees = options.minimumAngleDegrees;
    if (options.defaultElementSize > 0.0) {
        // xfemm expresses the control as a length; Tangle applies an area limit.
        translated.maximumElementArea = options.defaultElementSize * options.defaultElementSize;
    }
    translated.forceMaximumElementArea = options.forceMaximumElementArea;
    translated.suppressExteriorSteinerPoints = options.suppressExteriorSteinerPoints;
    translated.suppressUnusedVertices = options.suppressUnusedVertices;
    translated.verbose = options.verbose;
    return translated;
}

std::size_t matchBoundaryProperty(const femm::FemmProblem &problem,
                                  const femm::mesh::BoundaryMatch &match)
{
    using Kind = femm::mesh::SourceEntityReference::Kind;
    if (match.first.kind == Kind::Segment && match.first.index < problem.linelist.size()) {
        const auto *segment = problem.linelist[match.first.index].get();
        if (segment && segment->hasBoundaryMarker())
            return static_cast<std::size_t>(segment->BoundaryMarker);
    } else if (match.first.kind == Kind::Arc && match.first.index < problem.arclist.size()) {
        const auto *arc = problem.arclist[match.first.index].get();
        if (arc && arc->hasBoundaryMarker())
            return static_cast<std::size_t>(arc->BoundaryMarker);
    }
    return std::numeric_limits<std::size_t>::max();
}

void applyConstraintPolicy(femm::mesh::MeshResult &result,
                           const femm::FemmProblem &problem,
                           const femm::mesh::MeshingRequest &request)
{
    using femm::mesh::SolverMesh;
    if (request.boundaryMatches.empty()) {
        if (!request.createPeriodicFieldConstraints)
            result.mesh.periodicConstraints.clear();
        return;
    }

    std::set<std::size_t> emitting;
    for (const auto &match : request.boundaryMatches)
        if (match.createFieldConstraint)
            emitting.insert(matchBoundaryProperty(problem, match));

    std::set<std::size_t> produced;
    for (const auto &match : result.boundaryMatches)
        produced.insert(match.boundaryProperty);

    std::vector<SolverMesh::PeriodicConstraint> constraints;
    for (const auto &match : result.boundaryMatches) {
        if (emitting.count(match.boundaryProperty) == 0)
            continue;
        const std::size_t count =
            std::min(match.firstNodes.size(), match.secondNodes.size());
        for (std::size_t i = 0; i < count; ++i)
            constraints.push_back(
                {match.firstNodes[i], match.secondNodes[i], match.periodicity});
    }
    result.mesh.periodicConstraints = std::move(constraints);

    for (std::size_t property : emitting)
        if (produced.count(property) == 0)
            result.diagnostics.push_back(
                {femm::mesh::MeshDiagnosticSeverity::Warning,
                 "requested boundary match was not produced by the mesher; "
                 "the problem may not declare the pair as periodic",
                 "Tangle", 0});
}

} // namespace

TangleMesherBackend::TangleMesherBackend(Engine engine)
    : engine_(engine ? std::move(engine)
                     : Engine([](const std::string &path, const ::MeshOptions &options,
                                 ::Mesh &mesh) {
                           return tangle_mesh_fem(path, options, mesh);
                       }))
{
}

femm::mesh::MeshResult TangleMesherBackend::mesh(
        femm::FemmProblem &problem, const femm::mesh::MeshingRequest &request)
{
    femm::mesh::MeshResult result;
    if (!validateOptions(request.options, result.diagnostics)) {
        result.status = femm::mesh::MeshStatus::InvalidInput;
        return result;
    }
    if (problem.getTitle().empty())
        return failure(femm::mesh::MeshStatus::InvalidInput,
                       "Tangle requires a FemmProblem loaded from a FEMM file",
                       TANGLE_ERR_NO_FILE);

    const auto matchValidation = validateBoundaryMatches(problem, request);
    if (!matchValidation.valid()) {
        result.status = femm::mesh::MeshStatus::InvalidInput;
        result.diagnostics = matchValidation.diagnostics;
        return result;
    }

    ::Mesh tangleMesh;
    const int engineStatus =
        engine_(problem.getTitle(), translateOptions(request.options), tangleMesh);
    if (engineStatus != TANGLE_OK)
        return failure(statusFor(engineStatus), statusMessage(engineStatus), engineStatus);

    auto converted =
        convertTangleMesh(tangleMesh, femm::LengthConvMeters[problem.LengthUnits]);
    if (!converted.succeeded())
        return converted;

    applyConstraintPolicy(converted, problem, request);
    return converted;
}

} // namespace fmesher
