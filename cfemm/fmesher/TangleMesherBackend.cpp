#include "TangleMesherBackend.h"

#include "BoundaryMatchValidator.h"
#include "TangleMeshConverter.h"

#include "CArcSegment.h"
#include "CBoundaryProp.h"
#include "CMaterialProp.h"
#include "CSegment.h"
#include "FemmProblem.h"
#include "femmconstants.h"

#include <tangle_mesh.h>

#include <algorithm>
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

void templateDiagnostic(std::vector<femm::mesh::MeshDiagnostic> &diagnostics,
                        const std::string &message)
{
    diagnostics.push_back(
        {femm::mesh::MeshDiagnosticSeverity::Error, message, "Tangle", 0});
}

bool rotationalChainsMatch(const femm::mesh::SolverMesh &mesh,
                           const std::vector<femm::mesh::MeshIndex> &from,
                           const std::vector<femm::mesh::MeshIndex> &to,
                           double centerX, double centerY, double angleDegrees)
{
    using namespace femm::mesh;
    if (from.empty() || from.size() != to.size())
        return false;
    const RigidTransform2D transform =
        RigidTransform2D::rotationAbout(centerX, centerY, angleDegrees);
    for (std::size_t i = 0; i < from.size(); ++i) {
        if (from[i] >= mesh.nodes.size() || to[i] >= mesh.nodes.size())
            return false;
        double x = 0.0, y = 0.0;
        transform.applyPoint(mesh.nodes[from[i]].x, mesh.nodes[from[i]].y, x, y);
        const auto &target = mesh.nodes[to[i]];
        if (std::abs(x - target.x) > InstancedMeshWeldToleranceMetres ||
            std::abs(y - target.y) > InstancedMeshWeldToleranceMetres)
            return false;
    }
    return true;
}

struct TemplateBuild {
    bool ok = false;
    femm::mesh::InstancedMesh instanced;
    std::vector<femm::mesh::MeshDiagnostic> diagnostics;
};

TemplateBuild buildRotationalTemplate(const femm::FemmProblem &problem,
                                      const femm::mesh::MeshResult &tile,
                                      const femm::mesh::TemplateRequest &request)
{
    using namespace femm::mesh;
    TemplateBuild build;
    const auto fail = [&](const std::string &message) {
        templateDiagnostic(build.diagnostics, message);
        return build;
    };

    if (!std::isfinite(request.centerXMetres) || !std::isfinite(request.centerYMetres))
        return fail("rotational template centre must be finite");
    if (request.instanceCount < 2)
        return fail("a rotational template requires at least two instances");
    if (!std::isfinite(request.totalAngleDegrees) ||
        std::abs(request.totalAngleDegrees - 360.0) > 1e-9)
        return fail("a closed rotational template must cover 360 degrees");
    const double angleStep =
        request.totalAngleDegrees / static_cast<double>(request.instanceCount);
    if (!(angleStep > 0.0) || angleStep >= 360.0)
        return fail("rotational template instances overlap");
    if (!tile.mesh.airGaps.empty())
        return fail("an AGE inside a rotational template is not supported");
    if (tile.boundaryMatches.size() != 1)
        return fail("a rotational template requires exactly one matched seam pair");
    for (const auto &property : problem.blockproplist) {
        const auto *material = dynamic_cast<const femm::CMMaterialProp *>(property.get());
        if (material && material->mu_x != material->mu_y)
            return fail("anisotropic materials are not supported in a rotational template");
    }

    const auto &match = tile.boundaryMatches.front();
    if (match.firstNodes.empty() || match.firstNodes.size() != match.secondNodes.size())
        return fail("matched template seams have incompatible node counts");
    if (!request.seamBoundaryProperties.empty() &&
        std::find(request.seamBoundaryProperties.begin(),
                  request.seamBoundaryProperties.end(),
                  match.boundaryProperty) == request.seamBoundaryProperties.end())
        return fail("requested seam boundary property is not the matched template seam");

    MeshTemplate meshTemplate;
    meshTemplate.localMesh = tile.mesh;
    meshTemplate.localMesh.periodicConstraints.clear();
    meshTemplate.localMesh.airGaps.clear();
    meshTemplate.seams.push_back({"first", match.firstNodes});
    meshTemplate.seams.push_back({"second", match.secondNodes});

    // The leading seam is the one a positive step rotation maps onto the other.
    std::size_t leadSeam = 0;
    std::size_t trailSeam = 1;
    if (!rotationalChainsMatch(tile.mesh, match.firstNodes, match.secondNodes,
                               request.centerXMetres, request.centerYMetres, angleStep)) {
        if (rotationalChainsMatch(tile.mesh, match.secondNodes, match.firstNodes,
                                  request.centerXMetres, request.centerYMetres, angleStep)) {
            leadSeam = 1;
            trailSeam = 0;
        } else {
            return fail("matched template seams are not related by the declared rotation");
        }
    }

    build.instanced.templates.push_back(std::move(meshTemplate));
    for (std::size_t k = 0; k < request.instanceCount; ++k) {
        MeshInstance instance;
        instance.templateIndex = 0;
        instance.transform = RigidTransform2D::rotationAbout(
            request.centerXMetres, request.centerYMetres,
            static_cast<double>(k) * angleStep);
        instance.seamConnections.push_back(
            {trailSeam, (k + 1) % request.instanceCount, leadSeam,
             SeamOrientation::Forward});
        build.instanced.instances.push_back(std::move(instance));
    }
    build.ok = true;
    return build;
}

/**
 * Convert xfemm's in-memory problem into Tangle's FEMM-like record set. xfemm
 * stores zero-based boundary/circuit references (with -1 for none) and areas,
 * while the record set follows the .fem columns (one-based references, 0 for
 * none, and a mesh-size diameter), so the conversion normalises both.
 */
::FemProblem toTangleProblem(const femm::FemmProblem &problem)
{
    constexpr double Pi = 3.141592653589793238462643383279502884;
    ::FemProblem out;
    out.isMagnetics = problem.filetype == femm::FileType::MagneticsFile;
    out.minAngle = problem.MinAngle;
    out.doSmartMesh = problem.DoSmartMesh;

    for (const auto &node : problem.nodelist) {
        if (!node)
            continue;
        ::FemProblem::Node converted;
        converted.x = node->x;
        converted.y = node->y;
        converted.boundaryMarker = node->BoundaryMarker >= 0 ? node->BoundaryMarker + 1 : 0;
        converted.group = node->InGroup;
        out.nodes.push_back(converted);
    }
    for (const auto &segment : problem.linelist) {
        if (!segment)
            continue;
        ::FemProblem::Segment converted;
        converted.n0 = segment->n0;
        converted.n1 = segment->n1;
        converted.maxSideLength = segment->MaxSideLength;
        converted.boundaryMarker =
            segment->BoundaryMarker >= 0 ? segment->BoundaryMarker + 1 : 0;
        converted.hidden = segment->Hidden ? 1 : 0;
        converted.group = segment->InGroup;
        out.segments.push_back(converted);
    }
    for (const auto &arc : problem.arclist) {
        if (!arc)
            continue;
        ::FemProblem::Arc converted;
        converted.n0 = arc->n0;
        converted.n1 = arc->n1;
        converted.arcLength = arc->ArcLength;
        converted.maxSegDegrees = arc->MaxSideLength;
        converted.boundaryMarker = arc->BoundaryMarker >= 0 ? arc->BoundaryMarker + 1 : 0;
        converted.hidden = arc->Hidden ? 1 : 0;
        converted.group = arc->InGroup;
        out.arcs.push_back(converted);
    }
    for (const auto &property : problem.lineproplist) {
        ::FemProblem::Boundary converted;
        if (property) {
            converted.name = property->BdryName;
            converted.format = property->BdryFormat;
            converted.innerAngle = property->InnerAngle;
            converted.outerAngle = property->OuterAngle;
        }
        out.boundaries.push_back(converted);
    }
    for (const auto &labelPtr : problem.labellist) {
        const auto *label = dynamic_cast<const femm::CMBlockLabel *>(labelPtr.get());
        if (!label)
            continue;
        if (label->isHole()) {
            out.holes.push_back({label->x, label->y});
            continue;
        }
        ::FemProblem::Label converted;
        converted.x = label->x;
        converted.y = label->y;
        converted.blockType = label->BlockType + 1;
        converted.maxAreaDiameter =
            label->MaxArea > 0 ? std::sqrt(4.0 * label->MaxArea / Pi) : -1.0;
        converted.inCircuit = label->InCircuit >= 0 ? label->InCircuit + 1 : 0;
        converted.magDir = label->MagDir;
        converted.group = label->InGroup;
        converted.turns = label->Turns;
        converted.isExternal = label->IsExternal ? 1 : 0;
        out.labels.push_back(converted);
    }
    return out;
}

std::vector<femm::mesh::MeshDiagnostic> materializationDiagnostics(
    const std::vector<femm::mesh::MaterializationDiagnostic> &diagnostics)
{
    std::vector<femm::mesh::MeshDiagnostic> converted;
    converted.reserve(diagnostics.size());
    for (const auto &diagnostic : diagnostics)
        converted.push_back({femm::mesh::MeshDiagnosticSeverity::Error,
                             diagnostic.message, "Tangle",
                             static_cast<int>(diagnostic.category)});
    return converted;
}

} // namespace

TangleMesherBackend::TangleMesherBackend(Engine engine, InMemoryEngine inMemoryEngine)
    : engine_(engine ? std::move(engine)
                     : Engine([](const std::string &path, const ::MeshOptions &options,
                                 ::Mesh &mesh) {
                           return tangle_mesh_fem(path, options, mesh);
                       }))
    , inMemoryEngine_(inMemoryEngine
                          ? std::move(inMemoryEngine)
                          : InMemoryEngine([](const ::FemProblem &problem,
                                              const ::MeshOptions &options, ::Mesh &mesh) {
                                return tangle_mesh_fem(problem, options, mesh);
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
    // A titled problem is meshed from its FEMM file; a pathless problem is a
    // tile or model held in memory and is meshed through the record API.
    const bool inMemory = problem.getTitle().empty();

    const bool instancing = !request.templates.empty();
    if (instancing && request.templates.size() > 1) {
        result.status = femm::mesh::MeshStatus::Unsupported;
        result.diagnostics.push_back(
            {femm::mesh::MeshDiagnosticSeverity::Error,
             "only one rotational template is supported", "Tangle", 0});
        return result;
    }

    // A template tile is meshed with matched seams as topology only; the tile's
    // own periodic declarations are not field constraints.
    femm::mesh::MeshingRequest tileRequest = request;
    if (instancing) {
        tileRequest.createPeriodicFieldConstraints = false;
        tileRequest.templates.clear();
    }

    const auto matchValidation = validateBoundaryMatches(problem, tileRequest);
    if (!matchValidation.valid()) {
        result.status = femm::mesh::MeshStatus::InvalidInput;
        result.diagnostics = matchValidation.diagnostics;
        return result;
    }

    ::Mesh tangleMesh;
    const int engineStatus =
        inMemory
            ? inMemoryEngine_(toTangleProblem(problem), translateOptions(tileRequest.options),
                              tangleMesh)
            : engine_(problem.getTitle(), translateOptions(tileRequest.options), tangleMesh);
    if (engineStatus != TANGLE_OK)
        return failure(statusFor(engineStatus), statusMessage(engineStatus), engineStatus);

    auto converted =
        convertTangleMesh(tangleMesh, femm::LengthConvMeters[problem.LengthUnits]);
    if (!converted.succeeded())
        return converted;

    if (!instancing) {
        applyConstraintPolicy(converted, problem, request);
        return converted;
    }

    auto build = buildRotationalTemplate(problem, converted, request.templates.front());
    if (!build.ok) {
        converted.status = femm::mesh::MeshStatus::InvalidInput;
        converted.diagnostics = std::move(build.diagnostics);
        return converted;
    }
    auto materialized = build.instanced.materialize();
    if (!materialized.succeeded()) {
        converted.status = femm::mesh::MeshStatus::InvalidInput;
        converted.diagnostics = materializationDiagnostics(materialized.diagnostics);
        return converted;
    }
    converted.mesh = std::move(materialized.mesh);
    converted.instancing = std::move(materialized.provenance);
    converted.instancedTemplates = std::move(build.instanced);
    converted.boundaryMatches.clear();
    return converted;
}

} // namespace fmesher
