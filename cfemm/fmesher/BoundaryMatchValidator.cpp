#include "BoundaryMatchValidator.h"

#include "CArcSegment.h"
#include "CBoundaryProp.h"
#include "CSegment.h"
#include "FemmProblem.h"

#include <set>
#include <string>

namespace fmesher {
namespace {

using femm::mesh::MeshDiagnostic;
using femm::mesh::MeshDiagnosticSeverity;
using femm::mesh::SourceEntityReference;

void error(std::vector<MeshDiagnostic> &diagnostics, const std::string &message)
{
    diagnostics.push_back({MeshDiagnosticSeverity::Error, message, "BoundaryMatch", 0});
}

bool periodicProperty(const femm::FemmProblem &problem, int property,
                      femm::mesh::SolverMesh::Periodicity periodicity)
{
    if (property < 0 || static_cast<std::size_t>(property) >= problem.lineproplist.size())
        return false;
    const auto *boundary = problem.lineproplist[static_cast<std::size_t>(property)].get();
    if (!boundary)
        return false;
    using PeriodicityType = femm::CBoundaryProp::PeriodicityType;
    const PeriodicityType requested =
        periodicity == femm::mesh::SolverMesh::Periodicity::Antiperiodic
            ? PeriodicityType::AntiPeriodic
            : PeriodicityType::Periodic;
    return boundary->isPeriodic(requested);
}

std::size_t countSegments(const femm::FemmProblem &problem, int property)
{
    std::size_t count = 0;
    for (const auto &segment : problem.linelist)
        if (segment && segment->hasBoundaryMarker() && segment->BoundaryMarker == property)
            ++count;
    return count;
}

std::size_t countArcs(const femm::FemmProblem &problem, int property)
{
    std::size_t count = 0;
    for (const auto &arc : problem.arclist)
        if (arc && arc->hasBoundaryMarker() && arc->BoundaryMarker == property)
            ++count;
    return count;
}

struct ResolvedReference {
    bool valid = false;
    int boundaryProperty = -1;
    bool isArc = false;
};

ResolvedReference resolve(const femm::FemmProblem &problem,
                          const SourceEntityReference &reference)
{
    if (reference.kind == SourceEntityReference::Kind::Segment) {
        if (reference.index >= problem.linelist.size())
            return {};
        const auto *segment = problem.linelist[reference.index].get();
        if (!segment || !segment->hasBoundaryMarker())
            return {};
        return {true, segment->BoundaryMarker, false};
    }
    if (reference.index >= problem.arclist.size())
        return {};
    const auto *arc = problem.arclist[reference.index].get();
    if (!arc || !arc->hasBoundaryMarker())
        return {};
    return {true, arc->BoundaryMarker, true};
}

} // namespace

BoundaryMatchValidation validateBoundaryMatches(
    const femm::FemmProblem &problem,
    const femm::mesh::MeshingRequest &request)
{
    BoundaryMatchValidation result;
    std::set<std::size_t> usedProperties;
    for (std::size_t index = 0; index < request.boundaryMatches.size(); ++index) {
        const auto &match = request.boundaryMatches[index];
        const std::string prefix = "boundary match " + std::to_string(index) + ": ";
        if (match.first.kind != match.second.kind) {
            error(result.diagnostics, prefix + "cannot pair a line segment with an arc");
            continue;
        }
        if (match.first.kind == match.second.kind &&
            match.first.index == match.second.index) {
            error(result.diagnostics, prefix + "cannot match an entity with itself");
            continue;
        }
        const ResolvedReference first = resolve(problem, match.first);
        const ResolvedReference second = resolve(problem, match.second);
        if (!first.valid || !second.valid) {
            error(result.diagnostics,
                  prefix + "references a missing or markerless source entity");
            continue;
        }
        if (first.boundaryProperty != second.boundaryProperty) {
            error(result.diagnostics,
                  prefix + "source entities do not share one boundary property");
            continue;
        }
        const int property = first.boundaryProperty;
        if (!periodicProperty(problem, property, match.periodicity)) {
            error(result.diagnostics,
                  prefix + "boundary property does not declare the requested periodicity");
            continue;
        }
        if (!usedProperties.insert(static_cast<std::size_t>(property)).second) {
            error(result.diagnostics, prefix + "boundary property is matched more than once");
            continue;
        }
        const std::size_t count = first.isArc ? countArcs(problem, property)
                                              : countSegments(problem, property);
        if (count != 2) {
            error(result.diagnostics,
                  prefix + "boundary property must name exactly two " +
                      std::string(first.isArc ? "arcs" : "segments") + ", found " +
                      std::to_string(count));
        }
    }
    return result;
}

} // namespace fmesher
