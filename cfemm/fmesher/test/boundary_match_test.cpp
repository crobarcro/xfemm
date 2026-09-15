#include "BoundaryMatchValidator.h"
#include "TangleMesherBackend.h"
#include "TriangleMesherBackend.h"

#include "FemmProblem.h"
#include "FemmReader.h"

#include <iostream>
#include <memory>
#include <stdexcept>
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

BoundaryMatch periodicMatch(std::size_t first, std::size_t second, bool field)
{
    BoundaryMatch match;
    match.first = {SourceEntityReference::Kind::Segment, first};
    match.second = {SourceEntityReference::Kind::Segment, second};
    match.createFieldConstraint = field;
    match.periodicity = SolverMesh::Periodicity::Periodic;
    return match;
}

} // namespace

int main(int argc, char **argv)
{
    if (argc != 3)
        return fail("expected periodic and antiperiodic FEMM fixture paths");
    auto problem = readProblem(argv[1]);

    // B2: boundary references are validated before any backend runs.
    MeshingRequest request;
    request.boundaryMatches.push_back(periodicMatch(0, 2, false));
    if (!fmesher::validateBoundaryMatches(*problem, request).valid())
        return fail("a valid boundary match was rejected");

    MeshingRequest mixed = request;
    mixed.boundaryMatches.back().second = {SourceEntityReference::Kind::Arc, 0};
    if (fmesher::validateBoundaryMatches(*problem, mixed).valid())
        return fail("mixed line/arc boundary match was accepted");

    MeshingRequest stale = request;
    stale.boundaryMatches.back().second.index = 99;
    if (fmesher::validateBoundaryMatches(*problem, stale).valid())
        return fail("stale boundary match reference was accepted");

    MeshingRequest nonPeriodic;
    nonPeriodic.boundaryMatches.push_back(periodicMatch(0, 1, false));
    if (fmesher::validateBoundaryMatches(*problem, nonPeriodic).valid())
        return fail("non-periodic boundary match was accepted");

    MeshingRequest duplicated = request;
    duplicated.boundaryMatches.push_back(periodicMatch(2, 0, false));
    if (fmesher::validateBoundaryMatches(*problem, duplicated).valid())
        return fail("duplicated boundary property was accepted");

    // B1: the legacy overload and an equivalent request agree.
    {
        fmesher::TangleMesherBackend backend;
        auto legacy = backend.mesh(*problem, true);
        MeshingRequest equivalent;
        equivalent.createPeriodicFieldConstraints = true;
        auto requested = backend.mesh(*problem, equivalent);
        if (!legacy.succeeded() || !requested.succeeded())
            return fail("legacy or request meshing failed");
        if (legacy.mesh.nodes.size() != requested.mesh.nodes.size() ||
            legacy.mesh.elements.size() != requested.mesh.elements.size() ||
            legacy.mesh.periodicConstraints.size() != requested.mesh.periodicConstraints.size())
            return fail("legacy and request meshing disagree");
    }

    // B3/B4: topology-only matching, then the same match with field constraints.
    MeshBoundaryMatch topologyMatch;
    {
        fmesher::TangleMesherBackend backend;
        auto topology = backend.mesh(*problem, request);
        if (!topology.succeeded())
            return fail("topology-only meshing failed");
        if (!topology.mesh.periodicConstraints.empty())
            return fail("topology-only request emitted periodic field constraints");
        if (topology.boundaryMatches.size() != 1)
            return fail("topology-only request did not return one boundary match");
        topologyMatch = topology.boundaryMatches.front();
        if (topologyMatch.boundaryProperty != 1)
            return fail("boundary match reported the wrong boundary property");
        if (topologyMatch.firstNodes.empty() ||
            topologyMatch.firstNodes.size() != topologyMatch.secondNodes.size())
            return fail("boundary match correspondence is empty or unequal");
        for (auto node : topologyMatch.firstNodes)
            if (node >= topology.mesh.nodes.size())
                return fail("boundary match first node is invalid");
        for (auto node : topologyMatch.secondNodes)
            if (node >= topology.mesh.nodes.size())
                return fail("boundary match second node is invalid");

        MeshingRequest withField = request;
        withField.boundaryMatches.front().createFieldConstraint = true;
        auto constrained = backend.mesh(*problem, withField);
        if (!constrained.succeeded() || constrained.mesh.periodicConstraints.empty())
            return fail("explicit field constraint request did not emit constraints");
        if (constrained.mesh.periodicConstraints.size() != topologyMatch.firstNodes.size())
            return fail("field constraint count does not match the correspondence");
    }

    // Antiperiodic matches keep their type through topology-only conversion.
    {
        auto antiperiodic = readProblem(argv[2]);
        MeshingRequest antiRequest;
        antiRequest.boundaryMatches.push_back(periodicMatch(0, 2, true));
        antiRequest.boundaryMatches.back().periodicity =
            SolverMesh::Periodicity::Antiperiodic;
        fmesher::TangleMesherBackend backend;
        auto result = backend.mesh(*antiperiodic, antiRequest);
        if (!result.succeeded() || result.boundaryMatches.size() != 1)
            return fail("antiperiodic topology-only meshing failed");
        if (result.boundaryMatches.front().periodicity !=
            SolverMesh::Periodicity::Antiperiodic)
            return fail("antiperiodic boundary match lost its periodicity");
        if (result.mesh.periodicConstraints.empty())
            return fail("antiperiodic field constraint request emitted no constraints");
        for (const auto &constraint : result.mesh.periodicConstraints)
            if (constraint.periodicity != SolverMesh::Periodicity::Antiperiodic)
                return fail("antiperiodic constraint has the wrong type");
    }

    // B5: Triangle reports a clear capability diagnostic instead of emulating
    // topology-only matching with periodic physics.
    {
        fmesher::TriangleMesherBackend triangle;
        auto unsupported = triangle.mesh(*problem, request);
        if (unsupported.succeeded() || unsupported.status != MeshStatus::Unsupported)
            return fail("Triangle did not report unsupported boundary matching");
        if (unsupported.diagnostics.empty() ||
            unsupported.diagnostics.front().backendName != "Triangle")
            return fail("Triangle capability diagnostic is missing or mislabelled");
    }

    return 0;
}
