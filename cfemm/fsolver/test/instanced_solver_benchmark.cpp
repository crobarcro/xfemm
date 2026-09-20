/*
   Instanced-vs-conventional magnetostatic solver benchmark.

   Three methods are available:

     conventional  mesh the full .fem geometry with Tangle and solve it;
     materialized  mesh each tile once, expand the InstancedMesh into a
                   SolverMesh, and solve that;
     native        mesh each tile once, build a compact LogicalMeshView, and
                   solve through the compressed native assembly.

   Run one method per process (the shell runner does this) so peak RSS is not
   contaminated by an earlier method. Each row reports topology counts, stored
   and expanded mesh bytes, meshing/view/solve wall time, a field checksum for
   sanity, and the process peak resident set size.
*/

#include "AnalysisSession.h"
#include "FSolverAnalysisBackend.h"
#include "FemmReader.h"
#include "TangleMesherBackend.h"
#include "TiledModel.h"
#include "TiledModelJson.h"
#include "TiledModelMesher.h"
#include "mesh/LogicalMeshView.h"

#include <sys/resource.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;
using femm::AnalysisSession;
using femm::FSolverAnalysisBackend;
using femm::ModelDefinition;
using femm::TrialSolution;

struct Options {
    std::string method = "all";
    std::string tiled;
    std::string redraw;
    int repeats = 1;
    int instanceDivisor = 1;
    double innerAngle = 0.0;
    double outerAngle = 0.0;
    bool setAngle = false;
};

struct Report {
    std::string method;
    std::size_t nodes = 0;
    std::size_t elements = 0;
    std::size_t storedBytes = 0;
    std::size_t expandedBytes = 0;
    double meshSeconds = 0.0;
    double viewSeconds = 0.0;
    double solveSeconds = 0.0;
    double sumA2 = 0.0;
    std::size_t bandwidth = 0;
    long rssKb = 0;
};

double elapsed(Clock::time_point start)
{
    return std::chrono::duration<double>(Clock::now() - start).count();
}

long peakRssKb()
{
    rusage usage{};
    getrusage(RUSAGE_SELF, &usage);
    return usage.ru_maxrss; // KiB on Linux
}

std::size_t solverMeshBytes(const femm::mesh::SolverMesh &mesh)
{
    std::size_t bytes = sizeof(femm::mesh::SolverMesh);
    bytes += mesh.nodes.size() * sizeof(femm::mesh::SolverMesh::Node);
    bytes += mesh.elements.size() * sizeof(femm::mesh::SolverMesh::Element);
    bytes += mesh.edges.size() * sizeof(femm::mesh::SolverMesh::Edge);
    bytes += mesh.periodicConstraints.size() * sizeof(femm::mesh::SolverMesh::PeriodicConstraint);
    for (const auto &gap : mesh.airGaps) {
        bytes += sizeof(femm::mesh::SolverMesh::AirGap);
        bytes += gap.quadraturePoints.size() * sizeof(femm::mesh::SolverMesh::AirGapQuadraturePoint);
        bytes += gap.innerRing.size() * sizeof(femm::mesh::SolverMesh::AirGapRingPoint);
        bytes += gap.outerRing.size() * sizeof(femm::mesh::SolverMesh::AirGapRingPoint);
        bytes += gap.nodeIndices.size() * sizeof(femm::mesh::MeshIndex);
    }
    return bytes;
}

double fieldChecksum(const TrialSolution &trial)
{
    if (!trial.real)
        return 0.0;
    double sum = 0.0;
    for (double value : trial.real->nodal.magneticVectorPotential)
        sum += value * value;
    return sum;
}

double median(std::vector<double> values)
{
    if (values.empty())
        return 0.0;
    std::sort(values.begin(), values.end());
    return values[values.size() / 2];
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

/** Name of the first periodic/antiperiodic AGE boundary, if any. */
std::string airGapBoundaryName(const femm::tiled::TiledModel &model)
{
    for (const auto &property : model.boundaryProps) {
        const auto *boundary = dynamic_cast<const femm::CMBoundaryProp *>(property.get());
        if (boundary && (boundary->BdryFormat == 6 || boundary->BdryFormat == 7))
            return boundary->BdryName;
    }
    return {};
}

double runRepeats(AnalysisSession &session, const Options &options, double &checksum)
{
    std::vector<double> times;
    for (int repeat = 0; repeat < options.repeats; ++repeat) {
        const auto start = Clock::now();
        const TrialSolution trial = session.solve();
        times.push_back(elapsed(start));
        checksum = fieldChecksum(trial);
    }
    return median(times);
}

Report runTiled(const Options &options, const std::string &method)
{
    std::ifstream input(options.tiled);
    if (!input)
        throw std::runtime_error("could not open tiled model: " + options.tiled);
    std::ostringstream buffer;
    buffer << input.rdbuf();

    femm::tiled::TiledModel model;
    std::vector<femm::tiled::TiledDiagnostic> errors;
    if (!femm::tiled::loadTiledModelJson(buffer.str(), model, errors))
        throw std::runtime_error("could not parse tiled model: " + options.tiled);

    // Optionally shrink the machine to a sector by dividing every tile's repeat
    // count. The closure becomes periodic so the sector is a valid open model.
    if (options.instanceDivisor > 1) {
        for (auto &tile : model.tiles) {
            if (tile.repeat.count % options.instanceDivisor != 0)
                throw std::runtime_error("instance divisor does not divide tile count");
            tile.repeat.count /= options.instanceDivisor;
            tile.repeat.closure = femm::tiled::Closure::Periodic;
        }
        for (auto &override : model.overrides) {
            std::size_t count = 0;
            for (const auto &tile : model.tiles)
                if (tile.name == override.tile)
                    count = tile.repeat.count;
            if (count == 0)
                continue;
            if (override.circuit.size() > count)
                override.circuit.resize(count);
            if (override.turnScale.size() > count)
                override.turnScale.resize(count);
            if (override.magDir.size() > count)
                override.magDir.resize(count);
        }
    }

    const femm::tiled::TiledValidationResult validation = femm::tiled::validateTiledModel(model);
    if (!validation.succeeded())
        throw std::runtime_error("invalid tiled model: " +
                                 validation.diagnostics.front().message);

    Report report;
    report.method = method;

    fmesher::TangleMesherBackend backend;
    const auto meshStart = Clock::now();
    fmesher::TiledMeshResult meshed = fmesher::meshTiledModel(model, backend);
    report.meshSeconds = elapsed(meshStart);
    if (!meshed.ok)
        throw std::runtime_error("meshTiledModel failed: " +
                                 (meshed.diagnostics.empty()
                                      ? options.tiled
                                      : meshed.diagnostics.front().message));

    std::unique_ptr<femm::FemmProblem> owned = femm::tiled::buildTileProblem(model, 0);
    owned->nodelist.clear();
    owned->linelist.clear();
    owned->arclist.clear();
    owned->labellist.clear();
    owned->pathName.clear();

    auto solver = std::make_shared<FSolverAnalysisBackend>();
    AnalysisSession session(ModelDefinition(std::move(owned)), solver);
    session.setInstancedMesh(std::move(meshed.instanced));
    if (method == "native")
        session.setNativeInstanced(true);

    const std::string gapName = airGapBoundaryName(model);
    if (options.setAngle && !gapName.empty())
        session.setAirGapAngle(session.model().airGap(gapName), options.innerAngle,
                               options.outerAngle);

    if (method == "native") {
        // Build and cache the compact view, reporting its construction cost.
        const auto viewStart = Clock::now();
        session.ensureMesh();
        report.viewSeconds = elapsed(viewStart);
    }

    report.solveSeconds = runRepeats(session, options, report.sumA2);

    if (method == "native") {
        const auto view = session.logicalView();
        report.nodes = view ? view->nodeCount() : 0;
        report.elements = view ? view->elementCount() : 0;
        report.storedBytes = view ? view->storedByteCount() : 0;
        report.expandedBytes = view ? view->expandedByteCount() : 0;
        report.bandwidth = static_cast<std::size_t>(solver->nativeBandwidth());
    } else {
        const auto mesh = session.mesh();
        report.nodes = mesh ? mesh->nodes.size() : 0;
        report.elements = mesh ? mesh->elements.size() : 0;
        report.storedBytes = mesh ? solverMeshBytes(*mesh) : 0;
        report.expandedBytes = report.storedBytes;
        report.bandwidth = static_cast<std::size_t>(solver->solvedSolver().BandWidth);
    }
    report.rssKb = peakRssKb();
    return report;
}

Report runConventional(const Options &options)
{
    Report report;
    report.method = "conventional";

    auto solver = std::make_shared<FSolverAnalysisBackend>();
    AnalysisSession session(ModelDefinition(readProblem(options.redraw)),
                            std::make_shared<fmesher::TangleMesherBackend>(), solver);
    const auto meshStart = Clock::now();
    const auto mesh = session.ensureMesh();
    report.meshSeconds = elapsed(meshStart);
    report.nodes = mesh->nodes.size();
    report.elements = mesh->elements.size();
    report.storedBytes = solverMeshBytes(*mesh);
    report.expandedBytes = report.storedBytes;

    report.solveSeconds = runRepeats(session, options, report.sumA2);
    report.bandwidth = static_cast<std::size_t>(solver->solvedSolver().BandWidth);
    report.rssKb = peakRssKb();
    return report;
}

void printHeader()
{
    std::printf("%-14s %9s %9s %9s %9s %8s %8s %8s %8s %12s %8s\n", "method", "nodes",
                "elements", "store_MB", "expand_MB", "mesh_s", "view_s", "solve_s", "band", "sumA2",
                "rss_MB");
}

void printReport(const Report &report)
{
    std::printf("%-14s %9zu %9zu %9.2f %9.2f %8.2f %8.2f %8.2f %8zu %12.4g %8.1f\n",
                report.method.c_str(), report.nodes, report.elements,
                report.storedBytes / (1024.0 * 1024.0),
                report.expandedBytes / (1024.0 * 1024.0), report.meshSeconds,
                report.viewSeconds, report.solveSeconds, report.bandwidth, report.sumA2,
                report.rssKb / 1024.0);
}

void printUsage(const char *program)
{
    std::cout << "usage: " << program
              << " [--method conventional|materialized|native|all]"
                 " [--tiled FILE] [--redraw FILE] [--repeats N]"
                 " [--instance-divisor N] [--position INNER OUTER]\n";
}

} // namespace

int main(int argc, char **argv)
{
    Options options;
    options.tiled = "mfemm/testing/radial_machine/data/radial_machine_tiled.json";
    options.redraw = "mfemm/testing/radial_machine/data/radial_machine_redraw_01.fem";

    for (int i = 1; i < argc; ++i) {
        const std::string argument = argv[i];
        auto next = [&](const char *name) {
            if (i + 1 >= argc)
                throw std::runtime_error(std::string("missing value for ") + name);
            return std::string(argv[++i]);
        };
        if (argument == "--method")
            options.method = next("--method");
        else if (argument == "--tiled")
            options.tiled = next("--tiled");
        else if (argument == "--redraw")
            options.redraw = next("--redraw");
        else if (argument == "--repeats")
            options.repeats = std::stoi(next("--repeats"));
        else if (argument == "--instance-divisor")
            options.instanceDivisor = std::stoi(next("--instance-divisor"));
        else if (argument == "--position") {
            options.innerAngle = std::stod(next("--position"));
            options.outerAngle = std::stod(next("--position"));
            options.setAngle = true;
        } else if (argument == "--help" || argument == "-h") {
            printUsage(argv[0]);
            return 0;
        } else {
            std::cerr << "unknown argument: " << argument << '\n';
            printUsage(argv[0]);
            return 2;
        }
    }
    if (options.repeats < 1)
        options.repeats = 1;

    try {
        printHeader();
        if (options.method == "conventional") {
            printReport(runConventional(options));
        } else if (options.method == "materialized" || options.method == "native") {
            printReport(runTiled(options, options.method));
        } else if (options.method == "all") {
            printReport(runTiled(options, "materialized"));
            printReport(runTiled(options, "native"));
        } else {
            std::cerr << "unknown method: " << options.method << '\n';
            return 2;
        }
    } catch (const std::exception &error) {
        std::cerr << "benchmark failed: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
