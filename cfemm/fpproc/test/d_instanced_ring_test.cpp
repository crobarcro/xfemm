#include "AnalysisSession.h"
#include "FSolverAnalysisBackend.h"
#include "FemmReader.h"
#include "TangleMesherBackend.h"
#include "fpproc.h"
#include "mesh/SolverMeshValidator.h"

#include <tangle_mesh.h>

#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using femm::mesh::SolverMesh;

struct SamplePoint { double x; double y; };
struct Tolerance { double absolute; double relative; };

std::unique_ptr<femm::FemmProblem> readProblem(const std::string &path)
{
    auto problem = std::make_unique<femm::FemmProblem>(femm::FileType::MagneticsFile);
    auto view = std::shared_ptr<femm::FemmProblem>(problem.get(), [](femm::FemmProblem *) {});
    femm::MagneticsReader reader(view, std::cerr);
    if (reader.parse(path) != femm::F_FILE_OK)
        throw std::runtime_error("could not parse " + path);
    return problem;
}

/** Mesh the loaded problem as one tile and repeat it rotationally. */
class TemplateMesher final : public fmesher::MesherBackend {
public:
    TemplateMesher(femm::mesh::TemplateRequest request, std::shared_ptr<int> engineCalls)
        : request_(std::move(request)), engineCalls_(std::move(engineCalls))
    {
    }

    const char *name() const override { return "TangleTemplate"; }

    femm::mesh::MeshResult mesh(femm::FemmProblem &problem,
                                const femm::mesh::MeshingRequest &request) override
    {
        femm::mesh::MeshingRequest templated = request;
        templated.templates = {request_};
        fmesher::TangleMesherBackend tangle(
            [this](const std::string &path, const ::MeshOptions &options, ::Mesh &mesh) {
                ++(*engineCalls_);
                return tangle_mesh_fem(path, options, mesh);
            });
        return tangle.mesh(problem, templated);
    }

private:
    femm::mesh::TemplateRequest request_;
    std::shared_ptr<int> engineCalls_;
};

struct Run {
    std::unique_ptr<femm::AnalysisSession> session;
    std::shared_ptr<femm::FSolverAnalysisBackend> solver;
    std::unique_ptr<FPProc> postprocessor;
    std::shared_ptr<const SolverMesh> mesh;
};

Run runInstanced(const std::string &path, std::size_t instanceCount,
                 std::shared_ptr<int> engineCalls)
{
    femm::mesh::TemplateRequest request;
    request.centerXMetres = 0.0;
    request.centerYMetres = 0.0;
    request.instanceCount = instanceCount;
    request.totalAngleDegrees = 360.0;
    auto solver = std::make_shared<femm::FSolverAnalysisBackend>();
    auto session = std::make_unique<femm::AnalysisSession>(
        femm::ModelDefinition(readProblem(path)),
        std::make_shared<TemplateMesher>(request, engineCalls), solver);
    const auto mesh = session->ensureMesh();
    if (!femm::mesh::validateSolverMesh(*mesh).valid())
        throw std::runtime_error("instanced mesh is invalid");
    if (mesh->nodes.empty() || mesh->elements.empty() || mesh->edges.empty())
        throw std::runtime_error("instanced mesh is empty");
    session->solve();
    auto postprocessor = std::make_unique<FPProc>();
    if (!postprocessor->OpenDocument(session->model().problem(), solver->solvedSolver(),
                                     solver->solvedSystem()))
        throw std::runtime_error("instanced post-processing failed");
    for (auto &label : postprocessor->blocklist)
        label.IsSelected = true;
    return {std::move(session), solver, std::move(postprocessor), mesh};
}

Run runControl(const std::string &path)
{
    auto solver = std::make_shared<femm::FSolverAnalysisBackend>();
    auto session = std::make_unique<femm::AnalysisSession>(
        femm::ModelDefinition(readProblem(path)),
        std::make_shared<fmesher::TangleMesherBackend>(), solver);
    const auto mesh = session->ensureMesh();
    if (!femm::mesh::validateSolverMesh(*mesh).valid())
        throw std::runtime_error("control mesh is invalid");
    session->solve();
    auto postprocessor = std::make_unique<FPProc>();
    if (!postprocessor->OpenDocument(session->model().problem(), solver->solvedSolver(),
                                     solver->solvedSystem()))
        throw std::runtime_error("control post-processing failed");
    for (auto &label : postprocessor->blocklist)
        label.IsSelected = true;
    return {std::move(session), solver, std::move(postprocessor), mesh};
}

void compare(const char *quantity, double instanced, double control, Tolerance tolerance)
{
    const double difference = std::abs(instanced - control);
    const double permitted =
        tolerance.absolute + tolerance.relative * std::max(std::abs(instanced), std::abs(control));
    std::cout << quantity << ": instanced=" << instanced << " control=" << control
              << " abs-difference=" << difference << " tolerance=" << permitted << '\n';
    if (!std::isfinite(instanced) || !std::isfinite(control) || !std::isfinite(difference) ||
        !std::isfinite(permitted) || difference > permitted)
        throw std::runtime_error(std::string(quantity) + " differential failed");
}

} // namespace

int main(int argc, char **argv)
{
    if (argc != 4) {
        std::cerr << "expected tile, control, and instance count arguments\n";
        return 2;
    }
    try {
        const std::size_t instanceCount = static_cast<std::size_t>(std::stoul(argv[3]));
        auto engineCalls = std::make_shared<int>(0);
        auto instanced = runInstanced(argv[1], instanceCount, engineCalls);
        if (*engineCalls != 1)
            throw std::runtime_error("Tangle was not called exactly once for the tile");

        // The expanded mesh must contain one instance's worth of provenance.
        if (instanced.mesh->nodes.size() >= 1000000)
            throw std::runtime_error("instanced mesh is implausibly large");

        auto control = runControl(argv[2]);

        const std::vector<SamplePoint> samples{
            {0.075, 0.0}, {0.0, 0.075}, {-0.075, 0.0}, {0.0, -0.075},
            {0.069290964938346522, 0.028701257427381738}};
        for (const auto &point : samples) {
            CMPointVals iv, cv;
            if (!instanced.postprocessor->GetPointValues(point.x, point.y, iv) ||
                !control.postprocessor->GetPointValues(point.x, point.y, cv))
                throw std::runtime_error("sample point is outside a generated mesh");
            compare("A", iv.A.re, cv.A.re, {1e-3, 0.01});
            compare("B1", iv.B1.re, cv.B1.re, {1e-3, 0.03});
            compare("B2", iv.B2.re, cv.B2.re, {1e-3, 0.03});
        }

        const double instancedEnergy = instanced.postprocessor->BlockIntegral(2).re;
        const double controlEnergy = control.postprocessor->BlockIntegral(2).re;
        if (!std::isfinite(instancedEnergy) || !std::isfinite(controlEnergy) ||
            instancedEnergy <= 1e-12 || controlEnergy <= 1e-12)
            throw std::runtime_error("trivial or non-finite energy");
        compare("energy", instancedEnergy, controlEnergy, {1e-4, 0.02});

        std::cout << "instanced nodes=" << instanced.mesh->nodes.size()
                  << " elements=" << instanced.mesh->elements.size()
                  << " control nodes=" << control.mesh->nodes.size()
                  << " elements=" << control.mesh->elements.size() << '\n';
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
