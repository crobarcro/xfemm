#include "AnalysisSession.h"
#include "FSolverAnalysisBackend.h"
#include "fpproc.h"

#include "CBlockLabel.h"
#include "CCircuit.h"
#include "CMaterialProp.h"
#include "CPointProp.h"
#include "MesherBackend.h"

#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>

namespace {

using namespace femm;

/**
 * A zero-current "series" circuit is expanded by the session into one
 * sub-circuit per block label. The in-memory post-processor must still expose
 * the original circuit so its flux linkage can be formed from the coil labels,
 * exactly as the legacy .ans path does. Before the fix this returned NaN.
 */
class SquareMesher final : public fmesher::MesherBackend {
public:
    const char *name() const override { return "Square"; }
    femm::mesh::MeshResult mesh(femm::FemmProblem &,
                                const femm::mesh::MeshingRequest &) override
    {
        femm::mesh::MeshResult result;
        result.status = femm::mesh::MeshStatus::Success;
        result.mesh.nodes = {{0, 0, 2}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}};
        result.mesh.elements.push_back({{{0, 1, 2}}, 1});
        result.mesh.elements.push_back({{{0, 2, 3}}, 1});
        result.mesh.edges = {{0, 1, 0}, {1, 2, 0}, {2, 3, 0}, {3, 0, 0}, {0, 2, 0}};
        return result;
    }
};

std::unique_ptr<FemmProblem> makeProblem()
{
    auto problem = std::make_unique<FemmProblem>(FileType::MagneticsFile);

    auto point = std::make_unique<CMPointProp>();
    point->PointName = "fixed";
    point->A = 0;
    problem->nodeproplist.push_back(std::move(point));

    auto material = std::make_unique<CMMaterialProp>();
    material->BlockName = "copper";
    material->mu_x = material->mu_y = 1;
    material->H_c = 1;
    problem->blockproplist.push_back(std::move(material));

    auto circuit = std::make_unique<CMCircuit>();
    circuit->CircName = "phase";
    circuit->CircType = 1; // series
    circuit->Amps = 0;     // zero current: flux linkage is from the magnet
    problem->circproplist.push_back(std::move(circuit));

    auto label = std::make_unique<CMBlockLabel>();
    label->BlockType = 0;
    label->BlockTypeName = "copper";
    label->InCircuit = 0;
    label->Turns = 1;
    problem->labellist.push_back(std::move(label));
    return problem;
}

} // namespace

int main()
{
    try {
        auto solver = std::make_shared<FSolverAnalysisBackend>();
        AnalysisSession session(ModelDefinition(makeProblem()),
                                std::make_shared<SquareMesher>(), solver);
        session.solve();

        FPProc postprocessor;
        if (!postprocessor.OpenDocument(session.model().problem(),
                                        solver->solvedSolver(), solver->solvedSystem()))
            throw std::runtime_error("in-memory post-processing failed");

        const CComplex flux = postprocessor.GetFluxLinkage(0);
        if (!std::isfinite(flux.re) || !std::isfinite(flux.im)) {
            std::cerr << "non-finite series-circuit flux linkage: " << flux.re
                      << ", " << flux.im << '\n';
            return 1;
        }
        std::cout << "in-memory series-circuit flux linkage = " << flux.re << '\n';
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
