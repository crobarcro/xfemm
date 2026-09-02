#include "fpproc.h"

#include "fsolver.h"
#include "linsolve/LinearSystemBackend.h"

#include <stdexcept>
#include <cmath>

bool FPProc::OpenDocument(const femm::FemmProblem &problem, const FSolver &solver,
                          const femm::LinearSystemBackend<double> &solution)
{
    if (solution.dimension() != solver.NumNodes)
        throw std::invalid_argument("solution and solver node counts differ");

    NewDocument();
    Frequency = solver.Frequency;
    Depth = problem.Depth;
    Precision = problem.Precision;
    LengthUnits = problem.LengthUnits;
    problemType = problem.problemType;
    Coords = problem.Coords;
    ProblemNote = problem.comment;
    extRo = problem.extRo;
    extRi = problem.extRi;
    extZo = problem.extZo;
    PrevSoln = problem.previousSolutionFile;
    PrevType = problem.PrevType;
    bIncremental = problem.PrevType;

    for (const auto &item : problem.nodelist) nodelist.push_back(*item);
    for (const auto &item : problem.linelist) linelist.push_back(*item);
    for (const auto &item : problem.arclist) arclist.push_back(*item);
    blocklist = solver.labellist;
    nodeproplist = solver.nodeproplist;
    lineproplist = solver.lineproplist;
    blockproplist.clear();
    // CMMaterialProp's legacy copy constructor does not preserve the transient
    // MuMax flag. Avoid vector-growth copies after each explicitly prepared item.
    blockproplist.reserve(solver.blockproplist.size());
    for (const auto &source : solver.blockproplist) {
        blockproplist.push_back(static_cast<const femm::CMMaterialProp &>(source));
        auto &item = blockproplist.back();
        item.MuMax = problem.PrevType ? source.MuMax : 0.;
        item.Frequency = source.Frequency;
        item.mu_fdx = source.mu_fdx;
        item.mu_fdy = source.mu_fdy;
        item.clearSlopes();
        item.GetSlopes(2. * std::acos(-1.) * Frequency);
        if (item.BHpoints > 0 && item.slope.size() < static_cast<std::size_t>(item.BHpoints))
            throw std::runtime_error("could not prepare nonlinear material slopes");
    }
    circproplist = solver.circproplist;

    const double centimetresPerSourceUnit[] = {2.54, 0.1, 1., 100., 0.00254, 1.e-04};
    const double coordinateScale = centimetresPerSourceUnit[LengthUnits];
    meshnode.resize(solver.meshnode.size());
    for (std::size_t i = 0; i < solver.meshnode.size(); ++i) {
        meshnode[i].x = solver.meshnode[i].x / coordinateScale;
        meshnode[i].y = solver.meshnode[i].y / coordinateScale;
        meshnode[i].A = solution.rhs()[i];
    }

    meshelem.resize(solver.meshele.size());
    for (std::size_t i = 0; i < solver.meshele.size(); ++i) {
        static_cast<femmsolver::CMElement &>(meshelem[i]) = solver.meshele[i];
        meshelem[i].blk = blocklist[meshelem[i].lbl].BlockType;
    }

    // The solver stores AGE lengths in centimetres; FPProc uses metres here.
    agelist = solver.agelist;
    for (auto &age : agelist) {
        age.ri /= 100.;
        age.ro /= 100.;
        age.agc /= 100.;
    }
    NumAirGapElems = static_cast<int>(agelist.size());
    pmeshnode = &meshnode;
    pmeshelem = &meshelem;

    for (std::size_t i = 0; i < blocklist.size(); ++i) {
        const int circuit = blocklist[i].InCircuit;
        if (circuit < 0) {
            blocklist[i].Case = 1;
            blocklist[i].J = 0.;
        } else if (circproplist[circuit].Case == 0) {
            blocklist[i].Case = 0;
            blocklist[i].dVolts = circproplist[circuit].dVolts;
        } else {
            blocklist[i].Case = 1;
            blocklist[i].J = circproplist[circuit].J;
        }
    }

    return finalizeSolution();
}
