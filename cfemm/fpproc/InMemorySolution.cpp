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
    blockproplist.reserve(solver.blockproplist.size());
    for (const auto &source : solver.blockproplist) {
        blockproplist.push_back(static_cast<const femm::CMMaterialProp &>(source));
        auto &item = blockproplist.back();
        // MuMax also acts as the DC incremental-permeability flag, so clear it
        // for ordinary problems after copying the complete material state.
        item.MuMax = problem.PrevType ? source.MuMax : 0.;
        item.clearSlopes();
        item.GetSlopes(2. * std::acos(-1.) * Frequency);
        if (item.BHpoints > 0 && item.slope.size() < static_cast<std::size_t>(item.BHpoints))
            throw std::runtime_error("could not prepare nonlinear material slopes");
    }

    // Keep the post-processor's circuit list at the user-visible model level.
    // FSolver expands series circuits into one internal circuit per winding
    // block, with Amps multiplied by signed turns, and rewrites each label's
    // InCircuit to point at that internal entry.  Those expanded circuits are
    // an assembly detail and must not leak into getcircuitprops()/linkage
    // calculations.  Copy the original circuit identity/type, but retain the
    // solve-time state (notably a current changed through AnalysisSession).
    if (solver.circproplist.size() < problem.circproplist.size())
        throw std::runtime_error("solver circuit list is shorter than problem circuit list");
    circproplist.clear();
    circproplist.reserve(problem.circproplist.size());
    for (std::size_t i = 0; i < problem.circproplist.size(); ++i) {
        const auto *source = dynamic_cast<const femm::CMCircuit *>(problem.circproplist[i].get());
        if (!source)
            throw std::invalid_argument("magnetic problem contains a non-magnetic circuit property");
        femm::CMCircuit circuit = *source;
        const auto &solved = solver.circproplist[i];
        circuit.Amps = solved.Amps;
        circuit.dVolts = solved.dVolts;
        circuit.J = solved.J;
        circuit.dV = solved.dV;
        circuit.Case = solved.Case;
        circuit.OrigCirc = -1;
        circproplist.push_back(circuit);
    }

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
        const int solverCircuit = blocklist[i].InCircuit;
        if (solverCircuit < 0) {
            blocklist[i].Case = 1;
            blocklist[i].J = 0.;
            continue;
        }
        if (static_cast<std::size_t>(solverCircuit) >= solver.circproplist.size())
            throw std::runtime_error("block label refers to an invalid solver circuit");

        const auto &solvedCircuit = solver.circproplist[solverCircuit];
        if (solvedCircuit.Case == 0) {
            blocklist[i].Case = 0;
            blocklist[i].dVolts = solvedCircuit.dV;
        } else {
            blocklist[i].Case = 1;
            blocklist[i].J = solvedCircuit.J;
        }

        // Preserve the per-block solved excitation copied above, but restore
        // the public circuit index expected by FPProc's series-circuit voltage
        // and flux-linkage routines.  Expanded circuits remember their source
        // in OrigCirc; ordinary circuits keep their original index directly.
        const int originalCircuit = solvedCircuit.OrigCirc >= 0
                                  ? solvedCircuit.OrigCirc : solverCircuit;
        if (originalCircuit < 0 ||
            static_cast<std::size_t>(originalCircuit) >= circproplist.size())
            throw std::runtime_error("solver circuit does not map to a problem circuit");
        blocklist[i].InCircuit = originalCircuit;
    }

    return finalizeSolution();
}
