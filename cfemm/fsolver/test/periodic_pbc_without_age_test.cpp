#include "fsolver.h"

#include <iostream>

int main(int argc, char **argv)
{
    if (argc != 2) {
        std::cerr << "expected the legacy mesh fixture path without an extension\n";
        return 2;
    }

    FSolver solver;
    solver.PathName = argv[1];
    if (!solver.LoadProblemFile()) {
        std::cerr << "could not load the periodic FEM problem\n";
        return 1;
    }
    if (solver.LoadMesh(false) != NOERROR) {
        std::cerr << "could not load the legacy periodic mesh and .pbc file\n";
        return 1;
    }
    if (solver.NumPBCs <= 0 || solver.pbclist.empty()) {
        std::cerr << "legacy .pbc ordinary periodic pairs were not loaded\n";
        return 1;
    }
    if (solver.NumAirGapElems != 0 || !solver.agelist.empty()) {
        std::cerr << "legacy periodic .pbc with zero AGE records created an AGE\n";
        return 1;
    }
    return 0;
}
