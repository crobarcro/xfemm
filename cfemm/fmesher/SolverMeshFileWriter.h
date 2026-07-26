#ifndef FEMM_FMESHER_SOLVERMESHFILEWRITER_H
#define FEMM_FMESHER_SOLVERMESHFILEWRITER_H
#include "mesh/SolverMesh.h"
#include <string>
namespace femm { class FemmProblem; }
namespace fmesher {
class SolverMeshFileWriter {
public:
    static bool write(const femm::mesh::SolverMesh &, const femm::FemmProblem &,
                      const std::string &, int (*warn)(const char *, ...));
};
}
#endif
