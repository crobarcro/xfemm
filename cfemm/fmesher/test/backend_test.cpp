#include "MesherBackend.h"
#include "TriangleMesherBackend.h"
#include "fmesher.h"
#include "FemmReader.h"
#include <iostream>
#include <cstdio>
#include <memory>

int main(int argc, char **argv)
{
    if (argc != 2) return 2;
    fmesher::FMesher facade;
    facade.problem->filetype = femm::FileType::MagneticsFile;
    femm::MagneticsReader reader(facade.problem, std::cerr);
    if (reader.parse(argv[1]) != femm::F_FILE_OK) return 3;
    std::unique_ptr<fmesher::MesherBackend> backend(new fmesher::TriangleMesherBackend);
    femm::mesh::MeshingOptions options;
    const auto result = backend->mesh(*facade.problem, false, options);
    std::remove("xfemm-backend.node"); std::remove("xfemm-backend.edge");
    std::remove("xfemm-backend.ele"); std::remove("xfemm-backend.pbc");
    return result.succeeded() && !result.mesh.nodes.empty() && !result.mesh.elements.empty() ? 0 : 1;
}
