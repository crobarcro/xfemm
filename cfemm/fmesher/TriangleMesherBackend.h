#ifndef FEMM_FMESHER_TRIANGLEMESHERBACKEND_H
#define FEMM_FMESHER_TRIANGLEMESHERBACKEND_H
#include "MesherBackend.h"
#include <string>
namespace fmesher {
class TriangleMesherBackend final : public MesherBackend {
public:
    using MesherBackend::mesh;
    const char *name() const override { return "Triangle"; }
    femm::mesh::MeshResult mesh(femm::FemmProblem &,
                                const femm::mesh::MeshingRequest &) override;
private:
    friend class FMesher;
    int (*WarnMessage)(const char *, ...) = nullptr;
    int (*TriMessage)(const char *, ...) = nullptr;
    bool writePolyFiles = false;
    std::string compatibilityPath;
};
}
#endif
