#ifndef FEMM_FMESHER_TRIANGLEMESHERBACKEND_H
#define FEMM_FMESHER_TRIANGLEMESHERBACKEND_H
#include "MesherBackend.h"
#include <string>
namespace fmesher {
class TriangleMesherBackend final : public MesherBackend {
public:
    femm::mesh::MeshResult mesh(femm::FemmProblem &, bool periodic,
                                 const femm::mesh::MeshingOptions & = {}) override;
private:
    friend class FMesher;
    int (*WarnMessage)(const char *, ...) = nullptr;
    int (*TriMessage)(const char *, ...) = nullptr;
    bool writePolyFiles = false;
    std::string compatibilityPath;
};
}
#endif
