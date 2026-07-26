#ifndef FEMM_FMESHER_MESHERBACKEND_H
#define FEMM_FMESHER_MESHERBACKEND_H

#include "mesh/Meshing.h"

namespace femm { class FemmProblem; }

namespace fmesher {

/** Abstract mesher backend. Backend implementation details never escape this API. */
class MesherBackend {
public:
    virtual ~MesherBackend() = default;
    virtual femm::mesh::MeshResult mesh(femm::FemmProblem &problem,
                                         bool periodic,
                                         const femm::mesh::MeshingOptions &options = {}) = 0;
};

} // namespace fmesher
#endif
