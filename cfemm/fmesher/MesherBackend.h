#ifndef FEMM_FMESHER_MESHERBACKEND_H
#define FEMM_FMESHER_MESHERBACKEND_H

#include "mesh/Meshing.h"

namespace femm { class FemmProblem; }

namespace fmesher {

/** Abstract mesher backend. Backend implementation details never escape this API. */
class MesherBackend {
public:
    virtual ~MesherBackend() = default;

    /** Short backend name for diagnostics, e.g. "Triangle" or "Tangle". */
    virtual const char *name() const = 0;

    /** Mesh using an explicit backend-neutral request. */
    virtual femm::mesh::MeshResult mesh(femm::FemmProblem &problem,
                                         const femm::mesh::MeshingRequest &request) = 0;

    /**
     * Compatibility adapter for the pre-request API. The periodic flag maps to
     * the request-wide "emit field constraints" choice.
     */
    femm::mesh::MeshResult mesh(femm::FemmProblem &problem, bool periodic,
                                const femm::mesh::MeshingOptions &options = {})
    {
        femm::mesh::MeshingRequest request;
        request.options = options;
        request.createPeriodicFieldConstraints = periodic;
        return mesh(problem, request);
    }
};

} // namespace fmesher
#endif
