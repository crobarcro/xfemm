#ifndef FEMM_FMESHER_TANGLEMESHERBACKEND_H
#define FEMM_FMESHER_TANGLEMESHERBACKEND_H

#include "MesherBackend.h"

namespace fmesher {

/**
 * In-memory Tangle mesher adapter.
 *
 * Tangle deliberately exposes the same value-only SolverMesh boundary as the
 * Triangle adapter.  The current Tangle implementation uses the mature FEMM
 * PSLG preparation pipeline; keeping that detail behind this class prevents
 * callers from depending on a particular triangulation engine.
 */
class TangleMesherBackend final : public MesherBackend {
public:
    femm::mesh::MeshResult mesh(femm::FemmProblem &, bool periodic,
                                const femm::mesh::MeshingOptions & = {}) override;
};

} // namespace fmesher

#endif
