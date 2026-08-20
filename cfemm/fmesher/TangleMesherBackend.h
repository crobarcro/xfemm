#ifndef FEMM_FMESHER_TANGLEMESHERBACKEND_H
#define FEMM_FMESHER_TANGLEMESHERBACKEND_H

#include "MesherBackend.h"

#include <functional>
#include <string>

struct Mesh;

namespace fmesher {

/**
 * In-memory Tangle mesher adapter.
 *
 * Tangle deliberately exposes the same value-only SolverMesh boundary as the
 * Triangle adapter.  Its library entry point currently consumes the source
 * FEMM path and returns the generated mesh in memory; no Triangle code or
 * intermediate mesh files are used by this adapter.
 */
class TangleMesherBackend final : public MesherBackend {
public:
    using Engine = std::function<int(const std::string &, ::Mesh &)>;

    explicit TangleMesherBackend(Engine engine = {});

    femm::mesh::MeshResult mesh(femm::FemmProblem &, bool periodic,
                                const femm::mesh::MeshingOptions & = {}) override;

private:
    Engine engine_;
};

} // namespace fmesher

#endif
