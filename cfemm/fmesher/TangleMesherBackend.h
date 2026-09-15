#ifndef FEMM_FMESHER_TANGLEMESHERBACKEND_H
#define FEMM_FMESHER_TANGLEMESHERBACKEND_H

#include "MesherBackend.h"

#include <functional>
#include <string>

struct Mesh;
struct MeshOptions;

namespace fmesher {

/**
 * In-memory Tangle mesher adapter.
 *
 * Tangle deliberately exposes the same value-only SolverMesh boundary as the
 * Triangle adapter. Its library entry point currently consumes the source
 * FEMM path and returns the generated mesh in memory; no Triangle code or
 * intermediate mesh files are used by this adapter. The engine functor takes
 * Tangle's own MeshOptions so callers can exercise the option mapping without
 * changing the backend interface.
 */
class TangleMesherBackend final : public MesherBackend {
public:
    using Engine = std::function<int(const std::string &, const ::MeshOptions &, ::Mesh &)>;

    explicit TangleMesherBackend(Engine engine = {});

    using MesherBackend::mesh;
    const char *name() const override { return "Tangle"; }
    femm::mesh::MeshResult mesh(femm::FemmProblem &,
                                const femm::mesh::MeshingRequest &) override;

private:
    Engine engine_;
};

} // namespace fmesher

#endif
