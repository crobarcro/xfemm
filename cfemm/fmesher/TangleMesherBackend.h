#ifndef FEMM_FMESHER_TANGLEMESHERBACKEND_H
#define FEMM_FMESHER_TANGLEMESHERBACKEND_H

#include "MesherBackend.h"

#include <functional>
#include <string>

struct Mesh;
struct MeshOptions;
struct FemProblem;

namespace fmesher {

/**
 * In-memory Tangle mesher adapter.
 *
 * Tangle deliberately exposes the same value-only SolverMesh boundary as the
 * Triangle adapter. It can mesh either a source FEMM path or a FEMM-like
 * problem already held in memory; no Triangle code or intermediate mesh files
 * are used by this adapter. The engine functors take Tangle's own MeshOptions
 * so callers can exercise the option mapping without changing the backend
 * interface.
 */
class TangleMesherBackend final : public MesherBackend {
public:
    using Engine = std::function<int(const std::string &, const ::MeshOptions &, ::Mesh &)>;
    /** Meshes a FEMM-like problem converted from xfemm's FemmProblem. */
    using InMemoryEngine =
        std::function<int(const ::FemProblem &, const ::MeshOptions &, ::Mesh &)>;

    explicit TangleMesherBackend(Engine engine = {}, InMemoryEngine inMemoryEngine = {});

    using MesherBackend::mesh;
    const char *name() const override { return "Tangle"; }
    femm::mesh::MeshResult mesh(femm::FemmProblem &,
                                const femm::mesh::MeshingRequest &) override;

private:
    Engine engine_;
    InMemoryEngine inMemoryEngine_;
};

} // namespace fmesher

#endif
