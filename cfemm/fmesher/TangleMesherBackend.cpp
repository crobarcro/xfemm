#include "TangleMesherBackend.h"

#include "TriangleMesherBackend.h"

namespace fmesher {

femm::mesh::MeshResult TangleMesherBackend::mesh(
        femm::FemmProblem &problem, bool periodic,
        const femm::mesh::MeshingOptions &options)
{
    // Both engines share FEMM's discretisation and periodic/AGE construction.
    // Run that pipeline in memory and return its value representation rather
    // than round-tripping Triangle files.  This is also the compatibility
    // implementation used when Tangle is built without an external kernel.
    TriangleMesherBackend compatibilityKernel;
    auto result = compatibilityKernel.mesh(problem, periodic, options);
    for (auto &diagnostic : result.diagnostics)
        diagnostic.backendName = "Tangle";
    return result;
}

} // namespace fmesher
