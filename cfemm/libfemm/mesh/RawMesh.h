#ifndef FEMM_MESH_RAWMESH_H
#define FEMM_MESH_RAWMESH_H

#include <array>
#include <cstddef>
#include <limits>
#include <vector>

namespace femm {
namespace mesh {

/** Index into a mesh-owned vector. All mesh indices are zero-based. */
using MeshIndex = std::size_t;

/** Invalid value for every MeshIndex field in the mesh interchange types. */
constexpr MeshIndex InvalidMeshIndex = std::numeric_limits<MeshIndex>::max();

/** Backend output before FEMM markers and region attributes are decoded. */
struct RawMesh {
    struct Point {
        double x = 0.0;
        double y = 0.0;
        int marker = 0; ///< Backend marker; its meaning is deliberately not interpreted here.
    };

    struct Triangle {
        /** Zero-based indices into RawMesh::points; InvalidMeshIndex means unset. */
        std::array<MeshIndex, 3> nodes{{InvalidMeshIndex, InvalidMeshIndex, InvalidMeshIndex}};
        double regionAttribute = 0.0; ///< Uninterpreted backend region attribute.
    };

    struct Edge {
        /** Zero-based indices into RawMesh::points; InvalidMeshIndex means unset. */
        MeshIndex first = InvalidMeshIndex;
        MeshIndex second = InvalidMeshIndex;
        int marker = 0; ///< Uninterpreted backend edge marker.
    };

    std::vector<Point> points;
    std::vector<Triangle> triangles;
    std::vector<Edge> edges;
};

} // namespace mesh
} // namespace femm

#endif
