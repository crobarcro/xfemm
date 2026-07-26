#ifndef FEMM_MESH_SOLVERMESH_H
#define FEMM_MESH_SOLVERMESH_H

#include "RawMesh.h"

#include <array>
#include <string>
#include <vector>

namespace femm {
namespace mesh {

/** Marker/property index after decoding. Property indices are zero-based. */
using PropertyIndex = std::size_t;

/** Invalid value for every PropertyIndex field. It denotes "no property". */
constexpr PropertyIndex InvalidPropertyIndex = std::numeric_limits<PropertyIndex>::max();

/**
 * Mesh ready to be consumed by a FEMM solver.
 *
 * Coordinates, radii, and centres are expressed in metres, regardless of the
 * length unit used in the source problem. Angles are expressed in degrees.
 */
struct SolverMesh {
    struct Node {
        double x = 0.0; ///< x coordinate in metres.
        double y = 0.0; ///< y coordinate in metres.
        /** Zero-based point-property index, or InvalidPropertyIndex. */
        PropertyIndex pointProperty = InvalidPropertyIndex;
    };

    struct Element {
        /** Zero-based indices into SolverMesh::nodes; InvalidMeshIndex means unset. */
        std::array<MeshIndex, 3> nodes{{InvalidMeshIndex, InvalidMeshIndex, InvalidMeshIndex}};
        /** Validated zero-based index into the problem's block-label list. */
        PropertyIndex blockLabel = InvalidPropertyIndex;
    };

    struct BoundaryEdge {
        /** Zero-based indices into SolverMesh::nodes; InvalidMeshIndex means unset. */
        MeshIndex first = InvalidMeshIndex;
        MeshIndex second = InvalidMeshIndex;
        /** Zero-based boundary-property index, or InvalidPropertyIndex. */
        PropertyIndex boundaryProperty = InvalidPropertyIndex;
    };

    enum class Periodicity { Periodic, Antiperiodic };

    struct PeriodicNodePair {
        /** Zero-based indices into SolverMesh::nodes; InvalidMeshIndex means unset. */
        MeshIndex first = InvalidMeshIndex;
        MeshIndex second = InvalidMeshIndex;
        Periodicity periodicity = Periodicity::Periodic;
    };

    /** Interpolation data corresponding to one femm::CQuadPoint. */
    struct AirGapQuadraturePoint {
        /** Zero-based solver-node indices; InvalidMeshIndex means unset. */
        std::array<MeshIndex, 4> nodes{{InvalidMeshIndex, InvalidMeshIndex,
                                      InvalidMeshIndex, InvalidMeshIndex}};
        std::array<double, 4> weights{{0.0, 0.0, 0.0, 0.0}};
    };

    /** Value-only data sufficient to construct and populate CAirGapElement. */
    struct AirGap {
        std::string boundaryName;
        Periodicity periodicity = Periodicity::Periodic;
        std::size_t totalArcElements = 0;
        double totalArcLengthDegrees = 0.0;
        double innerRadius = 0.0; ///< Metres.
        double outerRadius = 0.0; ///< Metres.
        double innerAngleDegrees = 0.0;
        double outerAngleDegrees = 0.0;
        double innerShift = 0.0; ///< Fraction of an annular mesh element.
        double outerShift = 0.0; ///< Fraction of an annular mesh element.
        double centerX = 0.0; ///< Metres.
        double centerY = 0.0; ///< Metres.
        std::vector<AirGapQuadraturePoint> quadraturePoints;
        /** Zero-based indices into SolverMesh::nodes; InvalidMeshIndex means unset. */
        std::vector<MeshIndex> nodeIndices;
    };

    std::vector<Node> nodes;
    std::vector<Element> elements;
    std::vector<BoundaryEdge> boundaryEdges;
    std::vector<PeriodicNodePair> periodicNodePairs;
    std::vector<AirGap> airGaps;
};

} // namespace mesh
} // namespace femm

#endif
