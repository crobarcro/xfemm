#ifndef FEMM_MESH_SOLVERMESH_H
#define FEMM_MESH_SOLVERMESH_H

#include "RawMesh.h"

#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace femm {
namespace mesh {

/**
 * Value-based mesh transferred from a mesher to a FEMM solver.
 *
 * Coordinates, radii, and centres are expressed in metres, regardless of the
 * length unit used in the source problem. Angles are expressed in degrees. The
 * collections own all of their data and contain no solver solution state.
 *
 * Invariants:
 * - MeshIndex values are zero-based indices into `nodes`.
 * - Every Element has exactly three nodes and all node references are valid.
 * - Element region attributes are the raw Triangle attributes corresponding to
 *   solver block labels (Triangle writes block-label index + 1).
 * - Node and edge boundary markers retain Triangle's raw signed values. Marker
 *   decoding therefore remains an explicit operation at the consumer, e.g.
 *   `marker > 1 ? marker - 2 : -1` for a node point-property marker.
 */
struct SolverMesh {
    struct Node {
        double x = 0.0; ///< x coordinate in metres.
        double y = 0.0; ///< y coordinate in metres.
        /** Raw Triangle boundary marker; deliberately not a property index. */
        std::int32_t boundaryMarker = 0;
    };

    struct Element {
        /** Zero-based indices into SolverMesh::nodes; InvalidMeshIndex means unset. */
        std::array<MeshIndex, 3> nodes{{InvalidMeshIndex, InvalidMeshIndex, InvalidMeshIndex}};
        /** Raw Triangle region attribute (normally one-based block-label id). */
        std::int32_t regionAttribute = 0;
    };

    struct Edge {
        /** Zero-based indices into SolverMesh::nodes; InvalidMeshIndex means unset. */
        MeshIndex first = InvalidMeshIndex;
        MeshIndex second = InvalidMeshIndex;
        /** Raw Triangle boundary marker; deliberately not a property index. */
        std::int32_t boundaryMarker = 0;
    };

    enum class Periodicity { Periodic, Antiperiodic };

    /**
     * One node equivalence written in the primary section of a `.pbc` file.
     * These constraints are used by every periodic or antiperiodic model; they
     * are not specific to air-gap elements.
     */
    struct PeriodicConstraint {
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

    /** One virtual node in a reusable, full 360-degree AGE ring. */
    struct AirGapRingPoint {
        MeshIndex node = InvalidMeshIndex;
        /** Angular position in units of one annular element, before positioning. */
        double elementPosition = 0.0;
        /** Periodic sign (+1, or alternating +/-1 for an antiperiodic sector). */
        double weight = 1.0;
    };

    /**
     * Optional value-only air-gap data represented by the magnetic air-gap
     * extension following the periodic constraints in a `.pbc` file. A `.pbc`
     * file is not, in general, an air-gap file: periodic models without an air
     * gap have periodicConstraints and an empty airGaps collection.
     *
     * This is mesh topology/interpolation data, not a solver field quantity.
     */
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
        /** Position-independent ring topology used to regenerate quadraturePoints. */
        std::vector<AirGapRingPoint> innerRing;
        std::vector<AirGapRingPoint> outerRing;
        /** Zero-based indices into SolverMesh::nodes; InvalidMeshIndex means unset. */
        std::vector<MeshIndex> nodeIndices;
    };

    std::vector<Node> nodes;
    std::vector<Element> elements;
    std::vector<Edge> edges;
    /** Primary `.pbc` node-pair section, applicable to all periodic models. */
    std::vector<PeriodicConstraint> periodicConstraints;
    /** Optional magnetic air-gap extension of `.pbc`. */
    std::vector<AirGap> airGaps;
};

} // namespace mesh
} // namespace femm

#endif
