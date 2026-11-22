#ifndef FEMM_TILEDMODEL_H
#define FEMM_TILEDMODEL_H

#include "CBoundaryProp.h"
#include "CCircuit.h"
#include "CMaterialProp.h"
#include "CPointProp.h"
#include "femmenums.h"

#include <cstddef>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace femm {

class FemmProblem;

namespace tiled {

/**
 * A tile is a fundamental domain of a rotational tiling: the whole modelled
 * geometry is the tile repeated by a fixed pitch about a rotation centre.
 * This value model is deliberately independent of FemmProblem, the mesher, and
 * the file format so it can be built from JSON, from a .fem problem, or
 * programmatically.
 */

/** One tile vertex. boundaryMarker is zero-based into TiledModel::boundaryProps. */
struct TileNode {
    double x = 0.0;
    double y = 0.0;
    /** Zero-based boundary property index, or -1 for none. */
    int boundaryMarker = -1;
    int group = 0;
};

struct TileSegment {
    int n0 = 0;
    int n1 = 0;
    double maxSideLength = -1.0;
    int boundary = -1;
    bool hidden = false;
    int group = 0;
};

struct TileArc {
    int n0 = 0;
    int n1 = 0;
    double arcLength = 0.0;
    double maxSegDegrees = 0.0;
    int boundary = -1;
    bool hidden = false;
    int group = 0;
};

/** A region label. References materials/circuits/boundaries by index. */
struct TileLabel {
    /** Stable name used by overrides; must be unique within the tile. */
    std::string name;
    double x = 0.0;
    double y = 0.0;
    /** Zero-based into TiledModel::materialProps; -1 means a hole/air. */
    int material = -1;
    /** Zero-based into TiledModel::circuitProps; -1 means no circuit. */
    int circuit = -1;
    int turns = 1;
    double magDir = 0.0;
    std::string magDirFctn;
    double maxArea = 0.0;
    int group = 0;
    bool hole = false;
};

struct TileGeometry {
    std::vector<TileNode> nodes;
    std::vector<TileSegment> segments;
    std::vector<TileArc> arcs;
    std::vector<TileLabel> labels;
};

enum class RepeatKind { Rotation };

/** How the modelled sector closes at the two ends of the instance sequence. */
enum class Closure {
    /** count * pitch == 360 degrees; the last instance welds to the first. */
    Closed,
    /** The end seams are identified by a periodic field constraint. */
    Periodic,
    /** The end seams are identified by an antiperiodic field constraint. */
    Antiperiodic
};

struct TileRepeat {
    RepeatKind kind = RepeatKind::Rotation;
    double centerXMetres = 0.0;
    double centerYMetres = 0.0;
    /** Number of instances; the pitch itself is derived from the geometry. */
    std::size_t count = 1;
    Closure closure = Closure::Closed;
};

struct Tile {
    std::string name;
    TileGeometry geometry;
    TileRepeat repeat;
    /** Boundary property name whose chains are the tile's two radial seams. */
    std::string seamBoundary;
    /** Boundary property name of the air-gap seam, if this tile has one. */
    std::string airGapBoundary;
};

struct TileCoupling {
    std::string innerTile;
    std::string outerTile;
    std::string boundary;
    double centerXMetres = 0.0;
    double centerYMetres = 0.0;
    double innerRadiusMetres = 0.0;
    double outerRadiusMetres = 0.0;
};

/** How a per-instance magnetisation array is interpreted. */
enum class MagnetisationMode {
    /** Added to the prototype MagDir after the instance transform rotation. */
    Delta,
    /** Used directly as the absolute magnetisation angle. */
    Absolute
};

/**
 * Per-label, per-instance physics. Empty arrays mean "not overridden".
 * Every non-empty array must have one entry per tile instance.
 */
struct TileLabelOverride {
    std::string tile;
    std::string label;
    std::vector<double> magDir;
    MagnetisationMode magDirMode = MagnetisationMode::Delta;
    std::vector<int> circuit;      // -1 = no circuit
    std::vector<double> turnScale;
};

/** The complete tiled magnetic problem: shared properties plus tiles. */
struct TiledModel {
    FileType fileType = FileType::MagneticsFile;
    double frequency = 0.0;
    double precision = 1e-8;
    double minAngle = 30.0;
    double depth = 1.0;
    LengthUnit lengthUnits = LengthMeters;
    ProblemType problemType = PLANAR;
    CoordsType coords = CART;

    std::vector<std::unique_ptr<CMPointProp>> pointProps;
    std::vector<std::unique_ptr<CMBoundaryProp>> boundaryProps;
    std::vector<std::unique_ptr<CMMaterialProp>> materialProps;
    std::vector<std::unique_ptr<CMCircuit>> circuitProps;

    std::vector<Tile> tiles;
    std::vector<TileCoupling> couplings;
    std::vector<TileLabelOverride> overrides;
};

enum class TiledDiagnosticCategory {
    FileFormatError,
    UnsupportedVersion,
    DuplicateTileName,
    UnknownTile,
    UnknownLabel,
    UnknownBoundary,
    UnknownMaterial,
    UnknownCircuit,
    SeamBoundaryNotPeriodic,
    AirGapBoundaryNotAge,
    InvalidRepeatKind,
    InvalidCount,
    SeamArcUnsupported,
    SeamNotRadial,
    SeamCountMismatch,
    ClosureCountMismatch,
    GeometryOutsideSector,
    InvalidGeometryIndex,
    NonFiniteValue,
    OverrideLengthMismatch,
    OverrideValueOutOfRange,
    DuplicateOverride,
    CouplingTileMissing,
    CouplingBoundaryMismatch,
    CouplingCentreMismatch,
    CouplingRadiusInvalid
};

struct TiledDiagnostic {
    TiledDiagnosticCategory category = TiledDiagnosticCategory::NonFiniteValue;
    std::string tile;
    std::string object;
    std::string message;
};

/** Values derived from a valid tile's seam geometry. */
struct TileDerived {
    double pitchDegrees = 0.0;
    double leftSeamAngleDegrees = 0.0;
    double rightSeamAngleDegrees = 0.0;
    double totalAngleDegrees = 0.0;
};

struct TiledValidationResult {
    std::vector<TiledDiagnostic> diagnostics;
    /** Parallel to TiledModel::tiles; only meaningful when the tile is valid. */
    std::vector<TileDerived> tiles;

    bool succeeded() const { return diagnostics.empty(); }
};

/** Tolerance for straight-radial-seam and centre checks, in metres. */
constexpr double TiledModelLengthToleranceMetres = 1e-9;
/** Tolerance for the derived pitch and coverage checks, in degrees. */
constexpr double TiledModelAngleToleranceDegrees = 1e-7;

/**
 * Validate the tiling semantics: seam boundaries are straight radial lines
 * through the rotation centre, the pitch is derived from the two seam chains,
 * count*pitch covers exactly the declared closure, no tile geometry lies
 * outside its sector, overrides have one entry per instance, and couplings are
 * consistent. Returns one entry per tile in \c tiles.
 */
TiledValidationResult validateTiledModel(const TiledModel &model);

/** Find a tile by name, or nullptr. */
const Tile *findTile(const TiledModel &model, const std::string &name);

/**
 * Build a standalone magnetic FemmProblem for one tile: the shared properties
 * plus the tile's geometry. The problem deliberately has no pathName, so a
 * mesher can treat it as an in-memory problem (Tangle's record entry point).
 */
std::unique_ptr<FemmProblem> buildTileProblem(const TiledModel &model,
                                              std::size_t tileIndex);

/** Index of a named boundary/material/circuit/point property, or -1. */
int boundaryIndexByName(const TiledModel &model, const std::string &name);
int materialIndexByName(const TiledModel &model, const std::string &name);
int circuitIndexByName(const TiledModel &model, const std::string &name);

/** True when the boundary property at \p index is periodic/antiperiodic. */
bool boundaryIsPeriodic(const TiledModel &model, int index);
/** True when the boundary property at \p index is an air-gap element. */
bool boundaryIsAirGap(const TiledModel &model, int index);

} // namespace tiled
} // namespace femm

#endif
