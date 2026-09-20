#include "TiledModel.h"

#include "CArcSegment.h"
#include "CBlockLabel.h"
#include "CNode.h"
#include "CSegment.h"
#include "FemmProblem.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace femm {
namespace tiled {
namespace {

constexpr double Pi = 3.141592653589793238462643383279502884;

double wrap360(double angle)
{
    angle = std::fmod(angle, 360.0);
    if (angle < 0.0)
        angle += 360.0;
    return angle;
}

double angularDistance(double a, double b)
{
    const double d = std::fabs(wrap360(a - b));
    return d > 180.0 ? 360.0 - d : d;
}

bool finite(double value) { return std::isfinite(value); }

void addDiagnostic(TiledValidationResult &result, TiledDiagnosticCategory category,
                   const std::string &tile, const std::string &object,
                   const std::string &message)
{
    result.diagnostics.push_back({category, tile, object, message});
}

struct Point {
    double x = 0.0;
    double y = 0.0;
};

Point nodePoint(const TileGeometry &geometry, int index)
{
    const TileNode &node = geometry.nodes[static_cast<std::size_t>(index)];
    return {node.x, node.y};
}

bool validIndex(int index, std::size_t size)
{
    return index >= 0 && static_cast<std::size_t>(index) < size;
}

/**
 * Cluster the seam segments by the outward radial direction of their midpoint.
 * Returns the cluster angles in degrees. A straight radial seam contributes all
 * of its collinear segments to one cluster.
 */
std::vector<double> seamClusterAngles(const TileGeometry &geometry,
                                      const std::vector<std::size_t> &seamSegments,
                                      const Point &center)
{
    constexpr double GroupToleranceDegrees = 1e-6;
    std::vector<double> angles;
    angles.reserve(seamSegments.size());
    for (std::size_t segment : seamSegments) {
        const TileSegment &s = geometry.segments[segment];
        const Point a = nodePoint(geometry, s.n0);
        const Point b = nodePoint(geometry, s.n1);
        const double mx = 0.5 * (a.x + b.x) - center.x;
        const double my = 0.5 * (a.y + b.y) - center.y;
        angles.push_back(wrap360(std::atan2(my, mx) * 180.0 / Pi));
    }

    std::vector<bool> used(angles.size(), false);
    std::vector<double> clusters;
    for (std::size_t i = 0; i < angles.size(); ++i) {
        if (used[i])
            continue;
        used[i] = true;
        clusters.push_back(angles[i]);
        for (std::size_t j = i + 1; j < angles.size(); ++j) {
            if (!used[j] && angularDistance(angles[i], angles[j]) <= GroupToleranceDegrees)
                used[j] = true;
        }
    }
    return clusters;
}

/** True when the point lies inside [leftAngle, leftAngle+pitch]. */
bool pointInsideSector(const Point &point, const Point &center, double leftAngle,
                       double pitch)
{
    const double dx = point.x - center.x;
    const double dy = point.y - center.y;
    if (std::hypot(dx, dy) <= TiledModelLengthToleranceMetres)
        return true;
    const double angle = std::atan2(dy, dx) * 180.0 / Pi;
    const double relative = wrap360(angle - leftAngle);
    const double slack = TiledModelAngleToleranceDegrees;
    return relative <= pitch + slack || relative >= 360.0 - slack;
}

/** True when every tile node lies inside [leftAngle, leftAngle+pitch]. */
bool coversSector(const TileGeometry &geometry, const Point &center,
                  double leftAngle, double pitch)
{
    for (const TileNode &node : geometry.nodes)
        if (!pointInsideSector({node.x, node.y}, center, leftAngle, pitch))
            return false;
    return true;
}

/**
 * A representative interior point. A region label is a guaranteed interior
 * seed, so prefer one; fall back to the node centroid. It disambiguates which
 * of the two sectors between the seams the tile occupies.
 */
Point referencePoint(const TileGeometry &geometry)
{
    if (!geometry.labels.empty()) {
        double x = 0.0;
        double y = 0.0;
        for (const TileLabel &label : geometry.labels) {
            x += label.x;
            y += label.y;
        }
        const double n = static_cast<double>(geometry.labels.size());
        return {x / n, y / n};
    }
    double x = 0.0;
    double y = 0.0;
    for (const TileNode &node : geometry.nodes) {
        x += node.x;
        y += node.y;
    }
    const double n = static_cast<double>(geometry.nodes.size());
    return n > 0.0 ? Point{x / n, y / n} : Point{};
}

} // namespace

int boundaryIndexByName(const TiledModel &model, const std::string &name)
{
    for (std::size_t i = 0; i < model.boundaryProps.size(); ++i)
        if (model.boundaryProps[i] && model.boundaryProps[i]->BdryName == name)
            return static_cast<int>(i);
    return -1;
}

int materialIndexByName(const TiledModel &model, const std::string &name)
{
    for (std::size_t i = 0; i < model.materialProps.size(); ++i)
        if (model.materialProps[i] && model.materialProps[i]->BlockName == name)
            return static_cast<int>(i);
    return -1;
}

int circuitIndexByName(const TiledModel &model, const std::string &name)
{
    for (std::size_t i = 0; i < model.circuitProps.size(); ++i)
        if (model.circuitProps[i] && model.circuitProps[i]->CircName == name)
            return static_cast<int>(i);
    return -1;
}

bool boundaryIsPeriodic(const TiledModel &model, int index)
{
    if (!validIndex(index, model.boundaryProps.size()))
        return false;
    const CMBoundaryProp *prop = model.boundaryProps[static_cast<std::size_t>(index)].get();
    return prop && prop->isPeriodic(CBoundaryProp::PeriodicityType::Any);
}

bool boundaryIsAirGap(const TiledModel &model, int index)
{
    if (!validIndex(index, model.boundaryProps.size()))
        return false;
    const CMBoundaryProp *prop = model.boundaryProps[static_cast<std::size_t>(index)].get();
    return prop && (prop->BdryFormat == 6 || prop->BdryFormat == 7);
}

const Tile *findTile(const TiledModel &model, const std::string &name)
{
    for (const Tile &tile : model.tiles)
        if (tile.name == name)
            return &tile;
    return nullptr;
}

TiledValidationResult validateTiledModel(const TiledModel &model)
{
    TiledValidationResult result;
    result.tiles.resize(model.tiles.size());

    // Shared property tables and problem header.
    if (!finite(model.depth) || model.depth <= 0.0)
        addDiagnostic(result, TiledDiagnosticCategory::NonFiniteValue, {}, {},
                      "problem depth must be finite and positive");
    if (!finite(model.frequency) || !finite(model.precision) || !finite(model.minAngle))
        addDiagnostic(result, TiledDiagnosticCategory::NonFiniteValue, {}, {},
                      "problem header contains a non-finite value");

    // Tile names must be unique.
    for (std::size_t i = 0; i < model.tiles.size(); ++i) {
        for (std::size_t j = i + 1; j < model.tiles.size(); ++j) {
            if (!model.tiles[i].name.empty() && model.tiles[i].name == model.tiles[j].name)
                addDiagnostic(result, TiledDiagnosticCategory::DuplicateTileName,
                              model.tiles[i].name, {},
                              "tile name is used more than once");
        }
    }

    for (std::size_t tileIndex = 0; tileIndex < model.tiles.size(); ++tileIndex) {
        const Tile &tile = model.tiles[tileIndex];
        TileDerived &derived = result.tiles[tileIndex];
        const Point center{tile.repeat.centerXMetres, tile.repeat.centerYMetres};

        if (tile.name.empty())
            addDiagnostic(result, TiledDiagnosticCategory::UnknownTile, {}, {},
                          "tile name must not be empty");
        if (tile.repeat.kind != RepeatKind::Rotation)
            addDiagnostic(result, TiledDiagnosticCategory::InvalidRepeatKind, tile.name, {},
                          "only rotational repeat is supported");
        if (!finite(tile.repeat.centerXMetres) || !finite(tile.repeat.centerYMetres))
            addDiagnostic(result, TiledDiagnosticCategory::NonFiniteValue, tile.name, {},
                          "rotation centre must be finite");

        // Geometry index and finiteness checks.
        const std::size_t nodeCount = tile.geometry.nodes.size();
        for (const TileNode &node : tile.geometry.nodes)
            if (!finite(node.x) || !finite(node.y))
                addDiagnostic(result, TiledDiagnosticCategory::NonFiniteValue, tile.name, {},
                              "tile node coordinate is not finite");
        for (const TileSegment &s : tile.geometry.segments)
            if (!validIndex(s.n0, nodeCount) || !validIndex(s.n1, nodeCount))
                addDiagnostic(result, TiledDiagnosticCategory::InvalidGeometryIndex, tile.name, {},
                              "segment references an out-of-range node");
        for (const TileArc &a : tile.geometry.arcs)
            if (!validIndex(a.n0, nodeCount) || !validIndex(a.n1, nodeCount))
                addDiagnostic(result, TiledDiagnosticCategory::InvalidGeometryIndex, tile.name, {},
                              "arc references an out-of-range node");
        for (const TileLabel &label : tile.geometry.labels) {
            if (label.name.empty())
                addDiagnostic(result, TiledDiagnosticCategory::UnknownLabel, tile.name, {},
                              "tile label name must not be empty");
            if (!finite(label.x) || !finite(label.y) || !finite(label.magDir) ||
                !finite(label.maxArea))
                addDiagnostic(result, TiledDiagnosticCategory::NonFiniteValue, tile.name,
                              label.name, "tile label contains a non-finite value");
            if (label.material >= 0 &&
                !validIndex(label.material, model.materialProps.size()))
                addDiagnostic(result, TiledDiagnosticCategory::UnknownMaterial, tile.name,
                              label.name, "label references an unknown material");
            if (label.circuit >= 0 &&
                !validIndex(label.circuit, model.circuitProps.size()))
                addDiagnostic(result, TiledDiagnosticCategory::UnknownCircuit, tile.name,
                              label.name, "label references an unknown circuit");
        }
        // Label names must be unique within the tile (they key overrides).
        for (std::size_t i = 0; i < tile.geometry.labels.size(); ++i)
            for (std::size_t j = i + 1; j < tile.geometry.labels.size(); ++j)
                if (!tile.geometry.labels[i].name.empty() &&
                    tile.geometry.labels[i].name == tile.geometry.labels[j].name)
                    addDiagnostic(result, TiledDiagnosticCategory::UnknownLabel, tile.name,
                                  tile.geometry.labels[i].name,
                                  "tile label name is used more than once");

        // Seam derivation.
        const int seamIndex = boundaryIndexByName(model, tile.seamBoundary);
        if (seamIndex < 0) {
            addDiagnostic(result, TiledDiagnosticCategory::UnknownBoundary, tile.name,
                          tile.seamBoundary, "seam boundary property does not exist");
        } else if (!boundaryIsPeriodic(model, seamIndex)) {
            addDiagnostic(result, TiledDiagnosticCategory::SeamBoundaryNotPeriodic, tile.name,
                          tile.seamBoundary, "seam boundary must be periodic/antiperiodic");
        }

        std::vector<std::size_t> seamSegments;
        for (std::size_t s = 0; s < tile.geometry.segments.size(); ++s)
            if (seamIndex >= 0 && tile.geometry.segments[s].boundary == seamIndex)
                seamSegments.push_back(s);
        for (const TileArc &a : tile.geometry.arcs)
            if (seamIndex >= 0 && a.boundary == seamIndex)
                addDiagnostic(result, TiledDiagnosticCategory::SeamArcUnsupported, tile.name, {},
                              "a seam must be built from straight radial segments");

        // A radial seam segment must lie on a line through the rotation centre.
        for (std::size_t segment : seamSegments) {
            const TileSegment &s = tile.geometry.segments[segment];
            if (!validIndex(s.n0, nodeCount) || !validIndex(s.n1, nodeCount))
                continue;
            const Point a = nodePoint(tile.geometry, s.n0);
            const Point b = nodePoint(tile.geometry, s.n1);
            const double dx = b.x - a.x;
            const double dy = b.y - a.y;
            const double length = std::hypot(dx, dy);
            if (length <= TiledModelLengthToleranceMetres) {
                addDiagnostic(result, TiledDiagnosticCategory::NonFiniteValue, tile.name, {},
                              "degenerate seam segment");
                continue;
            }
            const double offset = std::fabs(dx * (a.y - center.y) - dy * (a.x - center.x)) / length;
            if (offset > TiledModelLengthToleranceMetres)
                addDiagnostic(result, TiledDiagnosticCategory::SeamNotRadial, tile.name, {},
                              "seam segment does not pass through the rotation centre");
        }

        const std::vector<double> clusters = seamClusterAngles(tile.geometry, seamSegments, center);
        if (clusters.size() != 2) {
            addDiagnostic(result, TiledDiagnosticCategory::SeamCountMismatch, tile.name, {},
                          "a tile requires exactly two distinct radial seams");
            continue;
        }

        const double count = static_cast<double>(tile.repeat.count);
        // Choose the assignment whose pitch contains the tile geometry. The
        // interior reference point resolves the reflex-sector ambiguity.
        const Point reference = referencePoint(tile.geometry);
        double leftAngle = clusters[0];
        double rightAngle = clusters[1];
        double pitch = wrap360(rightAngle - leftAngle);
        if (!pointInsideSector(reference, center, leftAngle, pitch)) {
            leftAngle = clusters[1];
            rightAngle = clusters[0];
            pitch = wrap360(rightAngle - leftAngle);
        }
        if (!coversSector(tile.geometry, center, leftAngle, pitch)) {
            addDiagnostic(result, TiledDiagnosticCategory::GeometryOutsideSector, tile.name, {},
                          "tile geometry does not lie within its seam sector");
            continue;
        }
        if (!(pitch > 0.0)) {
            addDiagnostic(result, TiledDiagnosticCategory::SeamCountMismatch, tile.name, {},
                          "derived tile pitch must be positive");
            continue;
        }

        const double total = count * pitch;
        if (tile.repeat.closure == Closure::Closed) {
            if (tile.repeat.count < 2 ||
                std::fabs(total - 360.0) > TiledModelAngleToleranceDegrees)
                addDiagnostic(result, TiledDiagnosticCategory::ClosureCountMismatch, tile.name, {},
                              "a closed tile requires count*pitch == 360 degrees");
        } else {
            if (tile.repeat.count < 1 ||
                total >= 360.0 - TiledModelAngleToleranceDegrees)
                addDiagnostic(result, TiledDiagnosticCategory::ClosureCountMismatch, tile.name, {},
                              "an open tile requires 0 < count*pitch < 360 degrees");
        }

        derived.pitchDegrees = pitch;
        derived.leftSeamAngleDegrees = leftAngle;
        derived.rightSeamAngleDegrees = rightAngle;
        derived.totalAngleDegrees = total;

        // Optional air-gap boundary.
        if (!tile.airGapBoundary.empty()) {
            const int gapIndex = boundaryIndexByName(model, tile.airGapBoundary);
            if (gapIndex < 0)
                addDiagnostic(result, TiledDiagnosticCategory::UnknownBoundary, tile.name,
                              tile.airGapBoundary, "air-gap boundary property does not exist");
            else if (!boundaryIsAirGap(model, gapIndex))
                addDiagnostic(result, TiledDiagnosticCategory::AirGapBoundaryNotAge, tile.name,
                              tile.airGapBoundary, "air-gap boundary must be an AGE property");
        }
    }

    // Overrides.
    for (std::size_t i = 0; i < model.overrides.size(); ++i) {
        const TileLabelOverride &override = model.overrides[i];
        const Tile *tile = findTile(model, override.tile);
        if (!tile) {
            addDiagnostic(result, TiledDiagnosticCategory::UnknownTile, override.tile, override.label,
                          "override references an unknown tile");
            continue;
        }
        const TileLabel *label = nullptr;
        for (const TileLabel &candidate : tile->geometry.labels)
            if (candidate.name == override.label) {
                label = &candidate;
                break;
            }
        if (!label)
            addDiagnostic(result, TiledDiagnosticCategory::UnknownLabel, override.tile,
                          override.label, "override references an unknown tile label");

        const std::size_t count = tile->repeat.count;
        const auto checkLength = [&](std::size_t length, const char *field) {
            if (length != 0 && length != count)
                addDiagnostic(result, TiledDiagnosticCategory::OverrideLengthMismatch,
                              override.tile, override.label,
                              std::string(field) + " must have one entry per instance");
        };
        checkLength(override.magDir.size(), "magDir");
        checkLength(override.circuit.size(), "circuit");
        checkLength(override.turnScale.size(), "turnScale");

        for (double value : override.magDir)
            if (!finite(value))
                addDiagnostic(result, TiledDiagnosticCategory::NonFiniteValue, override.tile,
                              override.label, "magDir override is not finite");
        for (double value : override.turnScale)
            if (!finite(value))
                addDiagnostic(result, TiledDiagnosticCategory::NonFiniteValue, override.tile,
                              override.label, "turnScale override is not finite");
        for (int circuit : override.circuit)
            if (circuit < -1 || circuit >= static_cast<int>(model.circuitProps.size()))
                addDiagnostic(result, TiledDiagnosticCategory::OverrideValueOutOfRange,
                              override.tile, override.label,
                              "circuit override is out of range");

        for (std::size_t j = i + 1; j < model.overrides.size(); ++j)
            if (model.overrides[j].tile == override.tile &&
                model.overrides[j].label == override.label)
                addDiagnostic(result, TiledDiagnosticCategory::DuplicateOverride, override.tile,
                              override.label, "duplicate override for this label");
    }

    // Couplings.
    for (const TileCoupling &coupling : model.couplings) {
        const Tile *inner = findTile(model, coupling.innerTile);
        const Tile *outer = findTile(model, coupling.outerTile);
        if (!inner || !outer || coupling.innerTile == coupling.outerTile) {
            addDiagnostic(result, TiledDiagnosticCategory::CouplingTileMissing, coupling.innerTile,
                          coupling.outerTile, "coupling references unknown or equal tiles");
            continue;
        }
        const int gapIndex = boundaryIndexByName(model, coupling.boundary);
        if (gapIndex < 0)
            addDiagnostic(result, TiledDiagnosticCategory::UnknownBoundary, coupling.innerTile,
                          coupling.boundary, "coupling boundary property does not exist");
        else if (!boundaryIsAirGap(model, gapIndex))
            addDiagnostic(result, TiledDiagnosticCategory::AirGapBoundaryNotAge, coupling.innerTile,
                          coupling.boundary, "coupling boundary must be an AGE property");
        if (inner->airGapBoundary != coupling.boundary ||
            outer->airGapBoundary != coupling.boundary)
            addDiagnostic(result, TiledDiagnosticCategory::CouplingBoundaryMismatch,
                          coupling.innerTile, coupling.boundary,
                          "coupled tiles must declare this air-gap boundary");
        const double centreTolerance = TiledModelLengthToleranceMetres;
        if (std::fabs(inner->repeat.centerXMetres - coupling.centerXMetres) > centreTolerance ||
            std::fabs(inner->repeat.centerYMetres - coupling.centerYMetres) > centreTolerance ||
            std::fabs(outer->repeat.centerXMetres - coupling.centerXMetres) > centreTolerance ||
            std::fabs(outer->repeat.centerYMetres - coupling.centerYMetres) > centreTolerance)
            addDiagnostic(result, TiledDiagnosticCategory::CouplingCentreMismatch,
                          coupling.innerTile, coupling.outerTile,
                          "coupling centre must match both tiles' rotation centre");
        if (!(coupling.innerRadiusMetres > 0.0) ||
            !(coupling.outerRadiusMetres > coupling.innerRadiusMetres))
            addDiagnostic(result, TiledDiagnosticCategory::CouplingRadiusInvalid,
                          coupling.innerTile, coupling.outerTile,
                          "coupling radii must satisfy 0 < inner < outer");
    }

    return result;
}

std::unique_ptr<FemmProblem> buildTileProblem(const TiledModel &model,
                                              std::size_t tileIndex)
{
    if (tileIndex >= model.tiles.size())
        return nullptr;
    const Tile &tile = model.tiles[tileIndex];

    auto problem = std::make_unique<FemmProblem>(FileType::MagneticsFile);
    problem->Frequency = model.frequency;
    problem->Precision = model.precision;
    problem->MinAngle = model.minAngle;
    problem->Depth = model.depth;
    problem->LengthUnits = model.lengthUnits;
    problem->problemType = model.problemType;
    problem->Coords = model.coords;

    for (const auto &prop : model.pointProps)
        if (prop)
            problem->nodeproplist.push_back(std::make_unique<CMPointProp>(*prop));
    for (const auto &prop : model.boundaryProps)
        if (prop)
            problem->lineproplist.push_back(std::make_unique<CMBoundaryProp>(*prop));
    for (const auto &prop : model.materialProps)
        if (prop)
            problem->blockproplist.push_back(std::make_unique<CMMaterialProp>(*prop));
    for (const auto &prop : model.circuitProps)
        if (prop)
            problem->circproplist.push_back(std::make_unique<CMCircuit>(*prop));

    for (const TileNode &node : tile.geometry.nodes) {
        auto converted = std::make_unique<CNode>(node.x, node.y);
        converted->BoundaryMarker = node.boundaryMarker;
        converted->InGroup = node.group;
        problem->nodelist.push_back(std::move(converted));
    }
    for (const TileSegment &segment : tile.geometry.segments) {
        auto converted = std::make_unique<CSegment>();
        converted->n0 = segment.n0;
        converted->n1 = segment.n1;
        converted->MaxSideLength = segment.maxSideLength;
        converted->BoundaryMarker = segment.boundary;
        converted->Hidden = segment.hidden;
        converted->InGroup = segment.group;
        problem->linelist.push_back(std::move(converted));
    }
    for (const TileArc &arc : tile.geometry.arcs) {
        auto converted = std::make_unique<CArcSegment>();
        converted->n0 = arc.n0;
        converted->n1 = arc.n1;
        converted->ArcLength = arc.arcLength;
        converted->MaxSideLength = arc.maxSegDegrees;
        converted->BoundaryMarker = arc.boundary;
        converted->Hidden = arc.hidden;
        converted->InGroup = arc.group;
        problem->arclist.push_back(std::move(converted));
    }
    for (const TileLabel &label : tile.geometry.labels) {
        auto converted = std::make_unique<CMBlockLabel>();
        converted->x = label.x;
        converted->y = label.y;
        converted->BlockType = label.material;
        converted->InCircuit = label.circuit;
        converted->MaxArea = label.maxArea;
        converted->MagDir = label.magDir;
        converted->MagDirFctn = label.magDirFctn;
        converted->InGroup = label.group;
        converted->Turns = label.turns;
        converted->IsExternal = false;
        converted->IsDefault = false;
        problem->labellist.push_back(std::move(converted));
    }
    problem->updateLabelsFromIndex();
    problem->updateBlockMap();
    problem->updateCircuitMap();
    problem->updateLineMap();
    problem->updateNodeMap();
    return problem;
}

} // namespace tiled
} // namespace femm
