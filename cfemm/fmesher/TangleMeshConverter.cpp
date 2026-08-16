#include "TangleMeshConverter.h"

#include <tangle_mesh.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <set>
#include <utility>
#include <vector>

namespace fmesher {
namespace {

using femm::mesh::MeshDiagnosticSeverity;
using femm::mesh::MeshIndex;
using femm::mesh::MeshResult;
using femm::mesh::MeshStatus;
using femm::mesh::SolverMesh;

constexpr double Pi = 3.141592653589793238462643383279502884;

MeshResult failure(const std::string &message)
{
    MeshResult result;
    result.status = MeshStatus::ConversionFailure;
    result.diagnostics.push_back(
        {MeshDiagnosticSeverity::Error, message, "Tangle", 0});
    return result;
}

bool validIndex(int index, std::size_t count)
{
    return index >= 0 && static_cast<std::size_t>(index) < count;
}

bool finite(double value)
{
    return std::isfinite(value);
}

std::pair<int, int> edgeKey(int first, int second)
{
    return first < second ? std::make_pair(first, second)
                          : std::make_pair(second, first);
}

bool regionAttribute(double value, std::int32_t &converted)
{
    if (!finite(value) || value != std::trunc(value) ||
        value < static_cast<double>(std::numeric_limits<std::int32_t>::min()) ||
        value > static_cast<double>(std::numeric_limits<std::int32_t>::max()))
        return false;
    converted = static_cast<std::int32_t>(value);
    return true;
}

SolverMesh::Periodicity periodicity(int type)
{
    return type == 0 ? SolverMesh::Periodicity::Periodic
                     : SolverMesh::Periodicity::Antiperiodic;
}

struct RingPoint {
    MeshIndex node = femm::mesh::InvalidMeshIndex;
    double absoluteElementPosition = 0.0;
    double sign = 1.0;
};

bool buildAirGap(const ::Mesh &source, const AGEDef &definition,
                 double lengthScaleToMetres, SolverMesh::AirGap &target,
                 std::string &error)
{
    const std::size_t nodeCount = source.vertices.size();
    if (definition.format != 0 && definition.format != 1) {
        error = "AGE '" + definition.name + "' has an invalid periodicity type";
        return false;
    }
    if (definition.n <= 0 || definition.innerNodes.size() != definition.outerNodes.size() ||
        definition.innerNodes.size() != static_cast<std::size_t>(definition.n)) {
        error = "AGE '" + definition.name + "' has inconsistent ring node counts";
        return false;
    }
    if (!finite(definition.innerAngle) || !finite(definition.outerAngle) ||
        !finite(definition.ri) || !finite(definition.ro) ||
        !finite(definition.totalArcLength) || !finite(definition.cx) ||
        !finite(definition.cy) || definition.ri < 0.0 ||
        definition.ro <= definition.ri || definition.totalArcLength <= 0.0) {
        error = "AGE '" + definition.name + "' has invalid geometry";
        return false;
    }
    for (int node : definition.innerNodes)
        if (!validIndex(node, nodeCount)) {
            error = "AGE '" + definition.name + "' has an invalid inner-ring node";
            return false;
        }
    for (int node : definition.outerNodes)
        if (!validIndex(node, nodeCount)) {
            error = "AGE '" + definition.name + "' has an invalid outer-ring node";
            return false;
        }

    const double spacing = definition.totalArcLength / definition.n;
    const double fullElementCount = 360.0 / spacing;
    const double copyCount = 360.0 / definition.totalArcLength;
    const auto totalElements = static_cast<std::size_t>(std::llround(fullElementCount));
    const auto copies = static_cast<std::size_t>(std::llround(copyCount));
    const double tolerance = 1e-9;
    if (totalElements == 0 || copies == 0 ||
        std::abs(fullElementCount - totalElements) > tolerance ||
        std::abs(copyCount - copies) > tolerance ||
        totalElements != copies * static_cast<std::size_t>(definition.n) ||
        (definition.format == 1 && copies % 2 != 0)) {
        error = "AGE '" + definition.name + "' does not tile a full annular ring";
        return false;
    }

    std::vector<RingPoint> innerRing;
    std::vector<RingPoint> outerRing;
    innerRing.reserve(totalElements);
    outerRing.reserve(totalElements);
    for (std::size_t copy = 0; copy < copies; ++copy) {
        const double sign = definition.format == 1 && copy % 2 != 0 ? -1.0 : 1.0;
        const double innerRotation =
            (copy * definition.totalArcLength + definition.innerAngle) * Pi / 180.0;
        const double outerRotation =
            (copy * definition.totalArcLength + definition.outerAngle) * Pi / 180.0;
        const double innerCos = std::cos(innerRotation);
        const double innerSin = std::sin(innerRotation);
        const double outerCos = std::cos(outerRotation);
        const double outerSin = std::sin(outerRotation);
        for (int i = 0; i < definition.n; ++i) {
            const int innerNode = definition.innerNodes[static_cast<std::size_t>(i)];
            double dx = source.vertices[static_cast<std::size_t>(innerNode)].x - definition.cx;
            double dy = source.vertices[static_cast<std::size_t>(innerNode)].y - definition.cy;
            double angle = std::atan2(innerSin * dx + innerCos * dy,
                                      innerCos * dx - innerSin * dy) * 180.0 / Pi;
            if (angle < 0.0)
                angle += 360.0;
            innerRing.push_back(
                {static_cast<MeshIndex>(innerNode), angle / spacing, sign});

            const int outerNode = definition.outerNodes[static_cast<std::size_t>(i)];
            dx = source.vertices[static_cast<std::size_t>(outerNode)].x - definition.cx;
            dy = source.vertices[static_cast<std::size_t>(outerNode)].y - definition.cy;
            angle = std::atan2(outerSin * dx + outerCos * dy,
                               outerCos * dx - outerSin * dy) * 180.0 / Pi;
            if (angle < 0.0)
                angle += 360.0;
            outerRing.push_back(
                {static_cast<MeshIndex>(outerNode), angle / spacing, sign});
        }
    }
    const auto byPosition = [](const RingPoint &left, const RingPoint &right) {
        return left.absoluteElementPosition < right.absoluteElementPosition;
    };
    std::stable_sort(innerRing.begin(), innerRing.end(), byPosition);
    std::stable_sort(outerRing.begin(), outerRing.end(), byPosition);

    target.boundaryName = definition.name;
    target.periodicity = periodicity(definition.format);
    target.totalArcElements = static_cast<std::size_t>(definition.n);
    target.totalArcLengthDegrees = definition.totalArcLength;
    target.innerRadius = definition.ri * lengthScaleToMetres;
    target.outerRadius = definition.ro * lengthScaleToMetres;
    target.innerAngleDegrees = definition.innerAngle;
    target.outerAngleDegrees = definition.outerAngle;
    target.centerX = definition.cx * lengthScaleToMetres;
    target.centerY = definition.cy * lengthScaleToMetres;
    target.innerShift = innerRing.front().absoluteElementPosition;
    target.outerShift = outerRing.front().absoluteElementPosition;
    target.innerRing.reserve(totalElements);
    target.outerRing.reserve(totalElements);
    target.nodeIndices.reserve(totalElements * 2);
    for (const auto &point : innerRing) {
        target.nodeIndices.push_back(point.node);
        target.innerRing.push_back(
            {point.node,
             point.absoluteElementPosition - definition.innerAngle / spacing,
             point.sign});
    }
    for (const auto &point : outerRing) {
        target.nodeIndices.push_back(point.node);
        target.outerRing.push_back(
            {point.node,
             point.absoluteElementPosition - definition.outerAngle / spacing,
             point.sign});
    }
    target.quadraturePoints.reserve(static_cast<std::size_t>(definition.n + 1));
    for (int i = 0; i <= definition.n; ++i) {
        const std::size_t next = static_cast<std::size_t>(i) % totalElements;
        const std::size_t previous = next == 0 ? totalElements - 1 : next - 1;
        SolverMesh::AirGapQuadraturePoint point;
        point.nodes = {{innerRing[previous].node, innerRing[next].node,
                        outerRing[previous].node, outerRing[next].node}};
        point.weights = {{innerRing[previous].sign, innerRing[next].sign,
                          outerRing[previous].sign, outerRing[next].sign}};
        target.quadraturePoints.push_back(point);
    }
    return true;
}

} // namespace

MeshResult convertTangleMesh(const ::Mesh &source, double lengthScaleToMetres)
{
    if (!finite(lengthScaleToMetres) || lengthScaleToMetres <= 0.0)
        return failure("Tangle conversion requires a positive finite length scale");

    MeshResult result;
    result.status = MeshStatus::ConversionFailure;
    result.mesh.nodes.reserve(source.vertices.size());
    for (const auto &node : source.vertices) {
        if (!finite(node.x) || !finite(node.y))
            return failure("Tangle mesh contains a non-finite node coordinate");
        result.mesh.nodes.push_back(
            {node.x * lengthScaleToMetres, node.y * lengthScaleToMetres,
             static_cast<std::int32_t>(node.marker)});
    }

    result.mesh.elements.reserve(source.triangles.size());
    for (const auto &triangle : source.triangles) {
        for (int node : triangle.v)
            if (!validIndex(node, source.vertices.size()))
                return failure("Tangle triangle references an invalid node");
        std::int32_t attribute = 0;
        if (!regionAttribute(triangle.region_attrib, attribute))
            return failure("Tangle triangle has a non-integral or out-of-range region attribute");
        result.mesh.elements.push_back(
            {{{static_cast<MeshIndex>(triangle.v[0]),
               static_cast<MeshIndex>(triangle.v[1]),
               static_cast<MeshIndex>(triangle.v[2])}}, attribute});
    }

    std::map<std::pair<int, int>, int> segmentMarkers;
    for (const auto &segment : source.segments) {
        if (!validIndex(segment.v0, source.vertices.size()) ||
            !validIndex(segment.v1, source.vertices.size()))
            return failure("Tangle segment references an invalid node");
        const auto key = edgeKey(segment.v0, segment.v1);
        const auto existing = segmentMarkers.find(key);
        if (existing != segmentMarkers.end() && existing->second != segment.marker)
            return failure("Tangle mesh has conflicting markers for one edge");
        segmentMarkers[key] = segment.marker;
    }
    result.mesh.edges.reserve(source.edges.size());
    std::set<std::pair<int, int>> seenEdges;
    for (const auto &edge : source.edges) {
        if (!validIndex(edge.first, source.vertices.size()) ||
            !validIndex(edge.second, source.vertices.size()))
            return failure("Tangle edge references an invalid node");
        const auto key = edgeKey(edge.first, edge.second);
        if (!seenEdges.insert(key).second)
            return failure("Tangle mesh contains a duplicate edge");
        const auto marker = segmentMarkers.find(key);
        result.mesh.edges.push_back(
            {static_cast<MeshIndex>(edge.first), static_cast<MeshIndex>(edge.second),
             marker == segmentMarkers.end() ? 0 : marker->second});
    }

    result.mesh.periodicConstraints.reserve(source.pbc_pairs.size());
    for (const auto &pair : source.pbc_pairs) {
        if (!validIndex(pair.node_a, source.vertices.size()) ||
            !validIndex(pair.node_b, source.vertices.size()))
            return failure("Tangle periodic pair references an invalid node");
        if (pair.type != 0 && pair.type != 1)
            return failure("Tangle periodic pair has an invalid periodicity type");
        result.mesh.periodicConstraints.push_back(
            {static_cast<MeshIndex>(pair.node_a), static_cast<MeshIndex>(pair.node_b),
             periodicity(pair.type)});
    }

    result.mesh.airGaps.reserve(source.age_defs.size());
    for (const auto &definition : source.age_defs) {
        SolverMesh::AirGap airGap;
        std::string error;
        if (!buildAirGap(source, definition, lengthScaleToMetres, airGap, error))
            return failure(error);
        result.mesh.airGaps.push_back(std::move(airGap));
    }

    result.status = MeshStatus::Success;
    return result;
}

} // namespace fmesher
