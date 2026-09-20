/*
   Native compressed magnetostatic assembly for xfemm.

   This file assembles a planar, static magnetostatic problem directly from a
   compact mesh::LogicalMeshView. It mirrors the physics of FSolver::Static2D in
   static2d.cpp, but iterates logical instances instead of a materialised element
   list, so it never stores the fully expanded mesh. The materialised Static2D
   path remains the oracle used by the equivalence tests.

   Scope of the native path (Milestone G):
     - planar coordinates;
     - linear and nonlinear isotropic materials (B-H curve), including the
       lamination types handled by Static2D;
     - constant magnetisation direction (no Lua MagDirFctn);
     - prescribed-current circuits, point currents, and line boundary
       conditions 0 and 2;
     - ordinary periodic constraints and air-gap elements.

   Incremental/frozen previous solutions, functional magnetisation, harmonic
   problems, and axisymmetric coordinates are rejected with a zero return value
   so a caller can fall back to the materialised path.

   Per-logical-element permeability (mu1, mu2, v12) and element coordinates are
   solver working state, not part of the stored mesh; the stored topology stays
   compressed in the LogicalMeshView.
*/

#include "femmconstants.h"
#include "femmcomplex.h"
#include "fsolver.h"
#include "linsolve/LinearSystemBackend.h"
#include "lua.h"
#include "LuaInstance.h"
#include "mesh/LogicalMeshView.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

using femm::mesh::MeshIndex;

namespace {

double nativePower(double x, int y)
{
    return std::pow(x, static_cast<double>(y));
}

/** Build the solver-side AGE list from a logical view's remapped AGE topology. */
std::vector<femmsolver::CAirGapElement>
buildAirGapElements(const femm::mesh::LogicalMeshView &view)
{
    std::vector<femmsolver::CAirGapElement> ages;
    ages.reserve(view.airGaps().size());
    for (const auto &source : view.airGaps()) {
        if (source.quadraturePoints.size() != source.totalArcElements + 1)
            return {};
        femmsolver::CAirGapElement age;
        age.BdryName = source.boundaryName;
        age.BdryFormat =
            source.periodicity == femm::mesh::SolverMesh::Periodicity::Antiperiodic ? 1 : 0;
        age.InnerAngle = source.innerAngleDegrees;
        age.OuterAngle = source.outerAngleDegrees;
        age.ri = source.innerRadius * 100.0;
        age.ro = source.outerRadius * 100.0;
        age.totalArcLength = source.totalArcLengthDegrees;
        age.agc = CComplex(source.centerX * 100.0, source.centerY * 100.0);
        age.totalArcElements = static_cast<int>(source.totalArcElements);
        age.InnerShift = source.innerShift;
        age.OuterShift = source.outerShift;
        for (const auto &sourcePoint : source.quadraturePoints) {
            femm::CQuadPoint point;
            point.n0 = static_cast<int>(sourcePoint.nodes[0]);
            point.n1 = static_cast<int>(sourcePoint.nodes[1]);
            point.n2 = static_cast<int>(sourcePoint.nodes[2]);
            point.n3 = static_cast<int>(sourcePoint.nodes[3]);
            point.w0 = sourcePoint.weights[0];
            point.w1 = sourcePoint.weights[1];
            point.w2 = sourcePoint.weights[2];
            point.w3 = sourcePoint.weights[3];
            age.quadNode.push_back(point);
        }
        for (const auto &ringPoint : source.innerRing) {
            femm::CQuadPoint point;
            point.n0 = static_cast<int>(ringPoint.node);
            point.w0 = ringPoint.elementPosition;
            point.w1 = ringPoint.weight;
            age.innerRingTopology.push_back(point);
        }
        for (const auto &ringPoint : source.outerRing) {
            femm::CQuadPoint point;
            point.n0 = static_cast<int>(ringPoint.node);
            point.w0 = ringPoint.elementPosition;
            point.w1 = ringPoint.weight;
            age.outerRingTopology.push_back(point);
        }
        for (auto index : source.nodeIndices)
            age.nodeNums.push_back(static_cast<int>(index));
        ages.push_back(std::move(age));
    }
    return ages;
}

/**
 * Reposition AGE quadrature nodes for the requested inner/outer angles,
 * mirroring FSolverAnalysisBackend::positionAirGaps without touching the
 * stored mesh topology.
 */
void applyAirGapPositions(
    std::vector<femmsolver::CAirGapElement> &ages,
    const std::map<std::string, std::pair<double, double>> &positions)
{
    for (auto &age : ages) {
        const auto requested = positions.find(age.BdryName);
        if (requested == positions.end())
            continue;
        const double innerAngle = requested->second.first;
        const double outerAngle = requested->second.second;
        if (age.InnerAngle == innerAngle && age.OuterAngle == outerAngle)
            continue;
        if (age.innerRingTopology.empty() || age.outerRingTopology.empty()) {
            age.InnerAngle = innerAngle;
            age.OuterAngle = outerAngle;
            continue;
        }
        const double step = age.totalArcLength / age.totalArcElements;
        const auto positioned = [step](const std::vector<femm::CQuadPoint> &topology,
                                       double angle) {
            auto ring = topology;
            for (auto &point : ring) {
                point.w0 = std::fmod(point.w0 * step + angle, 360.0);
                if (point.w0 < 0)
                    point.w0 += 360.0;
                point.w0 /= step;
            }
            std::stable_sort(ring.begin(), ring.end(),
                             [](const femm::CQuadPoint &left, const femm::CQuadPoint &right) {
                                 return left.w0 < right.w0;
                             });
            return ring;
        };
        const auto inner = positioned(age.innerRingTopology, innerAngle);
        const auto outer = positioned(age.outerRingTopology, outerAngle);
        const int fullCount = static_cast<int>(inner.size());
        age.InnerShift = inner.front().w0;
        age.OuterShift = outer.front().w0;
        age.quadNode.clear();
        age.quadNode.reserve(age.totalArcElements + 1);
        for (int i = 0; i <= age.totalArcElements; ++i) {
            const int p1 = i == fullCount ? 0 : i;
            const int p0 = p1 == 0 ? fullCount - 1 : p1 - 1;
            femm::CQuadPoint q;
            q.n0 = inner[p0].n0;
            q.n1 = inner[p1].n0;
            q.n2 = outer[p0].n0;
            q.n3 = outer[p1].n0;
            q.w0 = inner[p0].w1;
            q.w1 = inner[p1].w1;
            q.w2 = outer[p0].w1;
            q.w3 = outer[p1].w1;
            age.quadNode.push_back(q);
        }
        age.InnerAngle = innerAngle;
        age.OuterAngle = outerAngle;
    }
}

} // namespace

int FSolver::Static2DNative(const femm::mesh::LogicalMeshView &view,
                            const std::vector<std::size_t> &instanceLabelBase,
                            femm::LinearSystemBackend<double> &L,
                            const std::map<std::string, std::pair<double, double>>
                                &airGapPositions,
                            const std::vector<std::size_t> &nodePermutation)
{
    if (ProblemType != femm::PLANAR)
        return 0;
    if (Frequency != 0)
        return 0;
    if (!previousSolutionFile.empty())
        return 0;

    const std::size_t elementCount = view.elementCount();
    const std::size_t nodeCount = view.nodeCount();
    if (NumNodes != static_cast<int>(nodeCount) || NumEls != static_cast<int>(elementCount))
        return 0;

    // Optional Cuthill-McKee permutation: nodePermutation[old] is the new
    // degree-of-freedom index. An empty permutation keeps global node order.
    const bool permuted = !nodePermutation.empty();
    if (permuted && nodePermutation.size() != nodeCount)
        return 0;
    const auto mapped = [&](std::size_t index) -> int {
        return static_cast<int>(permuted ? nodePermutation[index] : index);
    };

    const double c = PI * 4.e-05;
    const double units[] = {2.54, 0.1, 1., 100., 0.00254, 1.e-04};

    int defaultLabel = -1;
    for (int i = 0; i < NumBlockLabels; ++i) {
        GetFillFactor(i);
        if (labellist[i].IsDefault)
            defaultLabel = i;
    }

    // Resolve the solver label, block property, and side markers of every
    // logical element once. Side markers use the same encoding as LoadMesh.
    std::vector<int> elementLabel(elementCount, -1);
    std::vector<int> elementBlock(elementCount, -1);
    std::vector<std::array<int, 3>> elementSides(elementCount);
    for (std::size_t e = 0; e < elementCount; ++e) {
        const femm::mesh::ElementProvenance provenance = view.elementProvenance(e);
        const std::int32_t attribute = view.elementRegionAttribute(e);
        int label = attribute - 1;
        if (attribute > 0 && provenance.instanceIndex + 1 < instanceLabelBase.size())
            label = static_cast<int>(instanceLabelBase[provenance.instanceIndex] + attribute) - 1;
        if (label < 0)
            label = defaultLabel;
        if (label < 0 || label >= NumBlockLabels)
            return 0;
        elementLabel[e] = label;
        elementBlock[e] = labellist[label].BlockType;
        if (elementBlock[e] < 0 || elementBlock[e] >= NumBlockProps)
            return 0;

        const auto markers = view.elementEdgeMarkers(e);
        for (std::size_t j = 0; j < 3; ++j) {
            const int side = markers[j] < 0 ? -(markers[j] + 2) : -1;
            elementSides[e][j] = (side >= 0 && side < NumLineProps) ? side : -1;
        }
    }

    // Element coordinates and permuted connectivity are computed once and
    // reused across Newton iterations instead of re-querying the view.
    std::vector<std::array<double, 3>> elementX(elementCount);
    std::vector<std::array<double, 3>> elementY(elementCount);
    std::vector<std::array<int, 3>> elementNodes(elementCount);
    for (std::size_t e = 0; e < elementCount; ++e) {
        const auto nodes = view.elementNodes(e);
        double x[3];
        double y[3];
        view.elementCoordinates(e, x, y);
        for (int k = 0; k < 3; ++k) {
            elementX[e][k] = x[k] * 100.0;
            elementY[e][k] = y[k] * 100.0;
            elementNodes[e][k] = mapped(nodes[k]);
        }
    }

    // Circuit preprocessing: derive Case/J/dV from the prescribed currents.
    std::vector<double> circuitInt1(NumCircProps, 0.0);
    std::vector<double> circuitInt2(NumCircProps, 0.0);
    std::vector<double> circuitInt3(NumCircProps, 0.0);
    for (std::size_t e = 0; e < elementCount; ++e) {
        const int label = elementLabel[e];
        if (label < 0 || labellist[label].InCircuit < 0)
            continue;
        const double *x = elementX[e].data();
        const double *y = elementY[e].data();
        const double p0 = y[1] - y[2];
        const double p1 = y[2] - y[0];
        const double q0 = x[2] - x[1];
        const double q1 = x[0] - x[2];
        const double area = (p0 * q1 - p1 * q0) / 2.0;
        double cduct = blockproplist[elementBlock[e]].Cduct;
        if (labellist[label].bIsWound)
            cduct = 0;
        const int circuit = labellist[label].InCircuit;
        circuitInt1[circuit] += area;
        circuitInt2[circuit] += area * cduct;
        circuitInt3[circuit] += blockproplist[elementBlock[e]].J.re * area * 100.;
    }
    for (int i = 0; i < NumCircProps; ++i) {
        if (circproplist[i].CircType == 0) {
            if (circuitInt2[i] == 0) {
                circproplist[i].Case = 1;
                circproplist[i].J = circuitInt1[i] == 0.
                    ? 0.
                    : 0.01 * (circproplist[i].Amps.re - circuitInt3[i]) / circuitInt1[i];
            } else {
                circproplist[i].Case = 0;
                circproplist[i].dV =
                    -0.01 * (circproplist[i].Amps.re - circuitInt3[i]) / circuitInt2[i];
            }
        } else {
            circproplist[i].Case = 0;
            circproplist[i].dV = circproplist[i].dVolts.re;
        }
    }

    // Magnetisation direction per element. A functional direction is evaluated
    // once from the element centroid; it depends only on position, so this is
    // equivalent to Static2D's per-iteration evaluation.
    std::vector<double> elementMagDir(elementCount, 0.0);
    for (std::size_t e = 0; e < elementCount; ++e) {
        const int label = elementLabel[e];
        double direction = labellist[label].MagDir;
        const std::string &function = labellist[label].MagDirFctn;
        if (!function.empty()) {
            CComplex centroid = 0;
            for (int j = 0; j < 3; ++j)
                centroid += CComplex(elementX[e][j], elementY[e][j]);
            centroid = centroid / units[LengthUnits] / 3.;
            char magbuff[4096];
            SNPRINTF(magbuff, sizeof magbuff,
                     "x=%.17g\ny=%.17g\nr=x\nz=y\ntheta=%.17g\nR=%.17g\nreturn %s",
                     centroid.re, centroid.im, arg(centroid) * 180 / PI, abs(centroid),
                     function.c_str());
            std::string source = magbuff;
            lua_State *lua = theLua->getLuaState();
            const int top1 = lua_gettop(lua);
            const int error =
                theLua->doString(source, femm::LuaInstance::LuaStackMode::Unsafe);
            if (error != 0)
                return 0;
            if (lua_gettop(lua) != top1) {
                source = lua_tostring(lua, -1);
                if (source.empty())
                    return 0;
                direction = Re(lua_tonumber(lua, -1));
                lua_pop(lua, 1);
            }
        }
        elementMagDir[e] = direction;
    }

    // Air-gap element contributions. Positions are applied to the quadrature
    // nodes without rebuilding the stored coupling topology.
    std::vector<femmsolver::CAirGapElement> ages = buildAirGapElements(view);
    if (ages.size() != view.airGaps().size())
        return 0;
    applyAirGapPositions(ages, airGapPositions);

    std::vector<double> mu1(elementCount, 1.0);
    std::vector<double> mu2(elementCount, 1.0);
    std::vector<double> v12(elementCount, 0.0);
    std::vector<double> vOld(nodeCount, 0.0);

    bool linearFlag = true;
    int iteration = 0;
    double relax = Relax;
    double residual = 0.0;
    double lastResidual = 0.0;

    do {
        if (iteration > 0)
            L.wipe();

        // Air-gap element contributions.
        for (std::size_t i = 0; i < ages.size(); ++i) {
            double MG[10][10];
            const double dt = (PI / 180.) *
                              (ages[i].totalArcLength / ages[i].totalArcElements);
            const double K = 2. * (ages[i].ro - ages[i].ri) /
                             (dt * (ages[i].ro + ages[i].ri));
            const double Ki = 1. / K;
            double ci = ages[i].InnerShift;
            double co = ages[i].OuterShift;
            if (ci > co) {
                ci = ci - co;
                co = 0;
            } else {
                ci = 1 - co + ci;
                co = 1;
            }

            MG[0][0] = (5 * nativePower(-1 + ci, 2) * nativePower(ci, 4) * (K + Ki)) / 48.;
            MG[0][1] = -((-1 + ci) * nativePower(ci, 3) *
                         (5 * (-1 + ci * (-5 + 4 * ci)) * K +
                          (-5 + ci * (-19 + 14 * ci)) * Ki)) / 48.;
            MG[0][2] = ((-1 + ci) * nativePower(ci, 2) *
                        (5 * (2 + ci * (-1 - 9 * ci + 6 * nativePower(ci, 2))) * K +
                         (10 + ci * (1 + 3 * ci * (-7 + 4 * ci))) * Ki)) / 48.;
            MG[0][3] = -(nativePower(-1 + ci, 2) * nativePower(ci, 2) *
                         (5 * (-2 + ci * (-3 + 4 * ci)) * K +
                          (2 + ci * (-3 + 2 * ci)) * Ki)) / 48.;
            MG[0][4] = (nativePower(-1 + ci, 3) * nativePower(ci, 3) * (5 * K - Ki)) / 48.;
            MG[0][5] = ((-1 + ci) * nativePower(ci, 2) * (-1 + co) * nativePower(co, 2) *
                        (K - 5 * Ki)) / 48.;
            MG[0][6] = -((-1 + ci) * nativePower(ci, 2) * co *
                         ((-1 + co * (-5 + 4 * co)) * K +
                          (5 + (19 - 14 * co) * co) * Ki)) / 48.;
            MG[0][7] = ((-1 + ci) * nativePower(ci, 2) *
                        ((2 + co * (-1 - 9 * co + 6 * nativePower(co, 2))) * K -
                         (10 + co * (1 + 3 * co * (-7 + 4 * co))) * Ki)) / 48.;
            MG[0][8] = -((-1 + ci) * nativePower(ci, 2) * (-1 + co) *
                         ((-2 + co * (-3 + 4 * co)) * K +
                          (-2 + (3 - 2 * co) * co) * Ki)) / 48.;
            MG[0][9] = ((-1 + ci) * nativePower(ci, 2) * nativePower(-1 + co, 2) * co *
                        (K + Ki)) / 48.;
            MG[1][1] = (nativePower(ci, 2) *
                        (5 * nativePower(1 + (5 - 4 * ci) * ci, 2) * K +
                         (5 + ci * (38 + ci * (49 + 4 * ci * (-29 + 11 * ci)))) * Ki)) / 48.;
            MG[1][2] = (-5 * ci * (-1 + 2 * ci) * (-2 + 3 * (-1 + ci) * ci) *
                            (-1 + ci * (-5 + 4 * ci)) * K +
                        ci * (10 + ci * (39 - ci * (50 + ci * (85 + 6 * ci * (-23 + 8 * ci))))) *
                            Ki) / 48.;
            MG[1][3] = ((-1 + ci) * ci *
                        (5 * (2 + ci * (13 + ci * (3 + 16 * (-2 + ci) * ci))) * K +
                         (-2 + 5 * ci * (1 + ci * (3 + 4 * (-2 + ci) * ci))) * Ki)) / 48.;
            MG[1][4] = -(nativePower(-1 + ci, 2) * nativePower(ci, 2) *
                         (5 * (-1 + ci * (-5 + 4 * ci)) * K + Ki + ci * (-1 + 2 * ci) * Ki)) / 48.;
            MG[1][5] = -(ci * (-1 + co) * nativePower(co, 2) *
                         ((-1 + ci * (-5 + 4 * ci)) * K + (5 + (19 - 14 * ci) * ci) * Ki)) / 48.;
            MG[1][6] = (ci * co *
                        ((-1 + ci * (-5 + 4 * ci)) * (-1 + co * (-5 + 4 * co)) * K +
                         (-5 + ci * (-19 + 14 * ci) - 19 * co + ci * (-77 + 58 * ci) * co +
                          2 * (7 + (29 - 22 * ci) * ci) * nativePower(co, 2)) *
                             Ki)) / 48.;
            MG[1][7] = (-(ci * (-1 + ci * (-5 + 4 * ci)) *
                          (2 + co * (-1 - 9 * co + 6 * nativePower(co, 2))) * K) +
                        ci * (-10 + co * (-1 + 3 * (7 - 4 * co) * co) +
                              ci * (-38 + co + 99 * nativePower(co, 2) -
                                    60 * nativePower(co, 3)) +
                              nativePower(ci, 2) *
                                  (28 + 2 * co * (-1 + 3 * co * (-13 + 8 * co)))) *
                            Ki) / 48.;
            MG[1][8] = (ci * (-1 + co) *
                        ((-1 + ci * (-5 + 4 * ci)) * (-2 + co * (-3 + 4 * co)) * K +
                         (2 + co * (-3 + 2 * co) +
                          nativePower(ci, 2) * (4 + 2 * (9 - 10 * co) * co) +
                          ci * (-2 + co * (-21 + 22 * co))) *
                             Ki)) / 48.;
            MG[1][9] = -(ci * nativePower(-1 + co, 2) * co *
                         ((-1 + ci * (-5 + 4 * ci)) * K +
                          (-1 + ci - 2 * nativePower(ci, 2)) * Ki)) / 48.;
            MG[2][2] = (5 * nativePower(-2 + ci + 9 * nativePower(ci, 2) - 6 * nativePower(ci, 3), 2) * K +
                        (20 + (-1 + ci) * ci *
                                  (-4 + 3 * (-1 + ci) * ci * (-25 + 24 * (-1 + ci) * ci))) *
                            Ki) / 48.;
            MG[2][3] = (-5 * (4 + nativePower(ci, 2) *
                                      (-33 + ci * (18 + ci * (65 + 6 * ci * (-13 + 4 * ci))))) *
                            K +
                        (4 + nativePower(ci, 2) *
                                 (39 - ci * (30 + ci * (115 + 6 * ci * (-25 + 8 * ci))))) *
                            Ki) / 48.;
            MG[2][4] = (nativePower(-1 + ci, 2) * ci *
                        (5 * (2 + ci * (-1 - 9 * ci + 6 * nativePower(ci, 2))) * K +
                         (-2 + ci * (-5 + 3 * ci * (-5 + 4 * ci))) * Ki)) / 48.;
            MG[2][5] = ((-1 + co) * nativePower(co, 2) *
                        ((2 + ci * (-1 - 9 * ci + 6 * nativePower(ci, 2))) * K -
                         (10 + ci * (1 + 3 * ci * (-7 + 4 * ci))) * Ki)) / 48.;
            MG[2][6] = (-((2 + ci * (-1 - 9 * ci + 6 * nativePower(ci, 2))) * co *
                          (-1 + co * (-5 + 4 * co)) * K) +
                        co * (-10 - 38 * co + 28 * nativePower(co, 2) +
                              nativePower(ci, 2) * (21 + 99 * co - 78 * nativePower(co, 2)) +
                              ci * (-1 + co - 2 * nativePower(co, 2)) +
                              12 * nativePower(ci, 3) * (-1 + co * (-5 + 4 * co))) *
                            Ki) / 48.;
            MG[2][7] = ((2 + ci * (-1 - 9 * ci + 6 * nativePower(ci, 2))) *
                            (2 + co * (-1 - 9 * co + 6 * nativePower(co, 2))) * K -
                        (2 * (10 + co) + 6 * nativePower(co, 2) * (-7 + 4 * co) +
                         3 * nativePower(ci, 2) * (-14 + co * (5 + (55 - 36 * co) * co)) +
                         ci * (2 + co * (5 + 3 * (5 - 4 * co) * co)) +
                         12 * nativePower(ci, 3) *
                             (2 + co * (-1 - 9 * co + 6 * nativePower(co, 2)))) *
                            Ki) / 48.;
            MG[2][8] = (-((2 + ci * (-1 - 9 * ci + 6 * nativePower(ci, 2))) *
                          (2 + co - 7 * nativePower(co, 2) + 4 * nativePower(co, 3)) * K) +
                        (-1 + co) *
                            (4 + 2 * ci * (5 + 3 * (5 - 4 * ci) * ci) +
                             3 * (-2 + ci * (3 + (17 - 12 * ci) * ci)) * co +
                             2 * (2 + ci * (-7 + 3 * ci * (-11 + 8 * ci))) *
                                 nativePower(co, 2)) *
                            Ki) / 48.;
            MG[2][9] = (nativePower(-1 + co, 2) * co *
                        ((2 + ci * (-1 - 9 * ci + 6 * nativePower(ci, 2))) * K +
                         (2 + ci * (5 + 3 * (5 - 4 * ci) * ci)) * Ki)) / 48.;
            MG[3][3] = (nativePower(-1 + ci, 2) *
                        (5 * nativePower(2 + (3 - 4 * ci) * ci, 2) * K +
                         (20 + ci * (36 + ci * (-35 - 60 * ci + 44 * nativePower(ci, 2)))) *
                             Ki)) / 48.;
            MG[3][4] = -(nativePower(-1 + ci, 3) * ci *
                         (5 * (-2 + ci * (-3 + 4 * ci)) * K +
                          (-10 + ci * (-9 + 14 * ci)) * Ki)) / 48.;
            MG[3][5] = -((-1 + ci) * (-1 + co) * nativePower(co, 2) *
                         ((-2 + ci * (-3 + 4 * ci)) * K +
                          (-2 + (3 - 2 * ci) * ci) * Ki)) / 48.;
            MG[3][6] = ((-1 + ci) * co *
                        ((-2 + ci * (-3 + 4 * ci)) * (-1 + co * (-5 + 4 * co)) * K +
                         (2 + ci * (-3 + 2 * ci) - 2 * co + ci * (-21 + 22 * ci) * co +
                          2 * (2 + (9 - 10 * ci) * ci) * nativePower(co, 2)) *
                             Ki)) / 48.;
            MG[3][7] = (-((2 + ci - 7 * nativePower(ci, 2) + 4 * nativePower(ci, 3)) *
                          (2 + co * (-1 - 9 * co + 6 * nativePower(co, 2))) * K) +
                        (-1 + ci) *
                            (4 + 2 * co * (5 + 3 * (5 - 4 * co) * co) +
                             ci * (-6 + 3 * co * (3 + (17 - 12 * co) * co)) +
                             2 * nativePower(ci, 2) *
                                 (2 + co * (-7 + 3 * co * (-11 + 8 * co)))) *
                            Ki) / 48.;
            MG[3][8] = ((-1 + ci) * (-1 + co) *
                        ((-2 + ci * (-3 + 4 * ci)) * (-2 + co * (-3 + 4 * co)) * K +
                         (-20 + 3 * ci * (1 + 2 * co) * (-6 + 5 * co) + 2 * co * (-9 + 14 * co) +
                          nativePower(ci, 2) * (28 + 30 * co - 44 * nativePower(co, 2))) *
                             Ki)) / 48.;
            MG[3][9] = -((-1 + ci) * nativePower(-1 + co, 2) * co *
                         ((-2 + ci * (-3 + 4 * ci)) * K + (10 + (9 - 14 * ci) * ci) * Ki)) / 48.;
            MG[4][4] = (5 * nativePower(-1 + ci, 4) * nativePower(ci, 2) * (K + Ki)) / 48.;
            MG[4][5] = (nativePower(-1 + ci, 2) * ci * (-1 + co) * nativePower(co, 2) *
                        (K + Ki)) / 48.;
            MG[4][6] = -(nativePower(-1 + ci, 2) * ci * co *
                         ((-1 + co * (-5 + 4 * co)) * K +
                          (-1 + co - 2 * nativePower(co, 2)) * Ki)) / 48.;
            MG[4][7] = (nativePower(-1 + ci, 2) * ci *
                        ((2 + co * (-1 - 9 * co + 6 * nativePower(co, 2))) * K +
                         (2 + co * (5 + 3 * (5 - 4 * co) * co)) * Ki)) / 48.;
            MG[4][8] = -(nativePower(-1 + ci, 2) * ci * (-1 + co) *
                         ((-2 + co * (-3 + 4 * co)) * K +
                          (10 + (9 - 14 * co) * co) * Ki)) / 48.;
            MG[4][9] = (nativePower(-1 + ci, 2) * ci * nativePower(-1 + co, 2) * co *
                        (K - 5 * Ki)) / 48.;
            MG[5][5] = (5 * nativePower(-1 + co, 2) * nativePower(co, 4) * (K + Ki)) / 48.;
            MG[5][6] = -((-1 + co) * nativePower(co, 3) *
                         (5 * (-1 + co * (-5 + 4 * co)) * K +
                          (-5 + co * (-19 + 14 * co)) * Ki)) / 48.;
            MG[5][7] = ((-1 + co) * nativePower(co, 2) *
                        (5 * (2 + co * (-1 - 9 * co + 6 * nativePower(co, 2))) * K +
                         (10 + co * (1 + 3 * co * (-7 + 4 * co))) * Ki)) / 48.;
            MG[5][8] = -(nativePower(-1 + co, 2) * nativePower(co, 2) *
                         (5 * (-2 + co * (-3 + 4 * co)) * K +
                          (2 + co * (-3 + 2 * co)) * Ki)) / 48.;
            MG[5][9] = (nativePower(-1 + co, 3) * nativePower(co, 3) * (5 * K - Ki)) / 48.;
            MG[6][6] = (nativePower(co, 2) *
                        (5 * nativePower(1 + (5 - 4 * co) * co, 2) * K +
                         (5 + co * (38 + co * (49 + 4 * co * (-29 + 11 * co)))) * Ki)) / 48.;
            MG[6][7] = (-5 * co * (-1 + 2 * co) * (-2 + 3 * (-1 + co) * co) *
                            (-1 + co * (-5 + 4 * co)) * K +
                        co * (10 + co * (39 - co * (50 + co * (85 + 6 * co * (-23 + 8 * co))))) *
                            Ki) / 48.;
            MG[6][8] = ((-1 + co) * co *
                        (5 * (2 + co * (13 + co * (3 + 16 * (-2 + co) * co))) * K +
                         (-2 + 5 * co * (1 + co * (3 + 4 * (-2 + co) * co))) * Ki)) / 48.;
            MG[6][9] = -(nativePower(-1 + co, 2) * nativePower(co, 2) *
                         (5 * (-1 + co * (-5 + 4 * co)) * K + Ki + co * (-1 + 2 * co) * Ki)) / 48.;
            MG[7][7] = (5 * nativePower(-2 + co + 9 * nativePower(co, 2) - 6 * nativePower(co, 3), 2) * K +
                        (20 + (-1 + co) * co *
                                  (-4 + 3 * (-1 + co) * co * (-25 + 24 * (-1 + co) * co))) *
                            Ki) / 48.;
            MG[7][8] = (-5 * (4 + nativePower(co, 2) *
                                      (-33 + co * (18 + co * (65 + 6 * co * (-13 + 4 * co))))) *
                            K +
                        (4 + nativePower(co, 2) *
                                 (39 - co * (30 + co * (115 + 6 * co * (-25 + 8 * co))))) *
                            Ki) / 48.;
            MG[7][9] = (nativePower(-1 + co, 2) * co *
                        (5 * (2 + co * (-1 - 9 * co + 6 * nativePower(co, 2))) * K +
                         (-2 + co * (-5 + 3 * co * (-5 + 4 * co))) * Ki)) / 48.;
            MG[8][8] = (nativePower(-1 + co, 2) *
                        (5 * nativePower(2 + (3 - 4 * co) * co, 2) * K +
                         (20 + co * (36 + co * (-35 - 60 * co + 44 * nativePower(co, 2)))) *
                             Ki)) / 48.;
            MG[8][9] = -(nativePower(-1 + co, 3) * co *
                         (5 * (-2 + co * (-3 + 4 * co)) * K +
                          (-10 + co * (-9 + 14 * co)) * Ki)) / 48.;
            MG[9][9] = (5 * nativePower(-1 + co, 4) * nativePower(co, 2) * (K + Ki)) / 48.;

            for (int k = 0; k < ages[i].totalArcElements; ++k) {
                int nn[10];
                double ww[10];
                const int arcElements = ages[i].totalArcElements;
                if ((k - 1) < 0) {
                    nn[0] = ages[i].quadNode[arcElements - 1].n0;
                    ww[0] = ages[i].quadNode[arcElements - 1].w0;
                } else {
                    nn[0] = ages[i].quadNode[k - 1].n0;
                    ww[0] = ages[i].quadNode[k - 1].w0;
                }
                nn[1] = ages[i].quadNode[k].n0;
                nn[2] = ages[i].quadNode[k].n1;
                nn[3] = ages[i].quadNode[k + 1].n1;
                ww[1] = ages[i].quadNode[k].w0;
                ww[2] = ages[i].quadNode[k].w1;
                ww[3] = ages[i].quadNode[k + 1].w1;
                if ((k + 2) > arcElements) {
                    nn[4] = ages[i].quadNode[1].n1;
                    ww[4] = ages[i].quadNode[1].w1;
                } else {
                    nn[4] = ages[i].quadNode[k + 2].n1;
                    ww[4] = ages[i].quadNode[k + 2].w1;
                }
                if ((k - 1) < 0) {
                    nn[5] = ages[i].quadNode[arcElements - 1].n2;
                    ww[5] = ages[i].quadNode[arcElements - 1].w2;
                } else {
                    nn[5] = ages[i].quadNode[k - 1].n2;
                    ww[5] = ages[i].quadNode[k - 1].w2;
                }
                nn[6] = ages[i].quadNode[k].n2;
                nn[7] = ages[i].quadNode[k].n3;
                nn[8] = ages[i].quadNode[k + 1].n3;
                ww[6] = ages[i].quadNode[k].w2;
                ww[7] = ages[i].quadNode[k].w3;
                ww[8] = ages[i].quadNode[k + 1].w3;
                if ((k + 2) > arcElements) {
                    nn[9] = ages[i].quadNode[1].n3;
                    ww[9] = ages[i].quadNode[1].w3;
                } else {
                    nn[9] = ages[i].quadNode[k + 2].n3;
                    ww[9] = ages[i].quadNode[k + 2].w3;
                }
                if ((k == 0) && (ages[i].BdryFormat == 1)) {
                    ww[0] = -ww[0];
                    ww[5] = -ww[5];
                }
                if (((k + 1) == arcElements) && (ages[i].BdryFormat == 1)) {
                    ww[4] = -ww[4];
                    ww[9] = -ww[9];
                }
                for (int ii = 0; ii < 10; ++ii)
                    nn[ii] = mapped(static_cast<std::size_t>(nn[ii]));
                for (int ii = 0; ii < 10; ++ii)
                    for (int jj = ii; jj < 10; ++jj)
                        L.add_to(MG[ii][jj] * ww[ii] * ww[jj], nn[ii], nn[jj]);
            }
        }

        // Volume element contributions.
        for (std::size_t e = 0; e < elementCount; ++e) {
            const std::array<int, 3> &n = elementNodes[e];
            const double *x = elementX[e].data();
            const double *y = elementY[e].data();
            const int label = elementLabel[e];
            const int block = elementBlock[e];

            double Me[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
            double be[3] = {0., 0., 0.};
            double Mx[3][3];
            double My[3][3];
            double Mxy[3][3];
            double Mn[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};

            double p[3];
            double q[3];
            double l[3];
            p[0] = y[1] - y[2];
            p[1] = y[2] - y[0];
            p[2] = y[0] - y[1];
            q[0] = x[2] - x[1];
            q[1] = x[0] - x[2];
            q[2] = x[1] - x[0];
            for (int j = 0, k = 1; j < 3; k++, j++) {
                if (k == 3)
                    k = 0;
                l[j] = std::sqrt(std::pow(x[k] - x[j], 2.) + std::pow(y[k] - y[j], 2.));
            }
            const double a = (p[0] * q[1] - p[1] * q[0]) / 2.;

            double K = (-1. / (4. * a));
            for (int j = 0; j < 3; ++j)
                for (int k = j; k < 3; ++k) {
                    Mx[j][k] = K * p[j] * p[k];
                    if (j != k)
                        Mx[k][j] = K * p[j] * p[k];
                }
            K = (-1. / (4. * a));
            for (int j = 0; j < 3; ++j)
                for (int k = j; k < 3; ++k) {
                    My[j][k] = K * q[j] * q[k];
                    if (j != k)
                        My[k][j] = K * q[j] * q[k];
                }
            K = (-1. / (4. * a));
            for (int j = 0; j < 3; ++j)
                for (int k = j; k < 3; ++k) {
                    Mxy[j][k] = K * (p[j] * q[k] + p[k] * q[j]);
                    if (j != k)
                        Mxy[k][j] = K * (p[j] * q[k] + p[k] * q[j]);
                }

            // Derivative boundary conditions.
            for (int j = 0; j < 3; ++j) {
                if (elementSides[e][j] >= 0 &&
                    lineproplist[elementSides[e][j]].BdryFormat == 2) {
                    const double coefficient =
                        -0.0001 * c * lineproplist[elementSides[e][j]].c0.re * l[j] / 6.;
                    const int k = (j + 1 == 3) ? 0 : j + 1;
                    Me[j][j] += coefficient * 2.;
                    Me[k][k] += coefficient * 2.;
                    Me[j][k] += coefficient;
                    Me[k][j] += coefficient;
                    const double source =
                        (lineproplist[elementSides[e][j]].c1.re * l[j] / 2.) * 0.0001;
                    be[j] += source;
                    be[k] += source;
                }
            }

            // Current density in the block.
            for (int j = 0; j < 3; ++j) {
                double t = 0;
                if (labellist[label].InCircuit >= 0) {
                    const int circuit = labellist[label].InCircuit;
                    if (circproplist[circuit].Case == 1)
                        t = circproplist[circuit].J.Re();
                    if (circproplist[circuit].Case == 0)
                        t = -circproplist[circuit].dV.Re() * blockproplist[block].Cduct;
                }
                K = -(blockproplist[block].J.re + t) * a / 3.;
                be[j] += K;
            }

            // Constant or functional magnetisation in the block.
            const double magDir = elementMagDir[e];
            for (int j = 0; j < 3; ++j) {
                const int k = (j + 1 == 3) ? 0 : j + 1;
                K = 0.0001 * blockproplist[block].H_c *
                    (std::cos(magDir * PI / 180.) * (x[k] - x[j]) +
                     std::sin(magDir * PI / 180.) * (y[k] - y[j])) / 2.;
                be[j] += K;
                be[k] += K;
            }

            // Permeability update. Iteration zero initialises mu1/mu2 from the
            // block's lamination and flags a nonlinear solve; later iterations
            // update mu1/mu2 from the B-H curve and build the Newton matrix Mn.
            if (iteration == 0) {
                if (blockproplist[block].LamType == 0) {
                    const double t = blockproplist[block].LamFill;
                    mu1[e] = blockproplist[block].mu_x * t + (1. - t);
                    mu2[e] = blockproplist[block].mu_y * t + (1. - t);
                } else if (blockproplist[block].LamType == 1) {
                    const double t = blockproplist[block].LamFill;
                    const double mu = blockproplist[block].mu_x;
                    mu1[e] = mu * t + (1. - t);
                    mu2[e] = mu / (t + mu * (1. - t));
                } else if (blockproplist[block].LamType == 2) {
                    const double t = blockproplist[block].LamFill;
                    const double mu = blockproplist[block].mu_y;
                    mu2[e] = mu * t + (1. - t);
                    mu1[e] = mu / (t + mu * (1. - t));
                } else {
                    mu1[e] = 1.0;
                    mu2[e] = 1.0;
                }
                if (blockproplist[block].BHpoints != 0)
                    linearFlag = false;
            } else {
                double B1 = 0.0;
                double B2 = 0.0;
                double mu = 0.0;
                double dv = 0.0;
                if (blockproplist[block].LamType == 0 &&
                    mu1[e] == mu2[e] && blockproplist[block].BHpoints > 0) {
                    for (int j = 0; j < 3; ++j) {
                        B1 += L.solution()[n[j]] * q[j];
                        B2 += L.solution()[n[j]] * p[j];
                    }
                    const double B = c * std::sqrt(B1 * B1 + B2 * B2) / (0.02 * a);
                    blockproplist[block].GetBHProps(B, mu, dv);
                    mu = 1. / (muo * mu);
                    mu1[e] = mu;
                    mu2[e] = mu;
                    double v[3];
                    for (int j = 0; j < 3; ++j) {
                        v[j] = 0.0;
                        for (int w = 0; w < 3; ++w)
                            v[j] += (Mx[j][w] + My[j][w]) * L.solution()[n[w]];
                    }
                    const double coefficient = -200. * c * c * c * dv / a;
                    for (int j = 0; j < 3; ++j)
                        for (int w = 0; w < 3; ++w)
                            Mn[j][w] = coefficient * v[j] * v[w];
                }
                if (blockproplist[block].LamType == 1 && blockproplist[block].BHpoints > 0) {
                    const double t = blockproplist[block].LamFill;
                    for (int j = 0; j < 3; ++j) {
                        B1 += L.solution()[n[j]] * q[j];
                        B2 += L.solution()[n[j]] * p[j] / t;
                    }
                    const double B = c * std::sqrt(B1 * B1 + B2 * B2) / (0.02 * a);
                    blockproplist[block].GetBHProps(B, mu, dv);
                    mu = 1. / (muo * mu);
                    mu1[e] = mu * t;
                    mu2[e] = mu / (t + mu * (1. - t));
                    double v[3];
                    double u[3];
                    for (int j = 0; j < 3; ++j) {
                        v[j] = 0.0;
                        u[j] = 0.0;
                        for (int w = 0; w < 3; ++w) {
                            v[j] += (My[j][w] / t + Mx[j][w]) * L.solution()[n[w]];
                            u[j] += (My[j][w] / t + t * Mx[j][w]) * L.solution()[n[w]];
                        }
                    }
                    const double coefficient = -100. * c * c * c * dv / a;
                    for (int j = 0; j < 3; ++j)
                        for (int w = 0; w < 3; ++w)
                            Mn[j][w] = coefficient * (v[j] * u[w] + v[w] * u[j]);
                }
                if (blockproplist[block].LamType == 2 && blockproplist[block].BHpoints > 0) {
                    const double t = blockproplist[block].LamFill;
                    for (int j = 0; j < 3; ++j) {
                        B1 += (L.solution()[n[j]] * q[j]) / t;
                        B2 += L.solution()[n[j]] * p[j];
                    }
                    const double B = c * std::sqrt(B1 * B1 + B2 * B2) / (0.02 * a);
                    blockproplist[block].GetBHProps(B, mu, dv);
                    mu = 1. / (muo * mu);
                    mu2[e] = mu * t;
                    mu1[e] = mu / (t + mu * (1. - t));
                    double v[3];
                    double u[3];
                    for (int j = 0; j < 3; ++j) {
                        v[j] = 0.0;
                        u[j] = 0.0;
                        for (int w = 0; w < 3; ++w) {
                            v[j] += (Mx[j][w] / t + My[j][w]) * L.solution()[n[w]];
                            u[j] += (Mx[j][w] / t + t * My[j][w]) * L.solution()[n[w]];
                        }
                    }
                    const double coefficient = -100. * c * c * c * dv / a;
                    for (int j = 0; j < 3; ++j)
                        for (int w = 0; w < 3; ++w)
                            Mn[j][w] = coefficient * (v[j] * u[w] + v[w] * u[j]);
                }
            }

            for (int j = 0; j < 3; ++j)
                for (int k = 0; k < 3; ++k) {
                    Me[j][k] += Mx[j][k] / mu2[e] + My[j][k] / mu1[e] +
                                Mxy[j][k] * v12[e] + Mn[j][k];
                    be[j] += Mn[j][k] * L.solution()[n[k]];
                }

            for (int j = 0; j < 3; ++j) {
                for (int k = j; k < 3; ++k)
                    L.add_to(-Me[j][k], n[j], n[k]);
                L.rhs()[n[j]] -= be[j];
            }
        }

        // Point currents and fixed point values.
        for (std::size_t i = 0; i < nodeCount; ++i) {
            const std::int32_t raw = view.nodeBoundaryMarker(static_cast<MeshIndex>(i));
            const int marker = raw > 1 ? raw - 2 : -1;
            if (marker < 0 || marker >= NumPointProps)
                continue;
            L.rhs()[mapped(i)] += 0.01 * nodeproplist[marker].J.re;
            if (nodeproplist[marker].J.re == 0 && nodeproplist[marker].J.im == 0)
                L.set_value(mapped(i), nodeproplist[marker].A.re / c);
        }

        // Fixed boundary conditions along segments.
        for (std::size_t e = 0; e < elementCount; ++e) {
            const std::array<int, 3> &nodes = elementNodes[e];
            for (int j = 0; j < 3; ++j) {
                const int side = elementSides[e][j];
                if (side < 0 || lineproplist[side].BdryFormat != 0)
                    continue;
                const int k = (j + 1 == 3) ? 0 : j + 1;
                for (int endpoint = 0; endpoint < 2; ++endpoint) {
                    const int node = nodes[endpoint == 0 ? j : k];
                    const double x = elementX[e][endpoint == 0 ? j : k];
                    const double y = elementY[e][endpoint == 0 ? j : k];
                    double value;
                    if (Coords == 0) {
                        const double scaledX = x / units[LengthUnits];
                        const double scaledY = y / units[LengthUnits];
                        value = lineproplist[side].A0 + scaledX * lineproplist[side].A1 +
                                scaledY * lineproplist[side].A2;
                    } else {
                        const double radius = std::sqrt(x * x + y * y);
                        const double theta =
                            (x == 0 && y == 0) ? 0.0 : std::atan2(y, x) / DEG;
                        value = lineproplist[side].A0 +
                                (radius / units[LengthUnits]) * lineproplist[side].A1 +
                                theta * lineproplist[side].A2;
                    }
                    value *= std::cos(lineproplist[side].phi * DEG);
                    L.set_value(node, value / c);
                }
            }
        }

        // Ordinary periodic and antiperiodic constraints.
        for (const auto &constraint : view.periodicConstraints()) {
            L.constrain_periodic(mapped(constraint.first), mapped(constraint.second),
                                 constraint.periodicity ==
                                     femm::mesh::SolverMesh::Periodicity::Antiperiodic);
        }

        for (std::size_t j = 0; j < nodeCount; ++j)
            vOld[j] = L.solution()[j];

        femm::SolveOptions options;
        options.warm_start = iteration > 0;
        if (L.solve(options).converged == false)
            return 0;

        if (!linearFlag) {
            double change = 0.0;
            double magnitude = 0.0;
            for (std::size_t j = 0; j < nodeCount; ++j) {
                change += (L.solution()[j] - vOld[j]) * (L.solution()[j] - vOld[j]);
                magnitude += L.solution()[j] * L.solution()[j];
            }
            if (magnitude == 0) {
                linearFlag = true;
            } else {
                lastResidual = residual;
                residual = std::sqrt(change / magnitude);
            }
            if (iteration > 5) {
                if (residual > lastResidual && relax > 0.125)
                    relax /= 2.;
                else
                    relax += 0.1 * (1. - relax);
                for (std::size_t j = 0; j < nodeCount; ++j)
                    L.solution()[j] = relax * L.solution()[j] + (1.0 - relax) * vOld[j];
            }
        }

        if (residual < 100. * Precision && iteration > 0)
            linearFlag = true;

        ++iteration;
    } while (linearFlag == false);

    for (std::size_t i = 0; i < nodeCount; ++i)
        L.rhs()[i] = L.solution()[i] * c;

    return 1;
}
