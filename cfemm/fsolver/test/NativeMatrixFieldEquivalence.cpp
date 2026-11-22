#include "NativeMatrixFieldEquivalence.h"

#include "CBlockLabel.h"
#include "CMaterialProp.h"
#include "mesh/LogicalMeshView.h"
#include "mesh/SolverMesh.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace xfemm_test {
namespace {

using Snapshot = femm::FSolverAnalysisBackend::StaticAssemblySnapshot;
using IndexPair = std::pair<std::size_t, std::size_t>;
using Matrix = std::map<IndexPair, double>;

MatrixFieldResult fail(const std::string &reason)
{
    MatrixFieldResult result;
    result.message = reason;
    return result;
}

double norm2(const std::vector<double> &v)
{
    double sum = 0;
    for (double x : v)
        sum += x * x;
    return std::sqrt(sum);
}

double relativeResidual(const Snapshot &snapshot, const std::vector<double> &solution)
{
    std::vector<double> ax(solution.size(), 0);
    for (const auto &entry : snapshot.upperTriangle) {
        ax[entry.row] += entry.value * solution[entry.column];
        if (entry.row != entry.column)
            ax[entry.column] += entry.value * solution[entry.row];
    }
    std::vector<double> residual(ax.size());
    for (std::size_t i = 0; i < residual.size(); ++i)
        residual[i] = ax[i] - snapshot.rhs[i];
    return norm2(residual) / std::max({norm2(ax), norm2(snapshot.rhs), 1e-30});
}

Matrix normalizedMatrix(const Snapshot &snapshot,
                        const std::vector<std::size_t> &localToReference,
                        double threshold)
{
    Matrix result;
    for (const auto &entry : snapshot.upperTriangle) {
        if (!std::isfinite(entry.value))
            throw std::runtime_error("non-finite matrix coefficient");
        if (entry.row < 0 || entry.column < entry.row ||
            entry.column >= snapshot.dimension)
            throw std::runtime_error("invalid upper-triangle entry");
        if (std::abs(entry.value) <= threshold)
            continue;
        const auto indices = std::minmax(
            localToReference[static_cast<std::size_t>(entry.row)],
            localToReference[static_cast<std::size_t>(entry.column)]);
        const IndexPair key(indices.first, indices.second);
        if (!result.emplace(key, entry.value).second)
            throw std::runtime_error("duplicate mapped matrix entry");
    }
    return result;
}

bool validSnapshot(const Snapshot &snapshot, const std::vector<double> &x)
{
    if (snapshot.dimension <= 0 || snapshot.capturedIterations == 0 ||
        snapshot.rhs.size() != static_cast<std::size_t>(snapshot.dimension) ||
        x.size() != snapshot.rhs.size())
        return false;
    for (double v : snapshot.rhs)
        if (!std::isfinite(v)) return false;
    for (double v : x)
        if (!std::isfinite(v)) return false;
    for (const auto &entry : snapshot.upperTriangle)
        if (entry.row < 0 || entry.column < entry.row ||
            entry.column >= snapshot.dimension || !std::isfinite(entry.value))
            return false;
    return true;
}

std::vector<double> solveVector(const femm::FSolverAnalysisBackend &solver)
{
    const auto &x = solver.solvedSystem().solution();
    return {x.begin(), x.end()};
}

struct ElementField { double bx = 0, by = 0, hx = 0, hy = 0; };

ElementField elementField(const femm::mesh::SolverMesh &mesh,
                          const femm::mesh::SolverMesh::Element &element,
                          const std::vector<double> &a,
                          femm::CMMaterialProp &material,
                          double magnetAngleDegrees)
{
    const auto &p0 = mesh.nodes[element.nodes[0]];
    const auto &p1 = mesh.nodes[element.nodes[1]];
    const auto &p2 = mesh.nodes[element.nodes[2]];
    const double dx1 = p1.x - p0.x;
    const double dy1 = p1.y - p0.y;
    const double dx2 = p2.x - p0.x;
    const double dy2 = p2.y - p0.y;
    const double area2 = dx1 * dy2 - dx2 * dy1;
    if (!std::isfinite(area2) || std::abs(area2) < 1e-24)
        throw std::runtime_error("degenerate element in field comparison");

    const double dA1 = a[element.nodes[1]] - a[element.nodes[0]];
    const double dA2 = a[element.nodes[2]] - a[element.nodes[0]];
    const double dAdx = (dA1 * dy2 - dA2 * dy1) / area2;
    const double dAdy = (dx1 * dA2 - dx2 * dA1) / area2;

    ElementField field;
    field.bx = dAdy;
    field.by = -dAdx;
    double mux = 0, muy = 0;
    material.GetMu(field.bx, field.by, mux, muy);
    const double mu0 = 4.0 * std::acos(-1.0) * 1e-7;
    if (!std::isfinite(mux) || !std::isfinite(muy) || mux <= 0 || muy <= 0)
        throw std::runtime_error("invalid constitutive permeability");
    const double angle = magnetAngleDegrees * std::acos(-1.0) / 180.;
    field.hx = field.bx / (mu0 * mux) - material.H_c * std::cos(angle);
    field.hy = field.by / (mu0 * muy) - material.H_c * std::sin(angle);
    if (!std::isfinite(field.bx) || !std::isfinite(field.by) ||
        !std::isfinite(field.hx) || !std::isfinite(field.hy))
        throw std::runtime_error("non-finite B/H element field");
    return field;
}

std::vector<double> canonicalPotential(const femm::mesh::SolverMesh &mesh,
                                       const femm::TrialSolution &trial,
                                       double coordinateTolerance)
{
    femm::RealNodalSolution canonical;
    canonical.x.reserve(mesh.nodes.size());
    canonical.y.reserve(mesh.nodes.size());
    canonical.magneticVectorPotential.assign(mesh.nodes.size(), 0);
    for (const auto &node : mesh.nodes) {
        canonical.x.push_back(node.x);
        canonical.y.push_back(node.y);
    }
    std::vector<std::size_t> toCanonical;
    std::string message;
    if (!associateNodalCoordinates(canonical, trial.real->nodal,
                                   coordinateTolerance, toCanonical, message))
        throw std::runtime_error("cannot map field to canonical mesh: " + message);

    std::vector<double> a(mesh.nodes.size());
    for (std::size_t i = 0; i < toCanonical.size(); ++i)
        a[toCanonical[i]] = trial.real->nodal.magneticVectorPotential[i];
    return a;
}

} // namespace

MatrixComparisonResult compareAssemblySnapshots(
    const Snapshot &reference, const Snapshot &native,
    const std::vector<std::size_t> &nativeToReference,
    const std::vector<double> &referenceSolution,
    const std::vector<double> &nativeSolution,
    const MatrixFieldOptions &options)
{
    MatrixComparisonResult result;
    if (!validSnapshot(reference, referenceSolution) ||
        !validSnapshot(native, nativeSolution) ||
        reference.dimension != native.dimension ||
        nativeToReference.size() != static_cast<std::size_t>(native.dimension)) {
        result.message = "missing/invalid or unequal assembled systems";
        return result;
    }
    const std::size_t n = nativeToReference.size();
    std::vector<bool> assigned(n, false);
    for (std::size_t p : nativeToReference) {
        if (p >= n || assigned[p]) {
            result.message = "invalid or duplicate native-to-reference DOF mapping";
            return result;
        }
        assigned[p] = true;
    }
    if (!std::isfinite(options.sparsityThreshold) || options.sparsityThreshold < 0 ||
        !std::isfinite(options.maximumRelativeResidual) ||
        options.maximumRelativeResidual < 0) {
        result.message = "invalid matrix or residual tolerances";
        return result;
    }

    try {
        const std::vector<std::size_t> identity = [&] {
            std::vector<std::size_t> ids(n);
            for (std::size_t i = 0; i < n; ++i) ids[i] = i;
            return ids;
        }();
        const Matrix a = normalizedMatrix(reference, identity, options.sparsityThreshold);
        const Matrix b = normalizedMatrix(native, nativeToReference, options.sparsityThreshold);
        result.numericNonzeros = a.size();
        result.referenceRelativeResidual = relativeResidual(reference, referenceSolution);
        result.nativeRelativeResidual = relativeResidual(native, nativeSolution);
        if (!std::isfinite(result.referenceRelativeResidual) ||
            !std::isfinite(result.nativeRelativeResidual) ||
            result.referenceRelativeResidual > options.maximumRelativeResidual ||
            result.nativeRelativeResidual > options.maximumRelativeResidual) {
            std::ostringstream s;
            s << "equation residual exceeds " << options.maximumRelativeResidual
              << " (materialized=" << result.referenceRelativeResidual
              << ", native=" << result.nativeRelativeResidual << ")";
            result.message = s.str();
            return result;
        }

        if (options.compareMatricesAndRhs) {
            for (std::size_t i = 0; i < n; ++i) {
                const double expected = reference.rhs[nativeToReference[i]];
                const double actual = native.rhs[i];
                const double difference = std::abs(expected - actual);
                result.maximumRhsDifference =
                    std::max(result.maximumRhsDifference, difference);
                if (difference > options.rhsAbsoluteTolerance +
                    options.rhsRelativeTolerance * std::max(std::abs(expected), std::abs(actual))) {
                    result.message = "RHS differs at native DOF " + std::to_string(i);
                    return result;
                }
            }
            // Union of the sparse patterns detects both missing and extra couplings,
            // including those introduced by periodic constraints and AGE elements.
            std::set<IndexPair> keys;
            for (const auto &e : a) keys.insert(e.first);
            for (const auto &e : b) keys.insert(e.first);
            for (const auto &key : keys) {
                const auto ai = a.find(key), bi = b.find(key);
                if (ai == a.end() || bi == b.end()) {
                    ++result.sparsityMismatches;
                    continue;
                }
                const double difference = std::abs(ai->second - bi->second);
                result.maximumMatrixDifference =
                    std::max(result.maximumMatrixDifference, difference);
                if (difference > options.matrixAbsoluteTolerance +
                    options.matrixRelativeTolerance *
                        std::max(std::abs(ai->second), std::abs(bi->second))) {
                    result.message = "matrix coefficient differs at (" +
                                     std::to_string(key.first) + "," +
                                     std::to_string(key.second) + ")";
                    return result;
                }
            }
            if (result.sparsityMismatches) {
                result.message = "numerical matrix sparsity differs in " +
                                 std::to_string(result.sparsityMismatches) + " entries";
                return result;
            }
        }
    } catch (const std::exception &e) {
        result.message = e.what();
        return result;
    }

    result.passed = true;
    std::ostringstream s;
    s << "PASS nnz=" << result.numericNonzeros
      << " max_matrix_error=" << result.maximumMatrixDifference
      << " max_rhs_error=" << result.maximumRhsDifference
      << " residual_ref=" << result.referenceRelativeResidual
      << " residual_native=" << result.nativeRelativeResidual;
    result.message = s.str();
    return result;
}

MatrixFieldResult compareMatrixFieldSessions(
    femm::AnalysisSession &referenceSession,
    femm::FSolverAnalysisBackend &referenceSolver,
    femm::AnalysisSession &nativeSession,
    femm::FSolverAnalysisBackend &nativeSolver,
    const MatrixFieldOptions &options)
{
    if (referenceSession.nativeInstanced() || !nativeSession.nativeInstanced() ||
        !referenceSession.instancedMesh() || !nativeSession.instancedMesh())
        return fail("expected independent materialized and native instanced sessions");

    try {
        referenceSolver.setStaticAssemblyCapture(true);
        nativeSolver.setStaticAssemblyCapture(true);
        const femm::TrialSolution reference = referenceSession.solve();
        const femm::TrialSolution native = nativeSession.solve();
        if (!reference.real || !native.real)
            return fail("expected two planar static real solutions");

        MatrixFieldResult result;
        result.nodal = compareNodalSolutions(reference, native, options.nodal);
        if (!result.nodal.passed)
            return fail("nodal A mismatch: " + result.nodal.message);

        const auto mesh = referenceSession.mesh();
        const auto view = nativeSession.logicalView();
        if (!mesh || !view || !view->valid() || nativeSession.materializationCount())
            return fail("invalid or unexpectedly materialized logical mesh");
        if (mesh->nodes.size() != view->nodeCount() ||
            mesh->elements.size() != view->elementCount() ||
            mesh->airGaps.size() != view->airGaps().size() ||
            mesh->periodicConstraints.size() != view->periodicConstraints().size())
            return fail("canonical topology counts differ");
        for (std::size_t i = 0; i < mesh->nodes.size(); ++i) {
            double x = 0, y = 0;
            view->nodeCoordinates(i, x, y);
            if (std::hypot(mesh->nodes[i].x - x, mesh->nodes[i].y - y) >
                options.nodal.coordinateTolerance)
                return fail("logical/materialized coordinate mismatch at " + std::to_string(i));
        }
        for (std::size_t i = 0; i < mesh->elements.size(); ++i) {
            if (mesh->elements[i].nodes != view->elementNodes(i))
                return fail("logical/materialized element connectivity mismatch at " +
                            std::to_string(i));
        }

        std::vector<std::size_t> nativeToReference;
        std::string mappingError;
        if (!associateNodalCoordinates(reference.real->nodal, native.real->nodal,
                                      options.nodal.coordinateTolerance,
                                      nativeToReference, mappingError))
            return fail("cannot map matrix DOFs: " + mappingError);

        result.matrix = compareAssemblySnapshots(
            referenceSolver.staticAssemblySnapshot(), nativeSolver.staticAssemblySnapshot(),
            nativeToReference, solveVector(referenceSolver), solveVector(nativeSolver), options);
        if (!result.matrix.passed)
            return fail("matrix/residual mismatch: " + result.matrix.message);

        if (options.compareFields) {
            const auto referenceA = canonicalPotential(*mesh, reference,
                                                       options.nodal.coordinateTolerance);
            const auto nativeA = canonicalPotential(*mesh, native,
                                                    options.nodal.coordinateTolerance);
            const auto &labels = referenceSession.prepared().labels;
            const auto &nativeLabels = nativeSession.prepared().labels;
            if (labels.size() != nativeLabels.size())
                return fail("prepared label counts differ");

            std::vector<femm::CMMaterialProp> materials;
            for (const auto &item : referenceSession.model().problem().blockproplist) {
                const auto *material = dynamic_cast<const femm::CMMaterialProp *>(item.get());
                if (!material)
                    return fail("non-magnetic material in G5 field test");
                materials.push_back(*material);
                if (materials.back().BHpoints > 0)
                    materials.back().GetSlopes(0.0);
            }

            for (std::size_t e = 0; e < mesh->elements.size(); ++e) {
                const auto &element = mesh->elements[e];
                const int label = element.regionAttribute - 1;
                if (label < 0 || static_cast<std::size_t>(label) >= labels.size())
                    return fail("element references an invalid prepared region label");
                const auto &region = labels[label];
                const auto &nativeRegion = nativeLabels[label];
                if (region.BlockType != nativeRegion.BlockType ||
                    region.MagDir != nativeRegion.MagDir ||
                    region.MagDirFctn != nativeRegion.MagDirFctn)
                    return fail("native and materialized region physics differ");
                const int block = region.BlockType;
                if (block < 0 || static_cast<std::size_t>(block) >= materials.size())
                    return fail("element references invalid material");
                if (options.compareH && !region.MagDirFctn.empty())
                    return fail("functional magnetization H comparison requires a field callback");

                auto &material = materials[block];
                const auto r = elementField(*mesh, element, referenceA, material, region.MagDir);
                const auto n = elementField(*mesh, element, nativeA, material, region.MagDir);
                result.referenceBScale = std::max(result.referenceBScale, std::hypot(r.bx, r.by));
                result.referenceHScale = std::max(result.referenceHScale, std::hypot(r.hx, r.hy));
                result.maximumBError = std::max(result.maximumBError,
                    std::hypot(r.bx - n.bx, r.by - n.by));
                result.maximumHError = std::max(result.maximumHError,
                    std::hypot(r.hx - n.hx, r.hy - n.hy));
                ++result.comparedElements;
            }
            if (result.comparedElements != mesh->elements.size())
                return fail("field comparison omitted elements");
            if (result.maximumBError > options.fieldAbsoluteTolerance +
                                        options.fieldRelativeTolerance * result.referenceBScale)
                return fail("element B differs: max error=" +
                            std::to_string(result.maximumBError));
            if (options.compareH &&
                result.maximumHError > options.hAbsoluteTolerance +
                                       options.hRelativeTolerance * result.referenceHScale)
                return fail("element H differs: max error=" +
                            std::to_string(result.maximumHError));
        }

        result.passed = true;
        std::ostringstream s;
        s << "PASS " << result.nodal.message << " | " << result.matrix.message
          << " | elements=" << result.comparedElements
          << " max_B_T=" << result.maximumBError
          << " max_H_A_per_m=" << result.maximumHError;
        result.message = s.str();
        return result;
    } catch (const std::exception &e) {
        return fail(std::string("matrix/field equivalence failed: ") + e.what());
    }
}

} // namespace xfemm_test
