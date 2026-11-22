#include "NativeEquivalence.h"

#include "mesh/LogicalMeshView.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace xfemm_test {
namespace {

using Bucket = std::pair<long long, long long>;

EquivalenceResult failure(const std::string &message)
{
    EquivalenceResult result;
    result.message = message;
    return result;
}

bool finite(double value) { return std::isfinite(value); }

bool validCoordinates(const femm::RealNodalSolution &nodal, double cell,
                      std::string &message)
{
    const std::size_t count = nodal.magneticVectorPotential.size();
    if (nodal.x.size() != count || nodal.y.size() != count) {
        message = "coordinate and solution vector lengths differ";
        return false;
    }
    const double safeMaximum = static_cast<double>(std::numeric_limits<long long>::max() / 4);
    for (std::size_t i = 0; i < count; ++i) {
        if (!finite(nodal.x[i]) || !finite(nodal.y[i]) ||
            !finite(nodal.magneticVectorPotential[i])) {
            message = "non-finite coordinate or nodal solution at index " + std::to_string(i);
            return false;
        }
        if (std::abs(nodal.x[i] / cell) > safeMaximum ||
            std::abs(nodal.y[i] / cell) > safeMaximum) {
            message = "coordinate exceeds the supported association range";
            return false;
        }
    }
    return true;
}

Bucket bucketFor(double x, double y, double cell)
{
    return {std::llround(x / cell), std::llround(y / cell)};
}

/** Find a unique reference node for each trial node. Search adjacent buckets
 * because coordinates close in Euclidean distance can straddle a rounding
 * boundary. Reject ambiguities rather than silently selecting the last node.
 */
bool nodeAssociation(const femm::RealNodalSolution &reference,
                     const femm::RealNodalSolution &trial, double tolerance,
                     std::vector<std::size_t> &trialToReference,
                     std::string &message)
{
    std::map<Bucket, std::vector<std::size_t>> buckets;
    for (std::size_t i = 0; i < reference.x.size(); ++i)
        buckets[bucketFor(reference.x[i], reference.y[i], tolerance)].push_back(i);

    std::vector<bool> used(reference.x.size(), false);
    trialToReference.resize(trial.x.size());
    for (std::size_t node = 0; node < trial.x.size(); ++node) {
        const Bucket cell = bucketFor(trial.x[node], trial.y[node], tolerance);
        std::size_t selected = reference.x.size();
        std::size_t candidates = 0;
        for (int dx = -1; dx <= 1; ++dx) {
            for (int dy = -1; dy <= 1; ++dy) {
                const auto found = buckets.find({cell.first + dx, cell.second + dy});
                if (found == buckets.end())
                    continue;
                for (const std::size_t candidate : found->second) {
                    const double xerr = reference.x[candidate] - trial.x[node];
                    const double yerr = reference.y[candidate] - trial.y[node];
                    if (std::hypot(xerr, yerr) <= tolerance) {
                        selected = candidate;
                        ++candidates;
                    }
                }
            }
        }
        if (candidates != 1) {
            message = "native node " + std::to_string(node) +
                " has " + std::to_string(candidates) +
                " geometric matches; explicit provenance is required for coincident DOFs";
            return false;
        }
        if (used[selected]) {
            message = "reference node " + std::to_string(selected) +
                " was matched more than once";
            return false;
        }
        used[selected] = true;
        trialToReference[node] = selected;
    }
    return true;
}

bool validOptions(const EquivalenceOptions &options)
{
    return finite(options.coordinateTolerance) && options.coordinateTolerance > 0 &&
           finite(options.absoluteTolerance) && options.absoluteTolerance >= 0 &&
           finite(options.relativeTolerance) && options.relativeTolerance >= 0 &&
           finite(options.minimumFieldMagnitude) && options.minimumFieldMagnitude >= 0;
}

} // namespace

EquivalenceResult compareNodalSolutions(const femm::TrialSolution &materialized,
                                         const femm::TrialSolution &native,
                                         const EquivalenceOptions &options)
{
    if (!validOptions(options))
        return failure("invalid equivalence tolerances");
    if (!materialized.real || !native.real)
        return failure("both solutions must contain real nodal values");

    const auto &reference = materialized.real->nodal;
    const auto &trial = native.real->nodal;
    if (reference.magneticVectorPotential.empty())
        return failure("reference nodal field is empty");
    if (reference.magneticVectorPotential.size() != trial.magneticVectorPotential.size())
        return failure("different reference and native node counts");

    std::string message;
    if (!validCoordinates(reference, options.coordinateTolerance, message))
        return failure("materialized: " + message);
    if (!validCoordinates(trial, options.coordinateTolerance, message))
        return failure("native: " + message);

    std::vector<std::size_t> association;
    if (!nodeAssociation(reference, trial, options.coordinateTolerance, association, message))
        return failure(message);

    EquivalenceResult result;
    result.matchedNodes = association.size();
    for (double a : reference.magneticVectorPotential)
        result.referenceScale = std::max(result.referenceScale, std::abs(a));
    if (options.requireNontrivialField &&
        result.referenceScale <= options.minimumFieldMagnitude)
        return failure("reference field is trivial; equivalence would not test excitation");

    if (options.removeGaugeOffset) {
        // Constant offset minimising the squared error; only enable for
        // fixtures without a prescribed vector-potential reference.
        double delta = 0;
        for (std::size_t i = 0; i < association.size(); ++i)
            delta += reference.magneticVectorPotential[association[i]] -
                     trial.magneticVectorPotential[i];
        result.appliedGaugeOffset = delta / static_cast<double>(association.size());
    }

    double sumSquared = 0;
    for (std::size_t i = 0; i < association.size(); ++i) {
        const std::size_t referenceNode = association[i];
        const double difference = (trial.magneticVectorPotential[i] +
                                   result.appliedGaugeOffset) -
                                  reference.magneticVectorPotential[referenceNode];
        const double error = std::abs(difference);
        sumSquared += difference * difference;
        if (error > result.maxAbsoluteError) {
            result.maxAbsoluteError = error;
            result.worstReferenceNode = referenceNode;
        }
    }
    result.rmsError = std::sqrt(sumSquared / static_cast<double>(association.size()));
    const double allowed = options.absoluteTolerance +
                           options.relativeTolerance * result.referenceScale;
    result.passed = finite(result.maxAbsoluteError) && result.maxAbsoluteError <= allowed;
    std::ostringstream explanation;
    explanation << (result.passed ? "PASS" : "FAIL")
                << " matched=" << result.matchedNodes
                << " max_abs_error=" << result.maxAbsoluteError
                << " rms_error=" << result.rmsError
                << " reference_scale=" << result.referenceScale
                << " limit=" << allowed
                << " gauge_offset=" << result.appliedGaugeOffset
                << " worst_reference_node=" << result.worstReferenceNode;
    result.message = explanation.str();
    return result;
}

EquivalenceResult compareSessionSolutions(femm::AnalysisSession &materialized,
                                           femm::AnalysisSession &native,
                                           const EquivalenceOptions &options)
{
    if (materialized.nativeInstanced() || !native.nativeInstanced())
        return failure("expected materialized and native session modes respectively");
    if (!materialized.instancedMesh() || !native.instancedMesh())
        return failure("both sessions need an instanced model");

    try {
        const femm::TrialSolution reference = materialized.solve();
        const femm::TrialSolution trial = native.solve();
        const auto mesh = materialized.mesh();
        const auto view = native.logicalView();
        if (!mesh || !view || !view->valid())
            return failure("missing materialized mesh or valid native logical view");
        if (mesh->nodes.size() != view->nodeCount() ||
            mesh->elements.size() != view->elementCount() ||
            mesh->periodicConstraints.size() != view->periodicConstraints().size() ||
            mesh->airGaps.size() != view->airGaps().size())
            return failure("logical and materialized topology counts differ");
        if (native.materializationCount() != 0)
            return failure("native path unexpectedly materialized a mesh");

        for (std::size_t i = 0; i < mesh->nodes.size(); ++i) {
            double x = 0, y = 0;
            view->nodeCoordinates(i, x, y);
            if (std::hypot(x - mesh->nodes[i].x, y - mesh->nodes[i].y) >
                options.coordinateTolerance)
                return failure("canonical logical/materialized node coordinates differ at " +
                               std::to_string(i));
        }
        for (std::size_t i = 0; i < mesh->elements.size(); ++i) {
            if (view->elementNodes(i) != mesh->elements[i].nodes)
                return failure("canonical logical/materialized connectivity differs at element " +
                               std::to_string(i));
        }
        return compareNodalSolutions(reference, trial, options);
    } catch (const std::exception &e) {
        return failure(std::string("equivalence session failed: ") + e.what());
    }
}

} // namespace xfemm_test
