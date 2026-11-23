#include "NativeCircuitEnergyEquivalence.h"

#include "CBlockLabel.h"
#include "CCircuit.h"
#include "fpproc.h"
#include "fsolver.h"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace xfemm_test {
namespace {

CircuitEnergyResult failure(const std::string &message)
{
    CircuitEnergyResult result;
    result.message = message;
    return result;
}

bool acceptable(double reference, double candidate, double absTol, double relTol)
{
    return std::isfinite(reference) && std::isfinite(candidate) &&
           std::abs(reference - candidate) <=
               absTol + relTol * std::max(std::abs(reference), std::abs(candidate));
}

void verifyDriveConsistency(const FSolver &reference, const FSolver &native,
                            std::size_t originalCircuitCount)
{
    if (reference.NumCircPropsOrig != static_cast<int>(originalCircuitCount) ||
        native.NumCircPropsOrig != static_cast<int>(originalCircuitCount) ||
        reference.labellist.size() != native.labellist.size() ||
        reference.blockproplist.size() != native.blockproplist.size())
        throw std::runtime_error("native/materialized circuit or region metadata counts differ");

    for (std::size_t i = 0; i < originalCircuitCount; ++i) {
        const auto &a = reference.circproplist[i];
        const auto &b = native.circproplist[i];
        if (a.CircName != b.CircName || a.CircType != b.CircType ||
            a.Amps.re != b.Amps.re || a.Amps.im != b.Amps.im)
            throw std::runtime_error("native/materialized drive state differs for circuit " +
                                     std::to_string(i));
    }

    for (std::size_t i = 0; i < reference.labellist.size(); ++i) {
        const auto &a = reference.labellist[i];
        const auto &b = native.labellist[i];
        if (a.BlockType != b.BlockType || a.Turns != b.Turns ||
            a.MagDir != b.MagDir || a.MagDirFctn != b.MagDirFctn)
            throw std::runtime_error("native/materialized region physics differs at label " +
                                     std::to_string(i));
        const auto originalCircuit = [](const FSolver &source, int driven) {
            if (driven < 0) return -1;
            if (static_cast<std::size_t>(driven) >= source.circproplist.size())
                throw std::runtime_error("invalid expanded circuit index");
            const int orig = source.circproplist[driven].OrigCirc;
            return orig >= 0 ? orig : driven;
        };
        if (originalCircuit(reference, a.InCircuit) !=
            originalCircuit(native, b.InCircuit))
            throw std::runtime_error("native/materialized circuit membership differs at label " +
                                     std::to_string(i));
    }
    for (std::size_t i = 0; i < reference.blockproplist.size(); ++i) {
        const auto &a = reference.blockproplist[i];
        const auto &b = native.blockproplist[i];
        if (a.BlockName != b.BlockName || a.mu_x != b.mu_x || a.mu_y != b.mu_y ||
            a.H_c != b.H_c || a.BHpoints != b.BHpoints ||
            a.Bdata != b.Bdata || a.Hdata != b.Hdata)
            throw std::runtime_error("native/materialized material differs at index " +
                                     std::to_string(i));
    }
}

void configureUserCircuits(FPProc &post, const femm::FemmProblem &problem,
                           const FSolver &source)
{
    if (post.blocklist.size() != source.labellist.size() ||
        source.NumCircPropsOrig != static_cast<int>(problem.circproplist.size()))
        throw std::runtime_error("post-processor circuit/label cardinality differs");

    // OpenDocument may use solver-expanded circuits for an instanced model.
    // Restore original circuit types and user-visible IDs without discarding
    // solver-prepared per-label current-density and voltage-gradient state.
    post.circproplist.clear();
    for (std::size_t i = 0; i < problem.circproplist.size(); ++i) {
        const auto *original = dynamic_cast<const femm::CMCircuit *>(
            problem.circproplist[i].get());
        if (!original)
            throw std::runtime_error("post-processor requires magnetic circuits");
        auto current = *original;
        current.Amps = source.circproplist[i].Amps;
        post.circproplist.push_back(current);
    }

    for (std::size_t i = 0; i < post.blocklist.size(); ++i) {
        const auto &sourceLabel = source.labellist[i];
        auto &label = post.blocklist[i];
        label.Turns = sourceLabel.Turns;
        label.BlockType = sourceLabel.BlockType;
        const int driven = sourceLabel.InCircuit;
        label.InCircuit = -1;
        label.Case = 1;
        label.J = 0.;
        if (driven >= 0) {
            if (static_cast<std::size_t>(driven) >= source.circproplist.size())
                throw std::runtime_error("invalid expanded circuit ID on label");
            const auto &expanded = source.circproplist[driven];
            label.InCircuit = expanded.OrigCirc >= 0 ? expanded.OrigCirc : driven;
            if (label.InCircuit < 0 ||
                static_cast<std::size_t>(label.InCircuit) >= post.circproplist.size())
                throw std::runtime_error("expanded circuit cannot map to user circuit");
            if (expanded.Case == 0) {
                label.Case = 0;
                label.dVolts = expanded.dVolts;
            } else {
                label.Case = 1;
                label.J = expanded.J;
            }
        }
        label.IsSelected = true;
    }
}

std::vector<double> nativeInReferenceOrder(const femm::TrialSolution &reference,
                                            const femm::TrialSolution &native,
                                            double tolerance)
{
    if (!reference.real || !native.real)
        throw std::runtime_error("missing real trial solution for G5.3");
    std::vector<std::size_t> nativeToReference;
    std::string diagnostic;
    if (!associateNodalCoordinates(reference.real->nodal, native.real->nodal,
                                   tolerance, nativeToReference, diagnostic))
        throw std::runtime_error("G5.3 node mapping failed: " + diagnostic);
    std::vector<double> values(reference.real->nodal.magneticVectorPotential.size());
    for (std::size_t i = 0; i < nativeToReference.size(); ++i)
        values[nativeToReference[i]] = native.real->nodal.magneticVectorPotential[i];
    return values;
}

} // namespace

CircuitEnergyMeasurements measureCircuitEnergy(
    const femm::FemmProblem &problem,
    const FSolver &expandedGeometry,
    const FSolver &driveState,
    const std::vector<double> &potentials)
{
    if (problem.problemType != femm::ProblemType::PLANAR || problem.Frequency != 0)
        throw std::invalid_argument("G5.3 only supports planar static problems");
    if (potentials.size() != expandedGeometry.meshnode.size() ||
        potentials.size() != static_cast<std::size_t>(expandedGeometry.NumNodes))
        throw std::invalid_argument("post-processor potential and mesh lengths differ");

    FPProc postprocessor;
    if (!postprocessor.OpenDocument(problem, expandedGeometry, potentials))
        throw std::runtime_error("FPProc could not load the solved field");
    configureUserCircuits(postprocessor, problem, driveState);

    CircuitEnergyMeasurements result;
    result.circuits.reserve(problem.circproplist.size());
    for (std::size_t i = 0; i < problem.circproplist.size(); ++i) {
        const CComplex linkage = postprocessor.GetFluxLinkage(static_cast<int>(i));
        const CComplex current = postprocessor.circproplist[i].Amps;
        if (!std::isfinite(linkage.re) || !std::isfinite(linkage.im) ||
            !std::isfinite(current.re) || !std::isfinite(current.im) ||
            std::abs(linkage.im) > 1e-9 * (1 + std::abs(linkage.re)) ||
            current.im != 0)
            throw std::runtime_error("non-finite or complex static circuit linkage/current " +
                                     std::to_string(i));
        result.circuits.push_back({current.re, linkage.re});
    }

    const CComplex energy = postprocessor.BlockIntegral(2);
    const CComplex coenergy = postprocessor.BlockIntegral(17);
    if (!std::isfinite(energy.re) || !std::isfinite(coenergy.re) ||
        !std::isfinite(energy.im) || !std::isfinite(coenergy.im) ||
        std::abs(energy.im) > 1e-9 * (1 + std::abs(energy.re)) ||
        std::abs(coenergy.im) > 1e-9 * (1 + std::abs(coenergy.re)))
        throw std::runtime_error("non-finite or complex static energy integral");
    result.energy = energy.re;
    result.coenergy = coenergy.re;
    return result;
}

CircuitEnergyResult compareCircuitEnergySessions(
    femm::AnalysisSession &referenceSession,
    femm::FSolverAnalysisBackend &referenceSolver,
    femm::AnalysisSession &nativeSession,
    femm::FSolverAnalysisBackend &nativeSolver,
    const CircuitEnergyOptions &options)
{
    if (!std::isfinite(options.fluxAbsoluteTolerance) ||
        !std::isfinite(options.fluxRelativeTolerance) ||
        !std::isfinite(options.energyAbsoluteTolerance) ||
        !std::isfinite(options.energyRelativeTolerance) ||
        !std::isfinite(options.minimumSignal) ||
        options.fluxAbsoluteTolerance < 0 || options.fluxRelativeTolerance < 0 ||
        options.energyAbsoluteTolerance < 0 || options.energyRelativeTolerance < 0 ||
        options.minimumSignal < 0)
        return failure("invalid circuit/energy equivalence options");

    try {
        CircuitEnergyResult result;
        femm::TrialSolution referenceTrial, nativeTrial;
        result.baseline = compareMatrixFieldSessions(
            referenceSession, referenceSolver, nativeSession, nativeSolver,
            options.field, &referenceTrial, &nativeTrial);
        if (!result.baseline.passed)
            return failure("G5.2 baseline failed: " + result.baseline.message);

        const auto &ref = referenceSolver.solvedSolver();
        const auto &native = nativeSolver.solvedSolver();
        const auto &problem = referenceSession.model().problem();
        verifyDriveConsistency(ref, native, problem.circproplist.size());

        const auto mapped = nativeInReferenceOrder(
            referenceTrial, nativeTrial, options.field.nodal.coordinateTolerance);
        result.reference = measureCircuitEnergy(
            problem, ref, ref, referenceTrial.real->nodal.magneticVectorPotential);
        result.native = measureCircuitEnergy(problem, ref, native, mapped);
        if (result.reference.circuits.size() != result.native.circuits.size())
            return failure("G5.3 circuit count differs");

        double maximumLinkage = 0;
        for (std::size_t i = 0; i < result.reference.circuits.size(); ++i) {
            const auto &a = result.reference.circuits[i];
            const auto &b = result.native.circuits[i];
            const double difference = std::abs(a.fluxLinkage - b.fluxLinkage);
            result.maximumFluxDifference = std::max(result.maximumFluxDifference, difference);
            maximumLinkage = std::max(maximumLinkage, std::abs(a.fluxLinkage));
            if (!acceptable(a.current, b.current, 0, 0))
                return failure("different circuit currents at index " + std::to_string(i));
            if (!acceptable(a.fluxLinkage, b.fluxLinkage,
                            options.fluxAbsoluteTolerance, options.fluxRelativeTolerance)) {
                std::ostringstream diagnostic;
                diagnostic << "flux linkage differs for circuit " << i
                           << ": materialized=" << a.fluxLinkage
                           << " native=" << b.fluxLinkage;
                return failure(diagnostic.str());
            }
        }

        result.energyDifference = std::abs(result.reference.energy - result.native.energy);
        result.coenergyDifference = std::abs(result.reference.coenergy - result.native.coenergy);
        if (!acceptable(result.reference.energy, result.native.energy,
                        options.energyAbsoluteTolerance, options.energyRelativeTolerance))
            return failure("magnetic energy differs: materialized=" +
                           std::to_string(result.reference.energy) +
                           " native=" + std::to_string(result.native.energy));
        if (!acceptable(result.reference.coenergy, result.native.coenergy,
                        options.energyAbsoluteTolerance, options.energyRelativeTolerance))
            return failure("magnetic coenergy differs: materialized=" +
                           std::to_string(result.reference.coenergy) +
                           " native=" + std::to_string(result.native.coenergy));
        if (options.requireNonzeroLinkage && maximumLinkage <= options.minimumSignal)
            return failure("all circuit flux linkages are trivial");
        if (options.requireNonzeroEnergy &&
            std::max(std::abs(result.reference.energy),
                     std::abs(result.reference.coenergy)) <= options.minimumSignal)
            return failure("energy and coenergy are both trivial");

        result.passed = true;
        std::ostringstream summary;
        summary << "PASS circuits=" << result.reference.circuits.size()
                << " max_flux_delta_Wb_turn=" << result.maximumFluxDifference
                << " energy_ref_J=" << result.reference.energy
                << " energy_native_J=" << result.native.energy
                << " coenergy_ref_J=" << result.reference.coenergy
                << " coenergy_native_J=" << result.native.coenergy;
        result.message = summary.str();
        return result;
    } catch (const std::exception &e) {
        return failure(std::string("G5.3 circuit/energy equivalence failed: ") + e.what());
    }
}

} // namespace xfemm_test
