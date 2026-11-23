#include "NativeForceTorqueEquivalence.h"

#include "fpproc.h"
#include "fsolver.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace xfemm_test {
namespace {

ForceTorqueResult failure(const std::string &message)
{
    ForceTorqueResult result;
    result.message = message;
    return result;
}

void validate(const ForceTorqueOptions &options)
{
    if (!std::isfinite(options.absoluteForceTolerance) ||
        !std::isfinite(options.relativeForceTolerance) ||
        !std::isfinite(options.absoluteTorqueTolerance) ||
        !std::isfinite(options.relativeTorqueTolerance) ||
        !std::isfinite(options.minimumForceSignal) ||
        !std::isfinite(options.minimumTorqueSignal) ||
        options.absoluteForceTolerance < 0 || options.relativeForceTolerance < 0 ||
        options.absoluteTorqueTolerance < 0 || options.relativeTorqueTolerance < 0 ||
        options.minimumForceSignal < 0 || options.minimumTorqueSignal < 0)
        throw std::invalid_argument("invalid G5.4 force/torque tolerances");
    if (!options.compareLorentz && !options.compareWeightedStress &&
        options.airGapBoundary.empty())
        throw std::invalid_argument("G5.4 requires at least one physical force/torque observable");
    if ((options.compareLorentz || options.compareWeightedStress) &&
        options.selectedLabels.empty())
        throw std::invalid_argument("G5.4 selected force/torque blocks cannot be empty");
    if (options.requireNonzeroLorentz && !options.compareLorentz)
        throw std::invalid_argument("G5.4 nontrivial Lorentz check not enabled");
    if (options.requireNonzeroWeightedStress && !options.compareWeightedStress)
        throw std::invalid_argument("G5.4 nontrivial weighted stress check not enabled");
    if (options.requireNonzeroAirGapTorque && options.airGapBoundary.empty())
        throw std::invalid_argument("G5.4 nontrivial AGE torque check not enabled");
}

double finiteReal(CComplex value, const char *name)
{
    if (!std::isfinite(value.re) || !std::isfinite(value.im) ||
        std::abs(value.im) > 1e-9 * (1 + std::abs(value.re)))
        throw std::runtime_error(std::string("invalid real static ") + name);
    return value.re;
}

std::vector<double> reorderedNativePotential(
    const femm::TrialSolution &reference, const femm::TrialSolution &native,
    double coordinateTolerance)
{
    if (!reference.real || !native.real)
        throw std::runtime_error("missing real G5.4 nodal field");
    std::vector<std::size_t> mapping;
    std::string explanation;
    if (!associateNodalCoordinates(reference.real->nodal, native.real->nodal,
                                   coordinateTolerance, mapping, explanation))
        throw std::runtime_error("G5.4 cannot map native DOFs: " + explanation);
    std::vector<double> potential(mapping.size(), 0.);
    for (std::size_t i = 0; i < mapping.size(); ++i)
        potential[mapping[i]] = native.real->nodal.magneticVectorPotential[i];
    return potential;
}

void checkDifference(const char *name, double reference, double native,
                     double absolute, double relative,
                     double &maximumDifference)
{
    if (!std::isfinite(reference) || !std::isfinite(native))
        throw std::runtime_error(std::string("nonfinite G5.4 ") + name);
    const double difference = std::abs(reference - native);
    maximumDifference = std::max(maximumDifference, difference);
    const double allowed = absolute +
        relative * std::max(std::abs(reference), std::abs(native));
    if (difference > allowed) {
        std::ostringstream diagnostic;
        diagnostic.precision(17);
        diagnostic << "G5.4 " << name << " differs: materialized=" << reference
                   << " native=" << native << " delta=" << difference
                   << " limit=" << allowed;
        throw std::runtime_error(diagnostic.str());
    }
}

} // namespace

ForceTorqueObservations measureForceTorque(
    const femm::FemmProblem &problem, const FSolver &expandedGeometry,
    const FSolver &driveState, const std::vector<double> &potentials,
    const ForceTorqueOptions &options)
{
    validate(options);
    FPProc postprocessor;
    prepareCircuitPostprocessor(postprocessor, problem, expandedGeometry,
                                driveState, potentials);
    // FPProc::OpenDocument may mark some labels as selected on import.
    // Always select exactly the requested physical regions for both paths.
    for (auto &label : postprocessor.blocklist)
        label.IsSelected = false;
    std::vector<bool> seen(postprocessor.blocklist.size(), false);
    for (std::size_t index : options.selectedLabels) {
        if (index >= seen.size() || seen[index])
            throw std::invalid_argument("invalid or duplicate G5.4 selected label " +
                                        std::to_string(index));
        seen[index] = true;
        postprocessor.blocklist[index].IsSelected = true;
    }

    ForceTorqueObservations result;
    if (options.compareLorentz) {
        result.lorentzFx = finiteReal(postprocessor.BlockIntegral(11), "Lorentz Fx");
        result.lorentzFy = finiteReal(postprocessor.BlockIntegral(12), "Lorentz Fy");
        result.lorentzTorque = finiteReal(postprocessor.BlockIntegral(15), "Lorentz torque");
    }
    if (options.compareWeightedStress) {
        // FEMM weighted-stress tensor integrals require a separately solved
        // scalar mask; without MakeMask a test can silently report zero.
        if (!postprocessor.MakeMask() || !postprocessor.bHasMask)
            throw std::runtime_error("FEMM weighted-stress selection/mask is invalid");
        result.stressFx = finiteReal(postprocessor.BlockIntegral(18), "stress Fx");
        result.stressFy = finiteReal(postprocessor.BlockIntegral(19), "stress Fy");
        result.stressTorque = finiteReal(postprocessor.BlockIntegral(22), "stress torque");
    }
    if (!options.airGapBoundary.empty()) {
        if (postprocessor.gapDCTorqueIntegral(options.airGapBoundary,
                                              result.airGapTorque) != FPProcError::NoError)
            throw std::runtime_error("FEMM AGE torque integral failed for boundary '" +
                                     options.airGapBoundary + "'");
        if (!std::isfinite(result.airGapTorque))
            throw std::runtime_error("non-finite AGE torque");
    }
    // All FEMM torques are about the shared global (0,0) origin.
    return result;
}

ForceTorqueResult compareForceTorqueSessions(
    femm::AnalysisSession &referenceSession,
    femm::FSolverAnalysisBackend &referenceSolver,
    femm::AnalysisSession &nativeSession,
    femm::FSolverAnalysisBackend &nativeSolver,
    const ForceTorqueOptions &options,
    femm::TrialSolution *referenceOutput,
    femm::TrialSolution *nativeOutput)
{
    try {
        validate(options);
        ForceTorqueResult result;
        femm::TrialSolution referenceTrial, nativeTrial;
        result.baseline = compareCircuitEnergySessions(
            referenceSession, referenceSolver, nativeSession, nativeSolver,
            options.baseline, &referenceTrial, &nativeTrial);
        if (!result.baseline.passed)
            return failure("G5.3 baseline failed: " + result.baseline.message);

        const auto &problem = referenceSession.model().problem();
        const auto &geometry = referenceSolver.solvedSolver();
        result.reference = measureForceTorque(
            problem, geometry, geometry,
            referenceTrial.real->nodal.magneticVectorPotential, options);
        result.native = measureForceTorque(
            problem, geometry, nativeSolver.solvedSolver(),
            reorderedNativePotential(referenceTrial, nativeTrial,
                                     options.baseline.field.nodal.coordinateTolerance),
            options);

        const auto compareForce = [&](const char *name, double r, double n) {
            checkDifference(name, r, n, options.absoluteForceTolerance,
                            options.relativeForceTolerance, result.maximumForceDifference);
        };
        const auto compareTorque = [&](const char *name, double r, double n) {
            checkDifference(name, r, n, options.absoluteTorqueTolerance,
                            options.relativeTorqueTolerance, result.maximumTorqueDifference);
        };
        if (options.compareLorentz) {
            compareForce("Lorentz Fx", result.reference.lorentzFx, result.native.lorentzFx);
            compareForce("Lorentz Fy", result.reference.lorentzFy, result.native.lorentzFy);
            compareTorque("Lorentz torque", result.reference.lorentzTorque,
                          result.native.lorentzTorque);
            if (options.requireNonzeroLorentz &&
                std::hypot(result.reference.lorentzFx, result.reference.lorentzFy) <=
                options.minimumForceSignal &&
                std::abs(result.reference.lorentzTorque) <= options.minimumTorqueSignal)
                throw std::runtime_error("Lorentz force and torque are both trivial");
        }
        if (options.compareWeightedStress) {
            compareForce("weighted stress Fx", result.reference.stressFx, result.native.stressFx);
            compareForce("weighted stress Fy", result.reference.stressFy, result.native.stressFy);
            compareTorque("weighted stress torque", result.reference.stressTorque,
                          result.native.stressTorque);
            if (options.requireNonzeroWeightedStress &&
                std::hypot(result.reference.stressFx, result.reference.stressFy) <=
                options.minimumForceSignal &&
                std::abs(result.reference.stressTorque) <= options.minimumTorqueSignal)
                throw std::runtime_error("weighted stress force and torque are both trivial");
        }
        if (!options.airGapBoundary.empty()) {
            compareTorque("air-gap torque", result.reference.airGapTorque,
                          result.native.airGapTorque);
            if (options.requireNonzeroAirGapTorque &&
                std::abs(result.reference.airGapTorque) <= options.minimumTorqueSignal)
                throw std::runtime_error("air-gap torque is trivial");
        }

        // G5.5 can verify constitutive H, including Lua magnetisation,
        // without invoking another solve or assuming native DOF ordering.
        if (referenceOutput) *referenceOutput = referenceTrial;
        if (nativeOutput) *nativeOutput = nativeTrial;
        result.passed = true;
        std::ostringstream summary;
        summary.precision(12);
        summary << "PASS max_force_delta_N=" << result.maximumForceDifference
                << " max_torque_delta_Nm=" << result.maximumTorqueDifference
                << " lorentz_F_ref=(" << result.reference.lorentzFx
                << "," << result.reference.lorentzFy
                << ") stress_F_ref=(" << result.reference.stressFx
                << "," << result.reference.stressFy
                << ") AGE_torque_ref=" << result.reference.airGapTorque;
        result.message = summary.str();
        return result;
    } catch (const std::exception &e) {
        return failure(std::string("force/torque equivalence failed: ") + e.what());
    }
}

} // namespace xfemm_test
