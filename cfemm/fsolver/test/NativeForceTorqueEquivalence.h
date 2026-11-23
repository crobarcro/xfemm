#ifndef XFEMM_TEST_NATIVE_FORCE_TORQUE_EQUIVALENCE_H
#define XFEMM_TEST_NATIVE_FORCE_TORQUE_EQUIVALENCE_H

#include "NativeCircuitEnergyEquivalence.h"

#include <cstddef>
#include <string>
#include <vector>

namespace xfemm_test {

/** Independent force/torque comparisons in SI units. Selection identifies
 * prepared per-instance block-label indices, not material indices.
 */
struct ForceTorqueOptions {
    CircuitEnergyOptions baseline;
    std::vector<std::size_t> selectedLabels;
    bool compareLorentz = false;         // FPProc BlockIntegral(11, 12, 15)
    bool compareWeightedStress = false;  // FPProc BlockIntegral(18, 19, 22)
    std::string airGapBoundary;         // FPProc gapDCTorqueIntegral
    double absoluteForceTolerance = 1e-6; // newtons
    double relativeForceTolerance = 1e-3;
    double absoluteTorqueTolerance = 1e-7; // newton metres
    double relativeTorqueTolerance = 1e-3;
    bool requireNonzeroLorentz = false;
    bool requireNonzeroWeightedStress = false;
    bool requireNonzeroAirGapTorque = false;
    double minimumForceSignal = 1e-6;
    double minimumTorqueSignal = 1e-7;
};

struct ForceTorqueObservations {
    double lorentzFx = 0;
    double lorentzFy = 0;
    double lorentzTorque = 0;
    double stressFx = 0;
    double stressFy = 0;
    double stressTorque = 0;
    double airGapTorque = 0;
};

struct ForceTorqueResult {
    bool passed = false;
    std::string message;
    CircuitEnergyResult baseline;
    ForceTorqueObservations reference;
    ForceTorqueObservations native;
    double maximumForceDifference = 0;
    double maximumTorqueDifference = 0;
};

/** Read physical loads from an independently solved field, using the
 * materialized geometry only as a test-side FPProc oracle. Weighted-stress
 * masks are solved independently for reference/native. Air-gap quantities
 * use the same specified AGE and relative rotor angle for both solutions.
 */
ForceTorqueObservations measureForceTorque(
    const femm::FemmProblem &problem, const FSolver &expandedGeometry,
    const FSolver &driveState, const std::vector<double> &potentials,
    const ForceTorqueOptions &options);

/** Execute the G5.1-G5.3 baseline (one solve per path), then compare
 * Lorentz/weighted-stress force and torque and/or air-gap torque.
 */
ForceTorqueResult compareForceTorqueSessions(
    femm::AnalysisSession &referenceSession,
    femm::FSolverAnalysisBackend &referenceSolver,
    femm::AnalysisSession &nativeSession,
    femm::FSolverAnalysisBackend &nativeSolver,
    const ForceTorqueOptions &options = {});

} // namespace xfemm_test

#endif
