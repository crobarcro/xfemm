#ifndef XFEMM_TEST_NATIVE_CIRCUIT_ENERGY_EQUIVALENCE_H
#define XFEMM_TEST_NATIVE_CIRCUIT_ENERGY_EQUIVALENCE_H

#include "NativeMatrixFieldEquivalence.h"

#include <cstddef>
#include <string>
#include <vector>

class FSolver;
namespace femm { class FemmProblem; }

namespace xfemm_test {

struct CircuitEnergyOptions {
    MatrixFieldOptions field;
    double fluxAbsoluteTolerance = 1e-9; // weber-turns
    double fluxRelativeTolerance = 1e-5;
    double energyAbsoluteTolerance = 1e-9; // joules
    double energyRelativeTolerance = 1e-5;
    bool requireNonzeroEnergy = true;
    bool requireNonzeroLinkage = true;
    double minimumSignal = 1e-12;
};

struct CircuitMeasurement {
    double current = 0;
    double fluxLinkage = 0;
};

struct CircuitEnergyMeasurements {
    std::vector<CircuitMeasurement> circuits;
    double energy = 0;
    double coenergy = 0;
};

struct CircuitEnergyResult {
    bool passed = false;
    std::string message;
    MatrixFieldResult baseline;
    CircuitEnergyMeasurements reference;
    CircuitEnergyMeasurements native;
    double maximumFluxDifference = 0;
    double energyDifference = 0;
    double coenergyDifference = 0;
};

/** Test-only adapter that runs the production FPProc integrals on the given
 * solved field, with a common expanded geometry and the *specified* solver's
 * own circuit/label state. This accepts solver-order potentials, not logical
 * node indices, and intentionally allocates no mesh in the native backend.
 *
 * In particular, the solver's expanded series sub-circuits are re-associated
 * with the original user circuit before GetFluxLinkage: otherwise instanced
 * series coils can incorrectly report zero flux linkage.
 */
CircuitEnergyMeasurements measureCircuitEnergy(
    const femm::FemmProblem &problem,
    const FSolver &expandedGeometry,
    const FSolver &driveState,
    const std::vector<double> &potentials);

/** G5.3: full independent native/materialized G5.2 comparison followed by
 * circuit linkage, magnetic energy and coenergy checks. Solves only once per
 * path. Throws no exceptions; returns actionable diagnostics.
 */
CircuitEnergyResult compareCircuitEnergySessions(
    femm::AnalysisSession &referenceSession,
    femm::FSolverAnalysisBackend &referenceSolver,
    femm::AnalysisSession &nativeSession,
    femm::FSolverAnalysisBackend &nativeSolver,
    const CircuitEnergyOptions &options = {});

} // namespace xfemm_test

#endif
