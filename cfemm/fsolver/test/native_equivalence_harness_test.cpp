#include "NativeEquivalence.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace {

femm::TrialSolution makeSolution(std::vector<double> x, std::vector<double> y,
                                  std::vector<double> potential)
{
    femm::TrialSolution solution;
    solution.real.emplace();
    solution.real->nodal.x = std::move(x);
    solution.real->nodal.y = std::move(y);
    solution.real->nodal.magneticVectorPotential = std::move(potential);
    return solution;
}

bool check(bool condition, const std::string &name)
{
    if (!condition)
        std::cerr << "FAIL: " << name << '\n';
    return condition;
}

} // namespace

int main()
{
    using xfemm_test::compareNodalSolutions;
    using xfemm_test::EquivalenceOptions;

    // The two solvers may return completely different node orderings.
    const auto reference = makeSolution({0., 1., 0.}, {0., 0., 1.}, {2., 3., 4.});
    const auto permuted = makeSolution({0., 0., 1.}, {1., 0., 0.}, {4., 2., 3.});
    EquivalenceOptions exact;
    exact.absoluteTolerance = 1e-12;
    exact.relativeTolerance = 0;
    if (!check(compareNodalSolutions(reference, permuted, exact).passed,
               "permuted fields agree")) return 1;

    auto perturbed = permuted;
    perturbed.real->nodal.magneticVectorPotential[1] += 1e-3;
    auto bad = compareNodalSolutions(reference, perturbed, exact);
    if (!check(!bad.passed && bad.matchedNodes == 3 &&
               bad.maxAbsoluteError > 9e-4 && bad.worstReferenceNode == 0,
               "field mismatch is detected and attributed")) return 1;

    auto shifted = permuted;
    for (double &a : shifted.real->nodal.magneticVectorPotential)
        a += 7.;
    if (!check(!compareNodalSolutions(reference, shifted, exact).passed,
               "strict potential comparison rejects a gauge shift")) return 1;
    exact.removeGaugeOffset = true;
    auto gauge = compareNodalSolutions(reference, shifted, exact);
    if (!check(gauge.passed && std::abs(gauge.appliedGaugeOffset + 7.) < 1e-12,
               "optional gauge alignment preserves equivalence")) return 1;

    // No coordinate map may silently overwrite coincident independent DOFs.
    const auto coincident = makeSolution({0., 0.}, {0., 0.}, {1., 2.});
    if (!check(!compareNodalSolutions(coincident, coincident, exact).passed,
               "ambiguous coincident nodes fail closed")) return 1;

    auto nan = permuted;
    nan.real->nodal.magneticVectorPotential[0] =
        std::numeric_limits<double>::quiet_NaN();
    if (!check(!compareNodalSolutions(reference, nan, exact).passed,
               "non-finite solution values are rejected")) return 1;

    auto missing = permuted;
    missing.real->nodal.x[0] = 0.25;
    if (!check(!compareNodalSolutions(reference, missing, exact).passed,
               "unmatched coordinates are rejected")) return 1;

    const auto zero = makeSolution({0., 1.}, {0., 0.}, {0., 0.});
    if (!check(!compareNodalSolutions(zero, zero, exact).passed,
               "zero-excitation reference is not an equivalence proof")) return 1;
    exact.requireNontrivialField = false;
    if (!check(compareNodalSolutions(zero, zero, exact).passed,
               "zero field can be compared when explicitly requested")) return 1;

    std::cout << "native equivalence harness self-tests passed\n";
    return 0;
}
