#include "NativeMatrixFieldEquivalence.h"

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace {
using Snapshot = femm::FSolverAnalysisBackend::StaticAssemblySnapshot;

bool require(bool condition, const std::string &name, const std::string &detail)
{
    if (!condition)
        std::cerr << "FAIL: " << name << ": " << detail << '\n';
    return condition;
}

Snapshot makeReference()
{
    Snapshot s;
    s.dimension = 3;
    s.capturedIterations = 1;
    s.upperTriangle = {{0, 0, 4}, {0, 1, -1}, {1, 1, 3},
                       {1, 2, -1}, {2, 2, 2}};
    s.rhs = {2, 2, 4};
    return s;
}

Snapshot makePermuted()
{
    // Native order [reference 2, reference 0, reference 1].
    Snapshot s;
    s.dimension = 3;
    s.capturedIterations = 1;
    s.upperTriangle = {{0, 0, 2}, {0, 2, -1}, {1, 1, 4},
                       {1, 2, -1}, {2, 2, 3}};
    s.rhs = {4, 2, 2};
    return s;
}
}

int main()
{
    const Snapshot a = makeReference();
    const Snapshot b = makePermuted();
    const std::vector<std::size_t> mapping{2, 0, 1};
    const std::vector<double> x{1, 2, 3}, y{3, 1, 2};

    xfemm_test::MatrixFieldOptions options;
    options.matrixAbsoluteTolerance = options.rhsAbsoluteTolerance = 1e-12;
    options.matrixRelativeTolerance = options.rhsRelativeTolerance = 0;
    options.maximumRelativeResidual = 1e-12;

    auto result = xfemm_test::compareAssemblySnapshots(a, b, mapping, x, y, options);
    if (!require(result.passed && result.numericNonzeros == 5 &&
                 result.maximumMatrixDifference == 0 &&
                 result.referenceRelativeResidual < 1e-12 &&
                 result.nativeRelativeResidual < 1e-12,
                 "permuted assembly", result.message)) return 1;

    auto wrongRhs = b;
    wrongRhs.rhs[0] += 0.1;
    result = xfemm_test::compareAssemblySnapshots(a, wrongRhs, mapping, x, y, options);
    if (!require(!result.passed && result.message.find("residual") != std::string::npos,
                 "residual detects perturbed RHS", result.message)) return 1;

    auto wrongCoefficient = b;
    wrongCoefficient.upperTriangle[0].value += 0.25;
    options.maximumRelativeResidual = 1;
    result = xfemm_test::compareAssemblySnapshots(a, wrongCoefficient, mapping, x, y, options);
    if (!require(!result.passed && result.message.find("matrix coefficient") != std::string::npos,
                 "matrix detects changed coefficient", result.message)) return 1;

    auto wrongSparsity = b;
    wrongSparsity.upperTriangle.push_back({0, 1, 0.25});
    result = xfemm_test::compareAssemblySnapshots(a, wrongSparsity, mapping, x, y, options);
    if (!require(!result.passed && result.sparsityMismatches == 1,
                 "matrix detects unexpected sparse edge", result.message)) return 1;

    const std::vector<std::size_t> duplicate{2, 2, 1};
    result = xfemm_test::compareAssemblySnapshots(a, b, duplicate, x, y, options);
    if (!require(!result.passed && result.message.find("mapping") != std::string::npos,
                 "matrix rejects ambiguous DOFs", result.message)) return 1;

    options.compareMatricesAndRhs = false;
    options.maximumRelativeResidual = 1e-12;
    result = xfemm_test::compareAssemblySnapshots(a, b, mapping, x, y, options);
    if (!require(result.passed, "residual-only comparison", result.message)) return 1;

    std::cout << "G5.2 matrix/residual self-tests passed\n";
    return 0;
}
