#include "femmcomplex.h"
#include "cspars.h"

#include <array>
#include <cmath>
#include <complex>
#include <iostream>
#include <string>

namespace {

constexpr std::size_t Dimension = 3;
using Complex = std::complex<double>;
using Vector = std::array<Complex, Dimension>;
using Matrix = std::array<std::array<Complex, Dimension>, Dimension>;

Vector multiply(const Matrix &matrix, const Vector &vector)
{
    Vector result{};
    for (std::size_t row = 0; row < Dimension; ++row)
        for (std::size_t column = 0; column < Dimension; ++column)
            result[row] += matrix[row][column] * vector[column];
    return result;
}

Matrix conjugate(const Matrix &matrix)
{
    Matrix result{};
    for (std::size_t row = 0; row < Dimension; ++row)
        for (std::size_t column = 0; column < Dimension; ++column)
            result[row][column] = std::conj(matrix[row][column]);
    return result;
}

Vector conjugate(const Vector &vector)
{
    Vector result{};
    for (std::size_t i = 0; i < Dimension; ++i) result[i] = std::conj(vector[i]);
    return result;
}

Vector operator+(Vector left, const Vector &right)
{
    for (std::size_t i = 0; i < Dimension; ++i) left[i] += right[i];
    return left;
}

Vector operator-(Vector left, const Vector &right)
{
    for (std::size_t i = 0; i < Dimension; ++i) left[i] -= right[i];
    return left;
}

void storeUpperTriangle(CBigComplexLinProb &problem, const Matrix &matrix, int kind)
{
    for (std::size_t row = 0; row < Dimension; ++row)
        for (std::size_t column = row; column < Dimension; ++column)
            problem.Put(CComplex(matrix[row][column].real(), matrix[row][column].imag()),
                    static_cast<int>(row), static_cast<int>(column), kind);
}

bool check(const std::string &name, const Vector &expected, const CComplex *actual)
{
    constexpr double tolerance = 1.e-12;
    for (std::size_t i = 0; i < Dimension; ++i) {
        const Complex value(actual[i].re, actual[i].im);
        if (std::abs(value - expected[i]) > tolerance) {
            std::cerr << name << '[' << i << "]: expected " << expected[i]
                      << ", got " << value << '\n';
            return false;
        }
    }
    return true;
}

} // namespace

int main()
{
    const Matrix symmetric{{
        {{{2., 1.}, {1., -2.}, {-3., .5}}},
        {{{1., -2.}, {4., -1.}, {.25, 3.}}},
        {{{-3., .5}, {.25, 3.}, {-2., 2.}}}
    }};
    const Matrix hermitian{{
        {{{1.5, 0.}, {2., 3.}, {-1., .75}}},
        {{{2., -3.}, {-2.5, 0.}, {4., -1.}}},
        {{{-1., -.75}, {4., 1.}, {3.5, 0.}}}
    }};
    const Matrix antiHermitian{{
        {{{0., 2.}, {1., -4.}, {2.5, .5}}},
        {{{-1., -4.}, {0., -3.}, {-2., 1.5}}},
        {{{-2.5, .5}, {2., 1.5}, {0., .25}}}
    }};
    const Matrix auxiliarySymmetric{{
        {{{-.5, 1.}, {3., .25}, {1.25, -2.}}},
        {{{3., .25}, {2., .5}, {-1., -1.5}}},
        {{{1.25, -2.}, {-1., -1.5}, {4., -3.}}}
    }};
    const Vector input{{{1., -1.}, {-.5, 2.}, {3., .25}}};
    std::array<CComplex, Dimension> x;
    std::array<CComplex, Dimension> actual;
    for (std::size_t i = 0; i < Dimension; ++i)
        x[i] = CComplex(input[i].real(), input[i].imag());

    CBigComplexLinProb problem;
    problem.Create(static_cast<int>(Dimension), 0, static_cast<int>(Dimension));
    storeUpperTriangle(problem, symmetric, 0);
    storeUpperTriangle(problem, hermitian, 1);
    storeUpperTriangle(problem, auxiliarySymmetric, 2);
    storeUpperTriangle(problem, antiHermitian, 3);

    bool ok = true;
    const std::array<Matrix, 4> matrices{{symmetric, hermitian, auxiliarySymmetric,
                                         antiHermitian}};
    const std::array<std::string, 4> names{{"symmetric", "Hermitian",
                                            "auxiliary symmetric", "anti-Hermitian"}};
    for (int kind = 0; kind < 4; ++kind) {
        problem.MultA(x.data(), actual.data(), kind);
        ok &= check("MultA " + names[kind], multiply(matrices[kind], input), actual.data());
        problem.MultConjA(x.data(), actual.data(), kind);
        ok &= check("MultConjA " + names[kind],
                multiply(conjugate(matrices[kind]), input), actual.data());
    }

    problem.MultA(x.data(), actual.data(), -1);
    const Vector combinedMinusOne = multiply(symmetric, input)
            + multiply(hermitian, input)
            + conjugate(multiply(conjugate(auxiliarySymmetric), input))
            + multiply(antiHermitian, input);
    ok &= check("MultA k=-1", combinedMinusOne, actual.data());

    problem.MultA(x.data(), actual.data(), -2);
    const Vector auxiliaryProduct = multiply(auxiliarySymmetric, input);
    const Vector combinedMinusTwo = multiply(symmetric, input)
            + multiply(conjugate(hermitian), input)
            + auxiliaryProduct + conjugate(auxiliaryProduct)
            - multiply(conjugate(antiHermitian), input);
    ok &= check("MultA k=-2", combinedMinusTwo, actual.data());

    return ok ? 0 : 1;
}
