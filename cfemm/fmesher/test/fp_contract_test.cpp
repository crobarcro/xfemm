// Checks that floating-point contraction is off for the code that feeds
// Triangle: Triangle itself (C) and libfemm (C++, the geometry code that
// computes the coordinates handed to Triangle). arm64 Clang and aarch64
// GCC contract a*b+c into a fused multiply-add by default. With libfemm
// fused, Triangle's divide-and-conquer triangulation loops forever on the
// periodic and AGE models; with Triangle fused, its results drift
// (a3_2_mesher_differential fails). This test fails in milliseconds instead
// of waiting for that hang.
//
// a = 1 + 2^-30, b = 1 - 2^-30, c = -1:
//   separate rounding: a*b = 1 - 2^-60 rounds to 1.0, so a*b + c == 0
//   fused:             a*b + c == -2^-60 exactly
//
// A third probe compiled with -ffp-contract=fast is the control: it shows
// whether this compiler/CPU fuses at all. Where it does not (e.g. x86-64
// without -mfma) the check passes trivially and the test says so.

#include <cmath>
#include <cstdio>

extern "C" double xfemm_fp_probe_triangle(double, double, double);
extern "C" double xfemm_fp_probe_femm(double, double, double);
extern "C" double xfemm_fp_probe_control(double, double, double);

namespace {

volatile double va = 1.0 + std::ldexp(1.0, -30);
volatile double vb = 1.0 - std::ldexp(1.0, -30);
volatile double vc = -1.0;

bool fused(double (*probe)(double, double, double))
{
    return probe(va, vb, vc) != 0.0;
}

} // namespace

int main()
{
    const bool controlFuses = fused(xfemm_fp_probe_control);
    std::printf("control (-ffp-contract=fast) fuses: %s\n", controlFuses ? "yes" : "no");
    if (!controlFuses)
        std::printf("note: this compiler/CPU does not fuse multiply-adds here, "
                    "so the checks below cannot fail on this platform\n");

    int failures = 0;
    if (fused(xfemm_fp_probe_triangle)) {
        std::printf("FAIL: the triangle target contracts multiply-adds; "
                    "it needs -ffp-contract=off\n");
        ++failures;
    } else {
        std::printf("ok: triangle target keeps multiply and add separate\n");
    }
    if (fused(xfemm_fp_probe_femm)) {
        std::printf("FAIL: the femm (libfemm) target contracts multiply-adds; "
                    "it needs -ffp-contract=off\n");
        ++failures;
    } else {
        std::printf("ok: femm (libfemm) target keeps multiply and add separate\n");
    }
    return failures == 0 ? 0 : 1;
}
