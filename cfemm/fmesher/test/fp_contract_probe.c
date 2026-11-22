/* One multiply-add, compiled with the options of a chosen target (see
 * CMakeLists.txt). If that target lets the compiler contract a*b+c into a
 * fused multiply-add, the result differs from the separately rounded one,
 * which fp_contract_test.cpp detects. Kept in its own translation unit so
 * the inputs cannot be constant-folded. The same file is compiled as C and
 * as C++ (fp_contract_probe.cpp includes it). */
#ifdef __cplusplus
extern "C"
#endif
double XFEMM_FP_PROBE_NAME(double a, double b, double c)
{
    return a * b + c;
}
