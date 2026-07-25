# Review of issue #18: sparse-solver performance

## Verdict

The performance claim is **plausible, but the proposed commit is not valid for
merge**. Replacing one heap allocation per matrix entry and linked-list traversal
with contiguous row storage can materially improve allocation cost and cache
locality in the matrix-vector product and SSOR preconditioner. Those operations
run on every PCG iteration, so a twofold solver-phase speedup on sufficiently
large matrices is credible. The issue does not, however, include a reproducible
benchmark, timings, mesh sizes, iteration counts, compiler settings, or output
comparisons. It therefore does not establish a measured twofold speedup for
xfemm as a whole.

Reviewed proposal:
[`Reoptimize-Systems/xfemm@f3dbaca`](https://github.com/Reoptimize-Systems/xfemm/commit/f3dbaca95c524c71742bac5597dea644d475edb8).

## Why the optimization direction is credible

The current representation allocates every non-diagonal `CEntry` separately and
links entries with `next`. `MultA` and both SSOR sweeps then pointer-chase across
each row. A row-wise contiguous representation removes most allocations, reduces
per-entry metadata, and improves spatial locality and hardware prefetching.
Since `PCGSolve` invokes both `MultA` and `MultPC` for every iteration, the
proposal targets the solver's hot path rather than incidental setup work.

The magnitude will vary. Small problems, low PCG iteration counts, assembly-heavy
or meshing-heavy workloads, and rows repeatedly inserted in the middle of a
`std::vector` will see less benefit. A contiguous format assembled once and then
frozen—preferably CSR, or row vectors with reserved capacity—is a sound direction.

## MATLAB, Octave, and MEX language-level compatibility

The language mode of a MEX source file and the MEX API/ABI version are separate
concerns. MATLAB's `-R2017b` and `-R2018a` switches select the complex-number API;
they do not select a C++ dialect. Likewise, `mex -setup C++` selects a compiler,
but the compiler's default dialect is not a promise that C++17 is enabled. A MEX
build that uses C++17 must pass the appropriate compiler option explicitly (for
example, `-std=c++17` with GCC/Clang or `/std:c++17` with MSVC) and must use a
compiler supported by that MATLAB release.

MathWorks' per-release
[supported-compiler tables](https://www.mathworks.com/support/requirements/supported-compilers.html)
are the authority for the accepted compiler versions on Windows, Linux, and
macOS. The compiler families accepted by MATLAB releases from 2023 onward have
C++17 implementations, so compiling a classic C/C++ Matrix API MEX file as
C++17 is technically reasonable on those releases. This does **not** mean that
MathWorks defines one universal "MEX C++ version," nor that `mex` enables C++17
by default. Users must check the table for their exact MATLAB release and
platform and configure it through
[`mex -setup C++`](https://www.mathworks.com/help/matlab/ref/mex.html).
Code using MATLAB's newer C++ Data API must additionally follow that API's
documented compiler and binary-compatibility requirements; xfemm's existing MEX
wrappers use the classic API, so that additional constraint does not currently
apply here.

GNU Octave has a more explicit progression in its build configuration:

- Octave 8 (2023) requires C++11 by default and optionally enables C++17 for
  `std::pmr::polymorphic_allocator`; its configure help explicitly says all
  libraries including `.oct` and `.mex` files must then compile with C++17.
- Octave 9 (2024) still has a mandatory C++11 fallback, while C++17 is required
  for a Qt 6 build and for the optional PMR configuration.
- Octave 10 (2025) makes C++17 mandatory.
- Octave 11 retains C++17 as its mandatory baseline.

These requirements can be inspected directly in the official GNU Octave
[`configure.ac` for Octave 8](https://github.com/gnu-octave/octave/blob/release-8-1-0/configure.ac#L303-L331),
[`configure.ac` for Octave 9](https://github.com/gnu-octave/octave/blob/release-9-4-0/configure.ac#L303-L376),
[`configure.ac` for Octave 10](https://github.com/gnu-octave/octave/blob/release-10-3-0/configure.ac#L303-L324),
and
[`configure.ac` for Octave 11](https://github.com/gnu-octave/octave/blob/release-11-3-0/configure.ac#L308-L324).
Octave's MEX compatibility layer is built through `mkoctfile --mex`; it uses the
compiler and flags recorded for that Octave installation. Consequently C++17 is
a safe baseline for Octave 10+, but Octave 8 and 9 installations cannot be
assumed to expose C++17 mode even though their underlying compiler often supports
it. Those installations require an explicit feature/compile test.

### Availability and MEX test on the review platform

The review environment is x86-64 Ubuntu 24.04 LTS (Noble). Its enabled native
Ubuntu `noble/universe` repository provides `octave` and `octave-dev` version
`8.4.0-1build5`; no Snap is needed. `octave-dev` supplies `mkoctfile` and the
headers needed to compile external `.oct` and `.mex` modules. The reproducible
installation command is:

```sh
apt-get update
apt-get install --no-install-recommends octave octave-dev
```

This native installation was exercised rather than merely inspected. With the
repository's GCC 13 toolchain, `mkoctfile --mex -std=c++17` successfully compiled
an external MEX source that calls `std::transform_reduce`. Octave 8.4 loaded the
resulting `/tmp/cxx17_mex.mex`, called it, and returned the expected dot product
of 32. This demonstrates that the non-Snap package available on this particular
platform can load external C++17 MEX files when the language flag is explicit.
It does not remove the need to test xfemm's full MEX build or other Octave 8/9
distributions, whose compilers and build options may differ.

For xfemm specifically, all five `mfemm/MMakefile_*.m` build descriptions now
request `-std=c++17`, including the environment passed to Octave. The obsolete
C++ language option was removed from `CFLAGS`, since a C++ dialect does not belong
in the C compiler flags. The fsolver makefile also relies on `mkoctfile` to supply
Octave's configured link directories and libraries instead of constructing an
incorrect `/usr/lib/octave/<version>` path. With those changes, a clean build of
all five MEX targets and the `mfemm_setup(..., 'RunTests', true)` functional tests
completed under the native Octave 8.4 package. MATLAB and other operating systems
remain outside this environment's tested matrix.

## Merge-blocking correctness and compatibility defects

1. **It requires the project to move to C++17.** The proposal uses
   `std::transform_reduce`, which is unavailable under the previous C++14
   configuration. A clean Release build of xfemm and its test suite succeeds
   after raising the required standard to C++17, so this repository now makes
   that requirement explicit. The MATLAB/Octave makefiles now request C++17 too,
   and the existing solver compiles and passes mfemm's tests with Octave 8.4.
2. **`Q` is never allocated.** `Create` previously allocated `Q`, and the heat and
   electrostatic solvers write permutation indices through it. The proposal
   leaves `int *Q` uninitialized, producing undefined behavior in those solvers.
3. **`PCGSolve(0)` no longer clears the initial solution.** The existing contract
   says the flag indicates whether a guess is present and explicitly zeroes `V`
   when it is false. The proposal removes that behavior. This can reuse stale
   values after a prior solve and change convergence or results.
4. **The new matrix routines assume every row has a first (diagonal) entry.**
   `MultA`, `MultPC`, and the singularity check use `M[i][0]`/`front()` without an
   emptiness check, while the proposed `Wipe` clears every row. Calling a matrix
   operation after `Wipe` but before every diagonal is reinserted becomes an
   out-of-bounds access.
5. **The patch mixes unrelated output changes into the solver rewrite.** Replacing
   the existing nonlinear-iteration reporting in planar and axisymmetric code is
   unnecessary for evaluating the storage optimization and complicates review.

The vector field ordering (`int` before `double`) probably reduces entry size on
common ABIs, but this should be verified rather than treated as portable fact.
Likewise, changing the convergence comparison from `sqrt(res/res_o) > Precision`
to `res/res_o > Precision*Precision` is algebraically reasonable only while the
preconditioned residual remains nonnegative and finite; regression tests should
cover breakdown and edge cases.

## Required validation for a replacement pull request

A narrower follow-up should preserve the public behavior and compile under the
project's required C++17 mode, then report, for baseline and candidate built with
identical Release flags:

- end-to-end meshing and solving time, plus solver-only time;
- representative small, medium, and large planar and axisymmetric models;
- node/element/nonzero counts and PCG/nonlinear iteration counts;
- at least 10 timed runs after warm-up, with median and dispersion;
- peak resident memory;
- solution comparisons (`.ans` quantities or solution vectors) at the configured
  tolerance;
- magnetic, electrostatic, and heat-flow regression tests, including repeated
  solves with and without an initial guess, periodic/antiperiodic boundaries, and
  a `Wipe`/reassembly cycle.

Useful implementation candidates are CSR after assembly, or sorted row vectors
with capacity reserved from mesh adjacency. Benchmark assembly separately:
binary-search insertion into a vector can shift a row repeatedly, so faster
iterations may otherwise hide slower construction.

## Reproduction notes

The review downloaded the issue's linked patch, applied it to a detached clean
worktree at the current repository revision, configured a Release build, and ran
the normal CMake build. The proposal failed with the former C++14 setting because
`std::transform_reduce` was unavailable. Raising the project requirement to
C++17 allowed both the unmodified project and the proposal to compile, confirming
that the language-version bump itself is viable. Static inspection still found
the `Q`, initial-guess, and empty-row problems above. No speed ratio is reported
because those correctness defects prevent a trustworthy comparison without
silently repairing the submitted implementation first.
