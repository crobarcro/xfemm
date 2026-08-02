# Python `FemmSession` interface plan

## Scope and compatibility target

The first Python API will be a single stateful magnetic-analysis class equivalent
to MATLAB's new `xfemm.femmsession`. It will not reproduce the older MATLAB
problem, property, or post-processor class hierarchy. The compatibility target
is the public method surface of `xfemm.femmsession`, including all of its useful
analysis methods. Post-processing is part of this one class, not a separately
constructed or user-visible `fpproc` object.

The proposed import and construction syntax is:

```python
from xfemm import FemmSession

with FemmSession("motor.fem") as session:
    session.setBackend("tangle")
    session.setCircuit("phase-a", "current", 10 + 0j)
    session.setAGEPosition("airgap", rotor_angle, 0.0)
    status = session.solve()
    b = session.getb([[0.0, 0.0], [1.0, 0.0]])
    session.accept()
```

`FemmSession` is the sole public Python class name. It follows Python class-naming
conventions and deliberately omits the `m` from mfemm because it is not
MATLAB-specific. MATLAB keeps the idiomatic package-qualified spelling
`xfemm.femmsession`. Existing MATLAB
camel-case method spellings are intentional. Adding snake-case aliases is not
part of the first release.

The initial package supports magnetic `.fem`/`.ans` files only. Electrostatic
and heat-flow sessions, model construction, Lua command emulation, and the
legacy multi-class API are explicitly out of scope.

The MATLAB class moves into its own package at `+xfemm/femmsession.m`, so its
fully qualified spelling is also `xfemm.femmsession(...)`. This is an intentional
breaking rename from the short-lived `mfemmsession`: the old class and alias are
not retained because the interface has not yet been released, and retaining the
MATLAB-derived `m` would give the wrong meaning to the shared API.

## Packaging and implementation architecture

Use a compiled CPython extension (preferably pybind11) plus a small pure-Python
facade:

```text
python/
  pyproject.toml             build metadata and NumPy dependency
  src/xfemm/__init__.py      public exports
  src/xfemm/session.py       validation, conveniences, persistence, context manager
  src/xfemm/_xfemm.*         compiled native extension
  tests/                     pytest compatibility and numerical tests
```

The native layer must call the same `AnalysisSession`, mesher, solver, and
`FPProc` C++ APIs as the MEX wrappers. Python and MATLAB are conversion and
lifetime-management layers only: field interpolation, smoothing, selection,
contours, integrals, mesh-derived values, circuits, and AGE calculations must
not be independently implemented in either wrapper. Before binding them, move the
reusable `SessionGateway` logic out of `session_interface_mex.cpp` into a
MATLAB-independent C++ class. Both bindings should use that class, preventing
the Python and MATLAB session semantics from drifting.

In-memory post-processing is a version-1 requirement, not a later optimization.
`solve()` must construct an immutable native solution snapshot directly from
the solved model, mesh, nodal values, circuits, and AGE data and install a native
post-processor view over that snapshot. Every point, contour, block, mesh,
circuit, and air-gap query is then immediately available on the same
`FemmSession` instance. The normal solve/query path must neither write nor parse
an `.ans` file, and it must work in a directory where the process has no write
permission.

The existing `MagneticSolutionSnapshot`/`FPProc` work demonstrates the desired
ownership model, but the solver still needs a file-format-independent transfer
object. Introduce a neutral `MagneticSolutionData` (name provisional) in a
common native library. It should own the problem description, mesh topology,
nodal solution, circuit results, and periodic/AGE coupling needed by
post-processing. `FSolverAnalysisBackend` populates it directly, and `FPProc`
constructs its immutable snapshot from it. The `.ans` reader and writer become
optional adapters to and from the same data rather than the bridge between the
solver and post-processor.

Build wheels with `scikit-build-core` and CMake so the extension shares xfemm's
existing build definitions. Start with CPython 3.10+ and NumPy 1.23+ on 64-bit
Linux, macOS, and Windows. Source builds remain supported; wheels should include
the required native code rather than require MATLAB, Octave, or a separately
installed xfemm.

Native session instances are not thread-safe. Different instances may be used
concurrently only after tests demonstrate that the underlying solver has no
process-global mutable state. Until then, document all solve/post-processing
calls as serialized. Release the GIL around long native meshing and solving
operations only when that condition is satisfied.

### Initial implementation status

The first implementation slice now lives under `python/` and exports only the
idiomatic `FemmSession` name. It binds model loading, Triangle/Tangle selection,
meshing, circuit/AGE/frequency/time setters, solve/result, accept/reject, and
lifecycle management directly to the existing session C++ APIs. Tests verify
that solving creates no `.ans` file.

The native in-memory bridge now initializes `FPProc` from the solved `FSolver`
state and reuses the same finalization path as `.ans` loading. Python point
values, smoothing, contour/line integrals, block selection/integrals, problem
and circuit information, and mesh counts are thin conversions around `FPProc`.
No `.ans` file is created in this path. The earlier independently implemented
planar interpolation was removed. Remaining mesh/group tables and air-gap
conveniences must likewise be added only as thin `FPProc` bindings. Explicit
`.ans` and accepted-state persistence are still pending.

The MATLAB `xfemm.femmsession` wrapper also initializes `FPProc` through
`FPProc::OpenDocument(problem, solver, solution)` in the session MEX gateway.
Its inherited post-processing methods dispatch to that session-owned processor,
so repeated solves do not create or parse temporary `.ans` files. Explicit
solution export remains the only session operation that writes an `.ans` file.

## Python data conventions

These conventions resolve MATLAB/Python differences once, rather than leaving
them to individual methods:

* Numeric array results are `numpy.ndarray`, never nested lists. Real output is
  `float64`, complex output is `complex128`, integer counts/indices are `int64`,
  and masks are `bool_`.
* Scalar mathematical results are NumPy scalars of the corresponding dtype.
  Methods that conceptually return a count (`mesh`, `numelements`, and similar)
  return Python `int`.
* Inputs accept Python scalars and array-like objects. Coordinate inputs are
  converted to contiguous `float64` arrays; complex circuit values accept real
  or complex scalar numbers.
* Point methods accept either `(x, y)` with equal-size broadcastable arrays or a
  single array whose last dimension is two. They flatten points in C order and
  return one row per point. This is deliberately idiomatic Python rather than a
  literal copy of MATLAB's column-major flattening.
* Every batch-shaped result uses samples on axis 0. Thus `getpointvalues` is
  `(n, 14)`, `getb`/`geth` are `(n, 2)`, `geta` is `(n,)`, and mesh tables retain
  their documented `(n, columns)` shape. Empty results have the same rank, such
  as `(0, 2)`, rather than collapsing to `(0,)`.
* Public mesh element and node indices remain **one-based**, matching MATLAB
  return values and accepted arguments. The wrapper converts only at the native
  boundary. This avoids silently changing the values in `getelements` and makes
  cross-language comparisons straightforward.
* Multiple MATLAB outputs become a Python tuple. A MATLAB struct becomes a
  dictionary with exactly the MATLAB field names. No output-count (`nargout`)
  behavior is emulated.
* Paths accept `str` and `os.PathLike`; UTF-8/path conversion is handled by the
  native boundary. Names accept `str`. Booleans are not accepted where a numeric
  index is expected.
* Native failures become specific exceptions rooted at `XfemmError`:
  `FileFormatError`, `SessionStateError`, `UnknownCircuitError`,
  `UnknownAirGapError`, and `SolverError`. Argument type/shape errors use normal
  `TypeError` or `ValueError`.

The row-oriented point-result decision is the largest intentional difference
from MATLAB, where `getpointvalues`, `getb`, and `geth` place quantities in rows
and points in columns. Returning the MATLAB orientation would technically be
closer but would conflict with the existing mesh methods and normal NumPy usage.
It must be called out prominently in migration documentation.

## Session lifecycle and state methods

| Method | Python signature | Return contract |
|---|---|---|
| constructor | `FemmSession(filename)` | New open session; raises on parse failure. |
| `close` | `close()` | `None`; idempotently releases native resources. |
| context manager | `__enter__()`, `__exit__(...)` | Returns self and always calls `close`. |
| `setBackend` | `setBackend(name)` | `None`; accepts case-insensitive `"triangle"` or `"tangle"`. |
| `mesh` | `mesh()` | Python `int` element count. |
| `setCircuit` | `setCircuit(name, constraint, value=0)` | `None`; constraint is `current`, `voltage`, `open`, or `coupled`. For the latter two, a supplied value is ignored for MATLAB compatibility. |
| `setAGEPosition` | `setAGEPosition(name, innerAngle, outerAngle)` | `None`; angles are scalar degrees. |
| `setFrequency` | `setFrequency(frequency)` | `None`; scalar hertz. |
| `setTime` | `setTime(time)` | `None`; scalar user time. |
| `solve` | `solve()` | Dictionary described below; refreshes the post-processor on success. |
| `result` | `result()` | Dictionary described below; arrays are independent read-only copies. |
| `accept` | `accept()` | `None`; latest trial becomes the initial state. |
| `reject` | `reject()` | `None`; discards latest trial and invalidates post-processing access to it. |
| `saveState` | `saveState(filename)` | `None`; writes a versioned, non-pickle archive. |
| `loadState` | `loadState(filename)` | `None`; validates schema and model identity before mutation. |
| `saveSolution` | `saveSolution(filename)` | `None`; explicitly persists the current trial as a legacy `.ans` file. |
| `loadSolution` | `loadSolution(filename)` | `None`; validates and installs a `.ans` solution matching the session model, without solving. |

`solve()` returns exactly these keys: `success` (`bool`), `id` (`int`), `time`
(`float`), `nodeCount`, `elementCount`, `meshGenerationCount`, `solveCount`,
`operatorAssemblyCount`, and `rightHandSideAssemblyCount` (all `int`). Solver
failure raises `SolverError`; `success` therefore reports successful completion
for MATLAB parity rather than replacing exceptions.

`result()` returns `id`, `time`, `current`, `fluxLinkage`, `terminalVoltage`,
`A`, `x`, and `y`. The five vector fields are one-dimensional arrays.
`current`, `fluxLinkage`, and `terminalVoltage` are `complex128`; unavailable
terminal voltage is complex NaN. `A` is `float64` for static solutions and must
become `complex128` when harmonic nodal results are exposed. `x` and `y` are
`float64`.

Unlike MATLAB MAT-files, Python state files should be a versioned ZIP/NPZ
container containing JSON metadata and opaque backend bytes. They must include
the state schema version, xfemm ABI/version, a SHA-256 identity of the parsed
model, `modelRevision`, `parameterRevision`, and time. `allow_pickle=False` is
mandatory when loading. Loading a mismatched model or ABI fails before changing
the current session. MATLAB/Python state-file interoperability is not promised
in version 1; the accepted-state semantics are interoperable, the serialization
containers are not.

Legacy solution I/O remains useful but is always explicit. `saveSolution()` is
invalid before a solve or load. `loadSolution()` imports through the `.ans`
adapter into the same immutable snapshot used by in-memory solves, validates
that its problem identity matches the session's loaded `.fem` model, and then
makes all query methods available. Neither method is called implicitly by
`solve()`. A future `fromSolution()` read-only constructor could be added if
users need to inspect an `.ans` file without its matching model, but it is not
required for version 1.

## Post-processing method surface

All methods below require a successful, unrejected trial. They retain the exact
MATLAB names.

### Point and circuit queries

| Method | Return |
|---|---|
| `getpointvalues(x, y=None)` | `(n, 14)` array, columns in the MATLAB-documented order `A, B1, B2, Sig, E, H1, H2, Je, Js, Mu1, Mu2, Pe, Ph, ff`. |
| `getb(x, y=None)` | `(n, 2)` array `[B1, B2]`. |
| `geth(x, y=None)` | `(n, 2)` array `[H1, H2]`. |
| `geta(x, y=None)` | `(n,)` array. |
| `getprobleminfo()` | Tuple `(problem_type, frequency, depth, length_unit)` with native numeric meanings preserved. |
| `getcircuitprops(circuitname)` | `(3,)` complex array `[current, voltage, flux_linkage]`. |
| `circuitRL(circuitname)` | Tuple `(R, L)` of complex NumPy scalars, preserving MATLAB division behavior (including NaN/Inf). |

### Contours, blocks, and integration

| Method | Return |
|---|---|
| `smoothon()`, `smoothoff()`, `smooth(flag)` | `None`; `flag` is case-sensitive `"on"` or `"off"` for parity. |
| `clearcontour()`, `addcontour(x, y)`, `newcontour(x, y)` | `None`; coordinate arrays must have the same flattened length, and a new contour requires at least two points. |
| `lineintegral(type)` | Always a tuple, with tuple length and scalar meanings matching the MATLAB integral type. |
| `selectblock(x, y, clearselected=False)` | `None`; scalar coordinate pair. |
| `groupselectblock(groupno, clearselected=False)` | `None`; scalar or array-like group numbers; an empty input selects every block. |
| `selectallblocks()`, `clearblock()` | `None`. |
| `blockintegral(type, x=None, y=None)` | Numeric scalar; optional coordinates clear and replace the selection first. |
| `gapintegral(boundname, type)` | Numeric scalar or tuple exactly matching the native integral type. |
| `totalfieldenergy()`, `totalfieldcoenergy()` | Numeric scalar. |

Integral type constants should additionally be offered through `IntEnum`
classes, while integers remain accepted for compatibility. The implementation
must create a golden table from the MEX behavior for every supported line,
block, and gap integral type; current MATLAB comments are incomplete and cannot
serve as the sole contract.

### Mesh and air-gap queries

| Method | Return |
|---|---|
| `nummeshnodes()`, `numelements()` | Python `int`. |
| `getvertices(n=None)` | `(n, 6)` `float64`; default is all elements. |
| `getelements(n=None)` | `(n, 7)` numeric table; node columns contain one-based integer-valued entries. A homogeneous NumPy array necessarily represents the whole table as `float64`. |
| `getcentroids(n=None)` | `(n, 2)` `float64`. |
| `getareas(n=None)`, `getvolumes(n=None)` | `(n,)` `float64`. |
| `numgroupelements(groupno)` | Array matching the input's shape with `int64` dtype; scalar input returns Python `int`. |
| `getgroupvertices(groupno)` | `(n, 6)` `float64`, concatenated in requested group order. |
| `getgroupelements(groupno)` | `(n, 7)` table, concatenated in requested group order. |
| `getgroupcentroids(groupno)` | `(n, 2)` `float64`. |
| `getgroupareas(groupno)`, `getgroupvolumes(groupno)` | `(n,)` `float64`. |
| `getgapb(bound_name, angles)` | `(n, 2)` complex array `[Br, Bt]`. |
| `getgapa(bound_name, angles)` | `(n, 2)` complex array preserving the two native outputs until their physical meanings are verified. |
| `numgapharmonics(bound_name)` | Python `int`. |
| `getgapharmonics(bound_name, n)` | Tuple of six `(k,)` complex arrays: `(acc, acs, brc, brs, btc, bts)`. |

The MATLAB `getvolumes` implementation appears to dispatch to `getareas`; the
Python implementation should call the native volume operation and a regression
test should establish whether MATLAB has a wrapper bug. Compatibility means
matching the documented physical result, not preserving an accidental dispatch
typo.

### Plotting

`plotBfield`, `plotHfield`, and `plotAfield` are public inherited methods, but
they are visualization conveniences rather than solver return APIs. Implement
them in a second milestone using Matplotlib and return `(figure, axes)`. Their
field sampling will use the query methods above. Absence of Matplotlib must not
prevent importing or using the solver; calling a plot method without the
optional `plot` dependency raises a clear `ImportError`.

## State machine and behavioral rules

The facade should enforce these observable states:

1. **Open/model loaded:** parameter setters and `mesh()` are valid; result and
   post-processing methods raise `SessionStateError`.
2. **Trial solved:** `result()` and post-processing are valid. A new setter or
   backend change marks that trial stale for subsequent solving but does not
   mutate its stored values.
3. **Accepted:** the trial remains queryable and is installed as the next
   initial state. Calling `accept()` without a trial is an error.
4. **Rejected:** the latest trial and its post-processor are invalidated.
5. **Closed:** every method except idempotent `close()` raises
   `SessionStateError`.

The object owns its native resources; copying and pickling the live object are
disabled. It may be moved only through ordinary Python reference assignment.
`saveState` persists accepted numerical state, not the live object.

## Weaknesses and ambiguities discovered

1. **“Equivalent NumPy arrays” does not define orientation.** MATLAB point
   results use quantities-by-points while Python convention is
   samples-by-features. This plan chooses row-oriented Python results and treats
   values, not memory layout, as the equivalence criterion.
2. **The MATLAB surface includes inherited methods.** A Python implementation
   could accidentally expose only the session gateway methods and force users
   into a second post-processor object. This plan instead puts every
   computational session and post-processing method on `FemmSession` itself and
   defers only plotting to milestone 2.
3. **The current MATLAB refresh is file-backed.** Reusing it would make hidden
   file I/O part of every Python solve and undermine optimization-loop use. The
   first release therefore requires a native solution-data transfer into the
   post-processor. This is more implementation work, but the repository's
   immutable snapshot support provides a foundation and `.ans` parity fixtures
   can verify the new adapter.
4. **MATLAB output arity is dynamic.** Python cannot reproduce `nargout`.
   Tuple-returning methods need one stable full tuple, particularly
   `lineintegral`, `gapintegral`, and `getprobleminfo`.
5. **Index bases naturally differ.** This plan favors strict cross-language
   compatibility and retains one-based public indices, despite NumPy's normal
   zero-based convention.
6. **Persistence is underspecified and unsafe if implemented with pickle.** A
   versioned, validated non-pickle format is required. Cross-language MAT-file
   compatibility is deferred.
7. **Complex harmonic behavior needs fixtures.** The MEX trial export currently
   exposes only real nodal `A`, while circuit and air-gap values can be complex.
   AC-model tests must define dtypes, NaN representation, and phase conventions.
8. **Some MATLAB documentation/code is inconsistent.** `getgapa` documents flux
   density, integral comments are incomplete, and `getvolumes` may call the area
   command. Native behavior plus numerical fixtures must resolve these rather
   than copying comments blindly.
9. **Thread safety and ownership have not been established.** pybind11 makes it
   easy to release the GIL, but doing so before auditing solver globals could
   introduce races.
10. **Distribution cost is significant.** Wheels embed an old Triangle/Lua and
    solver stack with multiple licenses. Wheel metadata and source notices need
    a license audit before publishing.

## Questions to settle before implementation

The following choices are proposed defaults, but should be explicitly approved:

1. **Point-array orientation:** approve Pythonic `(n, fields)` output rather
   than literal MATLAB `(fields, n)` output?
2. **Indexing:** retain one-based indices for compatibility (proposed), or make
   all Python indices zero-based and provide an opt-in compatibility mode?
3. **State files:** is semantic compatibility sufficient, or must Python read
   and write MATLAB `saveState` MAT-files in version 1?
4. **Plotting:** is deferring the three plotting helpers acceptable for the
   first usable release?
5. **Platform baseline:** are CPython 3.10+, NumPy 1.23+, and 64-bit desktop
   platforms acceptable?

## Implementation milestones and acceptance tests

1. **Contract fixtures:** run MATLAB/Octave against representative static,
   harmonic, planar, axisymmetric, periodic, and AGE problems. Record method
   shapes, dtypes, complex values, errors, and all integral arities.
2. **Shared native gateway:** extract session ownership/operations from the MEX
   file without changing MATLAB behavior; run existing C++ and MATLAB session
   tests.
3. **Direct solution snapshot:** introduce the neutral native solution data,
   populate it from `FSolverAnalysisBackend`, construct `FPProc` without `.ans`
   serialization, and prove parity against a persisted round trip. This C++
   bridge is a prerequisite for exposing any Python post-processing method.
4. **Python build and single-class lifecycle:** add packaging, native
   construction, cleanup, setters, mesh, solve/result, accept/reject, context
   management, and direct ownership of the snapshot-backed post-processor.
5. **Computational post-processing:** bind the existing `FPProc` point, contour,
   block, circuit, mesh, group, and air-gap methods on `FemmSession`; the binding
   contains only centralized NumPy input/output conversions.
6. **Persistence:** add schema validation, model/ABI identity checking, corrupt
   file tests, and atomic writes.
7. **Parity suite:** compare Python outputs element-by-element with stored
   MATLAB fixtures. Use exact checks for topology/counts and documented
   tolerances for floating/complex results. Include empty, scalar, batched,
   non-contiguous, invalid-name, invalid-index, pre-solve, rejected, and closed
   cases.
8. **Distribution and documentation:** build and install wheels in clean
   environments, audit bundled licenses, document the intentional orientation
   and indexing decisions, and add a MATLAB-to-Python migration example.
9. **Optional plotting and performance:** add the Matplotlib extra and benchmark
   solve/query loops. Add a test that traces filesystem operations and proves a
   normal in-memory solve and post-processing workflow creates no solution file.

The first release is complete when the computational surface above is present,
all arrays follow the declared shapes/dtypes, MATLAB fixture comparisons pass,
resources survive repeated create/solve/reject/close cycles without leaks, and
the package installs from both a wheel and an isolated source build without
MATLAB or Octave.
