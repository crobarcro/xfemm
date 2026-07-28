# High-level magnetic analysis sessions

`AnalysisSession` is the authority for a user-visible magnetic analysis. It
separates three kinds of data which legacy solver structures combine:

* `ModelDefinition` owns geometry, labels, materials, boundary and magnetic
  circuit-port definitions, and structural analysis settings.
* `SolveParameters` owns values that can vary between evaluations: circuit
  constraints, frequency, time, air-gap position, and an optional accepted
  state.
* `PreparedAnalysis` is disposable solver input derived from the other two.

The session also owns the solver mesh. A Triangle backend is selected by
default, and applications can inject another `MesherBackend` with `setMesher()`.
Backend-neutral `MeshingOptions`, the diagnostics from the most recent attempt,
and an immutable `shared_ptr<const SolverMesh>` remain available on the
session. `ensureMesh()` creates that mesh on demand. Each successful replacement
gets a new topology identity so solver backends can distinguish topology from
the changing operator and right-hand side.

The session exposes the model as const data. Callers change it through named
operations such as `setMaterialProperty`; circuit inputs, frequency, and AGE
positions likewise have named operations. Each operation invalidates the
smallest reusable layer that can safely be retained. `solve()` always calls
`synchronize()` before invoking its backend, while explicit synchronization is
available to separate preparation failures and timing from the solve itself.

## Magnetic circuit ports and future transient coupling

A model circuit is an electromagnetic port, not an instantaneous source and
not an entire external electrical network. Per-evaluation constraints can be
prescribed current, prescribed voltage, open circuit, or a coupled unknown.
`setCircuitCurrent` remains a convenience API; backends receive the complete
set of port constraints so an external or future XFEMM-managed circuit solver
can perform coupled iterations.

Solutions are trials. Calling `solve()` does not advance time history.
`acceptSolution()` explicitly converts a current trial into an immutable
`AcceptedState`, and rejects trials made obsolete by subsequent model or input
changes. This permits several current/position guesses at one time point and
supports step rejection by a future time integrator. Flux linkage is a
first-class port result; transient induced voltage must be calculated by the
coupled integration scheme rather than by applying the harmonic `2*pi*f` rule.

## Derived transformations

Series winding expansion is prepared state. The session creates one
`PreparedCircuit` mapping per winding label, including its signed turn count,
but never rewrites model circuit definitions or label circuit references.
Changing its source circuit, labels, turns, or constraint recreates that
mapping. Material slope data and positioned air-gap coupling are handled the
same way: they are recreated from authoritative inputs and are never copied
back into `FemmProblem`.

Concrete solvers implement `AnalysisSolverBackend`. The boundary deliberately
uses const model, parameter, and prepared views and receives the session's
immutable mesh plus its topology identity during synchronization. Circuit,
frequency, time, material, initial-state, and sliding-position updates preserve
the mesh; only mesher selection, meshing-option changes, and future geometry or
other mesh-affecting model operations invalidate it. A backend can retain
private matrix and nonlinear work state, but cannot turn those objects into
user-visible authority.
