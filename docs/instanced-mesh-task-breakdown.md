# Instanced mesh task breakdown

This is the working checklist for the
[instanced mesh implementation plan](instanced-mesh-implementation-plan.md).
Tasks are ordered by dependency, not merely grouped by topic. Complete a task only
when its deliverables and checks are in the same reviewed change; do not mark an
epic complete until all of its exit checks pass.

## Tracker conventions

- Status uses Markdown checkboxes: `[ ]` ready/not started, `[x]` complete.
- A task ID is stable and should appear in its branch name, commit/PR description,
  and new test names where practical.
- “Depends on” names hard prerequisites. Tasks at the same dependency level may
  proceed in parallel once their shared prerequisites are merged.
- Triangle is the differential/compatibility backend. New instancing behavior is
  implemented against Tangle unless a task explicitly says otherwise.
- Every behavior change needs a test that fails before the change and passes after
  it. Generated fixtures must record their generator inputs and external revision.

## Milestone A — A real Tangle backend

Milestone result: `TangleMesherBackend` executes pinned Tangle code in memory,
returns a validated `SolverMesh`, and is trustworthy enough to become the default.

### A1 — Define and pin the Tangle dependency

- [x] **A1.1: Choose the upstream integration contract.**
  - Depends on: none.
  - Determine from `dcm3c/tangle` whether xfemm can consume a CMake target directly
    or must temporarily provide a thin target around the pinned upstream sources.
  - Require a callable library entry point with no CLI `main`, no global process
    exit, and no mesh-file round trip.
  - Deliverable: a short decision record in the implementation PR explaining the
    selected target/API and the rejected alternatives.
- [x] **A1.2: Add reproducible CMake dependency resolution.**
  - Depends on: A1.1.
  - Add a pinned Tangle revision and support both a configured local source path
    (offline/development) and retrieval of that exact revision.
  - Do not silently fall back to a different Tangle revision or to Triangle when
    Tangle was explicitly requested.
  - Checks: clean online configure, clean offline/local-source configure, and a
    diagnostic for an unavailable explicitly requested provider.
- [x] **A1.3: Install and document Tangle licensing.**
  - Depends on: A1.2.
  - Include the pinned MIT license in source/binary packaging and identify Tangle
    in dependency documentation.
  - Check the generated install/package manifest contains the license.

### A2 — Convert real Tangle output

- [x] **A2.1: Add a pure Tangle-to-`SolverMesh` converter.**
  - Depends on: A1.2.
  - Convert nodes, triangle connectivity/region attributes, edges/markers, PBC
    pairs, and AGE definitions without accessing `FemmProblem` or the filesystem.
  - Validate signed/unsigned and zero-/one-based conversions at the boundary.
  - Checks: focused unit tests for ordinary, periodic-only, and AGE-bearing
    synthetic Tangle meshes, including invalid indices.
- [x] **A2.2: Execute Tangle from `TangleMesherBackend`.**
  - Depends on: A2.1 and the required upstream in-memory input API from A1.1.
  - Remove the current call to `TriangleMesherBackend` and invoke the real engine.
  - Map Tangle status values to `MeshStatus` and preserve actionable diagnostics.
  - Checks: an injected engine counter proves Tangle executed once; a scratch
    directory remains empty; diagnostics name `Tangle` only when Tangle ran.
- [ ] **A2.3: Exercise backend options.**
  - Depends on: A2.2.
  - Map all supported `MeshingOptions` explicitly. Reject or diagnose unsupported
    options rather than ignoring them.
  - Checks: minimum angle, area/size controls, boundary Steiner suppression,
    unused-vertex suppression, verbosity, and invalid values.

### A3 — Establish differential confidence

- [x] **A3.1: Build a backend-neutral mesh invariant checker.**
  - Depends on: A2.1.
  - Check finite coordinates, valid connectivity, positive element area,
    boundary-edge ownership, PBC indices/types, and every AGE node reference.
  - Reuse it for both Tangle and Triangle tests.
- [x] **A3.2: Add ordinary and periodic differential fixtures.**
  - Depends on: A2.2 and A3.1.
  - Run both engines on checked-in non-periodic, periodic, and antiperiodic FEMM
    problems. Do not require identical node/element ordering or triangulation.
  - Compare exact invariants where applicable and solver fields/energies within
    fixture-specific tolerances.
- [ ] **A3.3: Add AGE differential fixtures.**
  - Depends on: A3.2.
  - Compare AGE ring sizes, radii, centers, periodicity, valid quadrature, and a
    small angle sweep through the existing solver/session path.
  - Include a periodic `.pbc` fixture with zero AGE records to prevent conflation.
- [ ] **A3.4: Make Tangle the default backend.**
  - Depends on: A2.3, A3.2, and A3.3.
  - Change `AnalysisSession` default construction to `TangleMesherBackend` while
    preserving explicit Triangle injection.
  - Checks: default-selection test, explicit-backend tests, and full C++/CLI/MATLAB
    regression suites. Record any platform-specific extended suite separately.

**Milestone A exit check:** no Tangle-labeled path executes Triangle; Tangle is the
tested default; Triangle remains explicitly usable; neither backend creates files
through the in-memory backend API.

## Milestone B — Boundary matching as topology, not physics

Milestone result: callers can ask Tangle for matching ordered boundary chains
without implicitly imposing periodic or antiperiodic field constraints.

- [ ] **B1: Add `MeshingRequest` and compatibility adapter.**
  - Depends on: A3.4.
  - Add request options, boundary matches, and an explicit “emit field constraint”
    choice. Keep the old `(problem, periodic, options)` entry point as an adapter.
  - Checks: old and new calls produce equivalent `SolverMesh` results.
- [ ] **B2: Specify stable geometry references.**
  - Depends on: B1.
  - Select stable identifiers for matched source segments/arcs and reject stale,
    mixed line/arc, duplicated, or missing references before invoking Tangle.
  - Checks: geometry mutation and invalid-reference cases produce diagnostics.
- [ ] **B3: Expose Tangle topology-only paired refinement.**
  - Depends on: B1 and B2; expected upstream Tangle change.
  - Synchronize splitting and return ordered chain correspondence independently
    of PBC output semantics.
  - Checks: straight and arc chains, forward and reverse order, refinement caused
    from either side, and unequal-chain rejection.
- [ ] **B4: Convert matches to backend-neutral seam data.**
  - Depends on: B3.
  - Return ordered node references and orientation without adding them to
    `SolverMesh::periodicConstraints` unless explicitly requested.
  - Checks: the same matched geometry once with and once without field constraints.
- [ ] **B5: Characterize Triangle compatibility.**
  - Depends on: B4.
  - Keep legacy Triangle periodic behavior intact. If topology-only matching is
    unsupported there, return a clear capability diagnostic rather than emulating
    it with periodic physics.

**Milestone B exit check:** a Tangle mesh can contain two conforming matched chains
and zero periodic constraints; requesting periodic or antiperiodic physics emits
the expected typed node pairs.

## Milestone C — Instanced value model and materializer

Milestone result: hand-authored templates can be expanded deterministically into
an ordinary solver-compatible mesh without invoking a mesher.

- [ ] **C1: Add transform and template value types.**
  - Depends on: none; merge after B1 to avoid competing mesh API changes.
  - Add `RigidTransform2D`, `TemplateSeam`, `MeshTemplate`, `SeamConnection`,
    `MeshInstance`, and `InstancedMesh` under `cfemm/libfemm/mesh`.
  - Initially reject reflections and non-finite transforms.
- [ ] **C2: Validate template-local topology.**
  - Depends on: C1.
  - Validate local node/edge/element indices, unique seam nodes, ordered connected
    chains, transform orientation, and compatible seam cardinality.
  - Checks: one failure test per invariant with stable diagnostic categories.
- [ ] **C3: Build deterministic node provenance and seam welding.**
  - Depends on: C2.
  - Map `(template, instance, local node)` to global nodes. Weld only declared
    seam pairs; use transformed coordinates only to validate a declaration.
  - Checks: translated pair, reversed seam, closed annular ring, and close but
    undeclared nodes that must remain distinct.
- [ ] **C4: Materialize elements and edges.**
  - Depends on: C3.
  - Remap connectivity, reject degenerate/reversed elements, canonicalize duplicate
    internal seam edges, and define marker-conflict behavior.
  - Checks: exact counts and positive signed areas for repeated wedges.
- [ ] **C5: Remap periodic and AGE topology.**
  - Depends on: C3.
  - Remap PBC pairs, AGE rings, quadrature nodes, and AGE node-index lists through
    the checked provenance table while preserving periodicity and weights.
  - Checks: synthetic PBC-only and AGE templates plus invalid local references.
- [ ] **C6: Add reverse provenance and deterministic hashing.**
  - Depends on: C4 and C5.
  - Provide global-to-template/instance maps and stable topology/layout identities.
  - Checks: repeat runs hash identically; transform or seam changes alter the
    appropriate identity; physics-only metadata does not alter topology identity.
- [ ] **C7: Add sanitizer test configuration.**
  - Depends on: C2–C6.
  - Run the materializer unit suite with address and undefined-behavior sanitizers
    on a supported CI platform.

**Milestone C exit check:** a hand-authored annular wedge materializes to a closed,
valid `SolverMesh` with deterministic forward/reverse maps and correctly remapped
PBC/AGE data.

## Milestone D — Tangle-generated rotational templates

Milestone result: xfemm meshes one selected tile with Tangle, repeats it, and
solves the materialized mesh through the unchanged `FSolver`.

- [ ] **D1: Define `TemplateRequest`.**
  - Depends on: B4 and C1.
  - Identify source geometry, center of rotation, instance transforms, named seam
    chains, and supported region policy. Validate complete rotational coverage.
- [ ] **D2: Extract one tile PSLG for Tangle.**
  - Depends on: D1.
  - Preserve source markers, block-region attributes, holes, mesh controls, units,
    and stable mappings back to source entities.
  - Checks: isolated tile PSLG matches expected geometry and region seeds.
- [ ] **D3: Capture Tangle seam nodes and construct `MeshTemplate`.**
  - Depends on: D2 and B4.
  - Convert ordered boundary matches into template seams and prove Tangle is called
    exactly once per template.
- [ ] **D4: Materialize through `MeshResult`.**
  - Depends on: D3 and C6.
  - Construct instances, weld declared neighbor seams, materialize, and return the
    ordinary mesh plus additive provenance without breaking existing consumers.
- [ ] **D5: Add input diagnostics.**
  - Depends on: D4.
  - Cover incomplete rotations, overlap, seam mismatch, conflicting markers,
    anisotropic materials, reflections, and AGE topology inside a template.
- [ ] **D6: Solve a repeated-ring control.**
  - Depends on: D4 and D5.
  - Check exact topology invariants and compare sampled field, energy, and force
    values with a conventionally generated full ring at documented tolerances.

**Milestone D exit check:** a rotational template is triangulated once, expanded
deterministically, and solved by the unmodified solver with control-model agreement.

## Milestone E — Session and per-instance physics

Milestone result: an analysis session owns the instanced representation and can
change sources or AGE position without remeshing/materializing topology.

- [ ] **E1: Cache canonical instanced and materialized meshes.**
  - Depends on: D4.
  - Add separate template-topology, instance-layout, and materialized-topology
    identities to `AnalysisSession`.
- [ ] **E2: Implement cache invalidation rules.**
  - Depends on: E1.
  - Test geometry, transform, seam, material, circuit/current, initial state, and
    AGE-angle changes against meshing/materialization/topology-import counters.
- [ ] **E3: Resolve instance region overrides.**
  - Depends on: E1 and C6.
  - Apply circuit assignment, turns/current scale, and magnetization rotation while
    building prepared analysis; never mutate the shared template.
- [ ] **E4: Test independent field DOFs and physics.**
  - Depends on: E3.
  - Compare opposite coil sides and alternating magnet directions with explicitly
    drawn controls. Confirm instancing alone introduces no periodic constraint.
- [ ] **E5: Add C++ authoring/query API.**
  - Depends on: E2–E4.
  - Expose template creation, instance transforms/overrides, provenance queries,
    and capability diagnostics with API documentation.
- [ ] **E6: Extend post-processing provenance.**
  - Depends on: E5.
  - Allow selection and result attribution by template/instance while the legacy
    post-processor continues to consume global elements.
- [ ] **E7: Add MATLAB/MEX bindings.**
  - Depends on: E5 and E6.
  - Preserve all existing entry points and add one documented instanced smoke case.

**Milestone E exit check:** current and AGE angle sweeps reuse topology; transform
changes invalidate only required caches; instance-specific sources and results are
observable through C++ and MATLAB.

## Milestone F — RNFoundry machine validation

Milestone result: independent stator and rotor templates reproduce the existing
RNFoundry-derived machine results while reusing mesh topology across positions.

- [ ] **F1: Freeze generator provenance and tolerances.**
  - Depends on: E5.
  - Record the RNFoundry commit, design/options, xfemm commit, generated files,
    result schema, and justified tolerances.
- [ ] **F2: Generate comparable machine cases.**
  - Depends on: F1.
  - Produce conventional redraw, existing AGE/sliding, and instanced variants of
    the checked-in 12-pole, 36-slot design.
- [ ] **F3: Build independent stator and rotor templates.**
  - Depends on: F2.
  - Use a stator slot template and rotor pole/pole-pair template with independent
    counts; connect the domains only through AGE rings.
- [ ] **F4: Add one-position CI smoke comparison.**
  - Depends on: F3.
  - Compare winding flux linkage, coil flux-density magnitude, torque, topology
    invariants, number of Tangle calls, and solver topology imports.
- [ ] **F5: Add extended rotor-position sweep.**
  - Depends on: F4.
  - Compare redraw, sliding, and instanced results over all fixture positions;
    record solve time, meshing calls, node/element counts, and peak memory.
- [ ] **F6: Make fixture regeneration reproducible but optional.**
  - Depends on: F5.
  - Normal CI consumes checked-in fixtures without RNFoundry. Document the explicit
    maintenance command that regenerates and compares them.

**Milestone F exit check:** independent stator/rotor instance counts work, machine
observables meet tolerances, and position sweeps neither remesh templates nor
reimport solver topology.

## Milestone G — Native compressed solver path

Milestone result: the solver consumes logical instances without storing a fully
expanded mesh, with the materialized path retained as an oracle.

- [ ] **G1: Add a logical global-index/DOF view.**
  - Depends on: F4.
  - Map instance-local nodes to DOFs, honor welded seams, and keep unconnected
    instance unknowns independent.
- [ ] **G2: Build adjacency and ordering from logical instances.**
  - Depends on: G1.
  - Compare matrix dimensions, adjacency, and bandwidth/profile with materialized
    topology; handle disconnected components and explicit PBCs.
- [ ] **G3: Assemble elements through the logical view.**
  - Depends on: G2.
  - Support the Stage D/E isotropic magnetostatic scope first and retain a runtime
    materialized debug path.
- [ ] **G4: Integrate AGE coupling and per-instance sources.**
  - Depends on: G3.
  - Confirm AGE positioning updates coupling without rebuilding stored topology.
- [ ] **G5: Run full native/materialized equivalence suite.**
  - Depends on: G4 and F5.
  - Compare sparsity, residual, nodal solution modulo gauge, circuit quantities,
    fields, energy, forces, and torque for every earlier fixture.
- [ ] **G6: Verify compressed-memory scaling.**
  - Depends on: G5.
  - Demonstrate stored mesh memory scales with template size plus instance metadata
    rather than logical element count. Keep materialization available for at least
    one release cycle.

**Milestone G exit check:** native and materialized results agree and the measured
stored-mesh memory shows the intended asymptotic reduction.

## Milestone H — Optional optimisation and extension backlog

These tasks begin only after Milestone G. Each requires its own benchmark and an
equivalence test against the unoptimised native path.

- [ ] **H1:** cache rigidly rotated isotropic element contributions.
- [ ] **H2:** batch assembly by template and physics-override class.
- [ ] **H3:** design a versioned persistence format if user workflows require it.
- [ ] **H4:** transform anisotropic tensors correctly.
- [ ] **H5:** specify and implement reflection semantics.
- [ ] **H6:** add higher-level mfemm and RNFoundry authoring helpers.

## Recommended first work sequence

Start with this strictly ordered slice; it establishes the real default backend
before any instancing types depend on it:

1. A1.1 — agree the upstream Tangle library/input contract.
2. A1.2 and A1.3 — pin/build/package Tangle reproducibly.
3. A2.1 — implement and unit-test the value converter.
4. A2.2 and A2.3 — execute real Tangle and map options/diagnostics.
5. A3.1 — centralize mesh invariant checks.
6. A3.2 and A3.3 — establish ordinary, PBC, and AGE differential confidence.
7. A3.4 — make Tangle the default.
8. B1 — land `MeshingRequest` before developing topology-only seam matching.

Do not start by changing `AnalysisSession`'s default or writing `InstancedMesh`:
until A2/A3 prove the real engine and conversion path, those changes would build
new behavior on the current Triangle-delegating Tangle facade.
