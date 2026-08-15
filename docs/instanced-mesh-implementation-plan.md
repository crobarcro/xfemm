# Instanced mesh implementation plan

## Purpose

This plan turns the proposed rotary-machine meshing feature into incremental,
reviewable work.  The central idea is to mesh a small stator or rotor tile once,
repeat that topology under rigid transforms, and give each occurrence independent
field degrees of freedom.  Mesh repetition is therefore **not** the same thing as
periodic or antiperiodic field symmetry.

The design discussion that motivated the work is preserved in the
[shared ChatGPT conversation](https://chatgpt.com/share/6a7f0e14-117c-83eb-b8a4-3318d91d9154).
The first implementation deliberately materialises an ordinary `SolverMesh` so
that the magnetic solver remains unchanged.  The default meshing engine for this
work is [Tangle](https://github.com/dcm3c/tangle), a C++17 constrained Delaunay
mesher with native FEMM, synchronized periodic-boundary, arc, and AGE support.
Triangle remains a compatibility and differential-testing backend while Tangle
becomes the production path.  Compressed storage and solver optimisations come
only after equivalence has been demonstrated.

Repository reconnaissance for this plan used
[crobarcro/rnfoundry at `9585b267`](https://github.com/crobarcro/rnfoundry/commit/9585b26767175dc53a7451acb12c38068d45b256).
In particular, xfemm's fixture generator calls RNFoundry's
`slottedfemmprob_radial`, while the checked-in fixtures allow the regression suite
to run without cloning RNFoundry.  Recheck and deliberately update that pin when
Stage 4 begins; do not silently generate baselines from a moving default branch.

## Outcomes and non-goals

### Required outcomes

1. A template mesh can be repeated by rotation (and, in the data model, by any
   orientation-preserving rigid transform).
2. Adjacent instances have one conforming seam: coincident seam nodes are welded,
   while non-seam nodes and their unknowns remain independent.
3. Element regions, boundary markers, ordinary periodic constraints, and AGE
   ring references survive materialisation with zero-based indices.
4. Different stator and rotor templates can meet through the existing AGE
   coupling without requiring equal instance counts or conforming air-gap meshes.
5. Per-instance circuits, current signs, and magnetisation rotations can differ
   even when instances share mesh topology.
6. A conventionally meshed machine and a materialised instanced machine agree on
   topology invariants and selected magnetic results within stated tolerances.

### Deferred work

- Do not remove Triangle while introducing Tangle; retain it as a fallback and
  independent comparison backend through native-instancing validation.
- Do not require the legacy `FSolver` to consume compressed templates initially.
- Do not identify fields across instance seams as a consequence of instancing.
  Conventional periodic constraints remain an explicit, separate choice.
- Do not add an instancing extension to the `.fem`, `.node`, `.ele`, `.edge`, or
  `.pbc` formats until the in-memory API and semantics are stable.
- Do not promise reflection transforms in the first release.  Reflections reverse
  element orientation and magnetisation handedness and need separate semantics.

## Existing seams to preserve and exploit

The current pipeline already has an appropriate compatibility boundary:

```text
FemmProblem
  -> MesherBackend::mesh(...)
  -> TangleMesherBackend (default)
  -> Tangle in-memory meshing API
  -> SolverMesh
  -> AnalysisSession cache
  -> FSolverAnalysisBackend::synchronize()
  -> FSolver::LoadMesh()
```

`SolverMesh` is a value-only transfer object with zero-based connectivity.  It
holds nodes, elements, edges, ordinary periodic node pairs, and AGE topology.
`FSolverAnalysisBackend` imports that object only when the session topology
identity changes.  These properties make materialisation immediately upstream of
`SolverMesh` the least disruptive first insertion point.

Tangle's synchronized periodic-boundary splitting is the primary starting point
for template seams: paired chains are refined together and retain ordered
node-to-node correspondence.  The xfemm integration must generalise that facility
so a boundary match can request identical discretisation without automatically
creating a periodic field constraint.  The legacy periodic Triangle workflow,
which discovers boundary subdivisions through a trial mesh and then reconciles
paired segment/arc spacing, remains the behavioural reference for existing
models, not the implementation base for new instancing work.

The AGE data must remain a coupling boundary rather than be treated as evidence
of instancing.  Ordinary periodic pairs occupy the primary `.pbc` section and AGE
records are only an optional extension; both must be remapped independently when
a mesh is materialised.

### Tangle dependency and backend policy

The canonical Tangle project is housed at
[`dcm3c/tangle`](https://github.com/dcm3c/tangle).  Integration work must record an
exact upstream commit and its MIT license in xfemm's dependency metadata.  Prefer
a normal CMake library dependency fetched or supplied at configure time over
copying untracked snapshots of `tangle.cpp` into xfemm.  Release and offline builds
must have a documented way to use a pinned vendored source archive, and installed
license material must include Tangle's license.

This plan was checked against Tangle revision
[`a808c624`](https://github.com/dcm3c/tangle/commit/a808c624ec0584569e43593662f54890b602c6af),
whose public `tangle_mesh_fem` entry point returns an in-memory `Mesh` but currently
accepts a FEMM input path.  Stage 0 should pin the then-current reviewed revision
and coordinate any required library/CMake and in-memory-input API changes in the
Tangle repository.

The present `TangleMesherBackend` is only a compatibility facade: it delegates to
`TriangleMesherBackend` and relabels diagnostics.  Stage 0 must replace that
delegation with Tangle's real C++ API and convert Tangle's `Mesh`, PBC pairs, and
AGE definitions directly into xfemm's value-only `SolverMesh`.  No temporary mesh
files are allowed on the backend API path.  If the upstream API cannot yet accept
an in-memory `FemmProblem`/PSLG, make the required reusable input API an upstream
Tangle change rather than adding file round-tripping to xfemm.

Once differential tests pass, `AnalysisSession` and new instancing entry points
select `TangleMesherBackend` by default.  `TriangleMesherBackend` stays explicitly
selectable for compatibility, baseline generation, and fault isolation.  Backend
diagnostics must report the engine actually executed; a Triangle result must never
be labelled as Tangle.

## Proposed contracts

Names below are provisional, but each concept and ownership boundary is part of
the plan.

### Mesh representation (`cfemm/libfemm/mesh`)

```cpp
struct RigidTransform2D {
    double rotationDegrees;
    double translationXMetres;
    double translationYMetres;
};

struct TemplateSeam {
    std::string name;
    std::vector<MeshIndex> orderedNodes;
};

struct MeshTemplate {
    SolverMesh localMesh;
    std::vector<TemplateSeam> seams;
};

struct MeshInstance {
    std::size_t templateIndex;
    RigidTransform2D transform;
    // A seam on this instance is welded to a seam on another instance.
    std::vector<SeamConnection> seamConnections;
};

struct InstancedMesh {
    std::vector<MeshTemplate> templates;
    std::vector<MeshInstance> instances;
};
```

`InstancedMesh::materialize()` returns a `SolverMesh` plus provenance maps:

- `(template, local node, instance) -> global node`;
- `(template, local element, instance) -> global element`;
- global node/element -> template and instance, for diagnostics and post-processing.

Seam connections state orientation explicitly (`Forward` or `Reverse`) and are
validated before materialisation.  Welding uses those ordered index lists, not a
global coordinate-tolerance search.  Coordinates are checked within a documented
metre-scale tolerance only to reject a bad connection.  This avoids accidentally
welding close but physically distinct geometry.

During materialisation:

- connected seam nodes map to the same global node;
- interior and unconnected boundary nodes receive distinct global indices;
- degenerate and orientation-reversed elements are rejected;
- duplicate seam edges are canonicalised, while physical boundary markers are
  rejected if the two sides disagree;
- all periodic and AGE indices are remapped through the node provenance table;
- markers and mesher region attributes remain unchanged at this layer.

### Meshing request (`cfemm/libfemm/mesh/Meshing.h`)

Replace the ambiguous `bool periodic` parameter, compatibly and in a separate
change, with a request object:

```cpp
struct MeshingRequest {
    MeshingOptions options;
    bool createPeriodicFieldConstraints;
    std::vector<BoundaryMatch> boundaryMatches;
    std::vector<TemplateRequest> templates;
};
```

`BoundaryMatch` means “produce the same ordered boundary discretisation.”  A flag
on the match separately says whether to emit a periodic/antiperiodic field
constraint.  `TemplateRequest` identifies source geometry by stable entity/group
identifiers and declares its seams and transforms.  During migration, the old
`mesh(problem, bool, options)` overload constructs a `MeshingRequest`, so CLI and
MATLAB callers do not change.

Do not encode a template match by inventing a magnetic `BdryFormat`: electrostatic
and heat-flow formats differ, and matching mesh topology is not a physical
boundary condition.

### Physics metadata

Keep physics overrides out of `MeshTemplate`.  Add model-side metadata such as:

```cpp
struct InstanceRegionOverride {
    std::size_t sourceBlockLabel;
    std::optional<std::size_t> circuit;
    std::optional<double> magnetisationRotationDegrees;
    std::optional<double> currentScale;
};
```

This contract uses `std::optional`, consistent with xfemm's required C++17
language level and the optional-value conventions already used by
`AnalysisSession`.

Material identity, circuit assignment, turns/current sign, and magnetisation
direction are resolved while building `PreparedAnalysis`, using element
provenance.  A rigid transform rotates directional quantities; it must not mutate
the shared template.  Isotropic material support is the initial acceptance case.
Anisotropic materials are rejected until their tensor transform is implemented.

## Staged delivery

Every stage ends with a usable, tested state.  Later stages must not be started by
silently weakening an earlier stage's equivalence checks.
The issue-sized implementation order, dependencies, and completion checklist are
maintained in the companion
[instanced mesh task breakdown](instanced-mesh-task-breakdown.md).

### Stage 0 — Integrate Tangle and characterise boundary matching

**Goal:** make the real Tangle engine the tested default and expose its synchronized
boundary refinement independently of field periodicity before adding a new mesh
representation.

Work:

1. Pin and integrate `dcm3c/tangle` as a CMake library dependency, install its
   license, and replace the compatibility delegation in `TangleMesherBackend` with
   checked conversion between Tangle and `SolverMesh` value types.
2. Add focused tests for ordered segment and arc matches in both forward and
   reverse orientation, including unequal geometry rejection.
3. Add characterization tests proving that an ordinary periodic model produces
   node pairs but no AGE, and that the machine fixture produces AGE topology.
4. Extend Tangle's paired-chain refinement API with a topology-only boundary-match
   mode that returns ordered correspondence without emitting PBC field semantics.
   Keep the legacy Triangle workflow unchanged as the differential baseline.
5. Introduce `MeshingRequest` and retain the legacy overload as an adapter.
6. Compare Tangle and Triangle on the checked-in ordinary, periodic, and AGE
   fixtures.  Compare invariants and solver observables rather than requiring
   identical triangulations.
7. Switch `AnalysisSession`'s default construction to `TangleMesherBackend` only
   after those comparisons pass; retain explicit Triangle injection.

Exit criteria:

- Existing mesher, solver, femmcli, and MATLAB tests pass with Tangle as the
  default and also pass their designated Triangle compatibility configurations.
- Instrumentation proves that the Tangle backend executes Tangle rather than the
  Triangle compatibility kernel and creates no temporary files.
- Tangle can match two boundaries without emitting a field constraint.
- Invalid count, length, arc geometry, and orientation inputs produce diagnostics
  rather than partial meshes.

### Stage 1 — Value model and deterministic materialisation

**Goal:** prove mesh repetition independently of `FemmProblem` and any meshing
engine.

Work:

1. Add `MeshTemplate`, `TemplateSeam`, `MeshInstance`, `InstancedMesh`, provenance,
   validation, and `materialize()` under `cfemm/libfemm/mesh`.
2. Add small hand-authored meshes: two translated triangles, an annular wedge
   repeated by rotation, reversed seam ordering, a closed ring, and intentionally
   invalid seam declarations.
3. Test remapping of boundary markers, region attributes, ordinary periodic pairs,
   and synthetic AGE ring/quadrature/node-index data.
4. Add stable topology hashing based on templates, transforms, and seam maps (not
   pointer identity or floating-point object addresses).

Exit criteria:

- Repeating an annular wedge through 360 degrees has the predicted node, edge,
  and element counts, no seam cracks, no duplicate elements, and positive element
  orientation.
- Materialising the same input twice is deterministic.
- Address/undefined-behaviour sanitizers pass the new unit tests.

### Stage 2 — Tangle-generated templates

**Goal:** mesh one geometry tile and produce a solver-compatible full mesh.

Work:

1. Extend `TangleMesherBackend` to honour one `TemplateRequest` using Tangle's
   topology-only boundary matching from Stage 0.
2. Mesh only the selected tile geometry, capture each declared seam's ordered
   returned nodes, construct `InstancedMesh`, and materialise it.
3. Initially support one planar rotational template about a declared centre,
   same-template neighbour seams, isotropic regions, and no AGE inside a template.
4. Return the materialised `SolverMesh` in the existing `MeshResult`; optionally
   expose instancing/provenance in an additive field so existing consumers compile.
5. Add diagnostics for incomplete rotations, overlapping instances, unmatched
   seams, conflicting markers, and unsupported transforms/materials.

Exit criteria:

- A repeated full ring solves through the unmodified `FSolver`.
- Its mesh invariants and sampled field results agree with a conventionally
  generated full-ring control model at documented tolerances.
- Tangle triangulates the template once (verified with an injectable/counting test
  backend), even though the returned compatibility mesh is expanded.

### Stage 3 — Session integration and per-instance physics

**Goal:** make instancing a first-class in-memory analysis option without changing
the legacy solver import.

Work:

1. Let `AnalysisSession` retain the canonical `InstancedMesh`, its provenance, and
   a lazily cached materialised `SolverMesh`.
2. Separate template-topology identity, instance-layout identity, and materialised
   solver-topology identity so circuit/current changes do not remesh.
3. Resolve `InstanceRegionOverride` entries while rebuilding prepared materials
   and circuits.  Rotate magnetisation by the instance transform.
4. Expose a C++ API first.  Add MATLAB/MEX bindings only after that API has passed
   the radial-machine test; preserve existing calls.
5. Extend post-processing selection/provenance so a caller can identify an
   instance even though the legacy post-processor sees global elements.

Exit criteria:

- Two instances of one coil-slot template can use opposite turns/current signs
  and produce independent unknowns.
- Alternating permanent-magnet directions match an explicitly drawn control.
- Changing current or AGE angle does not remesh or rematerialise topology;
  changing an instance transform invalidates only the necessary caches.

### Stage 4 — Independent stator/rotor templates and RNFoundry validation

**Goal:** demonstrate the intended electrical-machine workflow.

Use the pinned RNFoundry revision recorded in the fixture generator/README when
this stage is implemented.  RNFoundry is a fixture-generation dependency only;
checked-in `.fem` inputs and compact reference results keep xfemm CI independent.
Start from its `slottedfemmprob_radial` machine model and xfemm's existing
`mfemm/testing/radial_machine` fixtures rather than introducing a second machine.

Work:

1. Generate three comparable cases from the existing 12-pole, 36-slot design:
   conventional redraw, current AGE/sliding mesh, and instanced mesh.
2. Build a stator slot template and rotor pole (or pole-pair, if required by
   magnet physics) template with independent counts.  Join each domain's internal
   seams, then couple stator and rotor only through the existing AGE rings.
3. Exercise multiple rotor positions without remeshing.  Record template and
   materialised node/element counts, meshing calls, topology imports, solve time,
   and peak resident memory; performance numbers are informative, not pass/fail,
   at this stage.
4. Compare gauge-invariant quantities already used by the fixture tests: winding
   flux linkage from opposing coil sides, coil flux-density magnitude, and torque.
   Add air-gap flux samples/harmonics if they are stable across mesh variants.
5. Check in the generator, its RNFoundry commit, all relevant options, fixture
   provenance, and tolerances.  Do not require RNFoundry in ordinary CI.

Exit criteria:

- The instanced and conventional cases agree within tolerances justified from the
  existing redraw-versus-sliding comparison.
- Stator and rotor instance counts are demonstrably independent.
- Rotor-position sweeps reuse both template meshes and update AGE coupling without
  a topology import.
- CI has a small one-position smoke test; the full sweep can run as an extended or
  release test.

### Stage 5 — Native compressed solver consumption

**Goal:** obtain the memory benefit that materialisation intentionally postpones.

Work:

1. Add a solver-side global-index view that maps `(instance, local node)` to a DOF
   and honours welded seams without duplicating template coordinates/connectivity.
2. Teach Cuthill–McKee/adjacency construction and assembly to iterate logical
   instance elements.  Keep an adapter that materialises for debugging and
   equivalence comparison.
3. Retain separate unknowns for non-welded instance nodes.  Periodic or
   antiperiodic algebraic identification is applied only from explicit
   `PeriodicConstraint` data.
4. Compare native and materialised matrix dimensions, sparsity pattern, residual,
   nodal solution (modulo gauge), forces, and torque.

Exit criteria:

- Native and materialised paths agree numerically on every Stage 1–4 case.
- Stored mesh memory scales primarily with template size plus instance metadata,
  not with the number of logical elements.
- The materialised compatibility path remains available behind a test/debug
  option for at least one release cycle.

### Stage 6 — Repeated-element optimisation and persistence

**Goal:** optimise only after the representation and solver path are proven.

Candidates:

- cache local element geometry/stiffness contributions for rigidly rotated,
  isotropic template elements;
- batch assembly by template and physics-override class;
- add a versioned in-memory/file representation if real users need persistence;
- extend to anisotropic tensor rotation, translations, and carefully specified
  reflections;
- add authoring helpers to mfemm and RNFoundry.

Each optimisation needs a benchmark showing a meaningful improvement and an
equivalence test against the unoptimised native path.

## Test matrix

| Level | Mandatory cases | Principal assertions |
| --- | --- | --- |
| Unit | transform, seam ordering, closed ring, invalid seam, index overflow | deterministic maps, positive elements, precise diagnostics |
| Mesher | Tangle/Triangle differential fixtures, line/arc match, periodic without AGE, AGE, one Tangle template | ordered conformity; correct engine diagnostics; field constraints emitted only when requested |
| Solver | independent instance fields, alternating coils/magnets, periodic control | residual and field/energy/force agreement |
| Session | current change, instance-layout change, AGE angle change | correct cache invalidation and import/meshing counters |
| MATLAB | current fixture API plus an instanced smoke case | no interface regression and correct result shapes |
| Machine | redraw, sliding, materialised instancing, native instancing | flux linkage, flux density, torque, topology reuse |

Numerical tolerances must be stored beside each fixture and justified.  Topological
properties (counts, connectivity validity, seam ownership, mapping determinism,
number of mesh calls) are exact assertions and must not use numerical tolerances.

## Pull-request sequence

Keep reviews bounded by using this order:

1. real Tangle dependency/backend integration and differential characterization;
2. `MeshingRequest` compatibility migration;
3. instanced value model and materialiser;
4. Tangle template integration;
5. session/provenance and physics overrides;
6. RNFoundry fixture generation and machine equivalence;
7. native solver view/assembly;
8. optional optimisations and persistence.

Each PR should state which exit criteria it satisfies, include before/after mesh
pipeline diagrams when contracts change, and avoid mixing numerical optimisation
with representation changes.

## Risks and explicit mitigations

- **Near-coincident nodes are welded incorrectly:** weld only declared ordered
  seams; use coordinates solely as validation.
- **Seam ownership loses a physical boundary:** define marker merge rules and
  reject conflicts rather than choosing one side.
- **Instancing accidentally imposes symmetry:** independent DOFs are the default;
  only explicit periodic constraints identify or sign-link fields.
- **Transform changes element orientation:** initially accept only
  orientation-preserving transforms and verify signed area.
- **Physics metadata is confused with mesher region attributes:** retain provenance and
  resolve per-instance physics in prepared analysis, not in `SolverMesh`.
- **AGE indices become stale after welding:** remap every ring, quadrature, and
  node-index reference through the same checked node map.
- **Legacy APIs are destabilised:** adapt old meshing calls to `MeshingRequest` and
  keep materialising at the solver boundary until Stage 5 is proven.
- **RNFoundry makes CI fragile:** pin only the external fixture generator; commit
  generated inputs/results and run normal tests without RNFoundry.
- **A full-machine benchmark is too slow:** use hand-authored unit cases and a
  one-position machine smoke test on every build; reserve full sweeps for extended
  testing.

## Definition of feature completion

The feature is complete—not merely prototyped—when Stage 5 passes: users can mesh
independent stator and rotor tiles once, assign per-instance electromagnetic
metadata, couple them through AGE, solve without expanding stored topology, and
obtain results equivalent to the conventional full mesh.  Stage 6 items are
optimisations and extensions, not prerequisites for that definition.
