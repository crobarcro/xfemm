# Decision 0002: extend Tangle's library API for options and boundary matches

- Status: accepted (temporary fork pin until upstreamed)
- Date: 2026-09-15
- Tasks: A2.3, B1–B5

## Context

Milestone A2.3 must map xfemm's `MeshingOptions` onto the real engine rather
than ignore them, and Milestone B needs matched boundary chains to be available
as topology without implicitly emitting periodic field constraints.

The pinned Tangle revision `a808c624` exposes only:

```cpp
int tangle_mesh_fem(const std::string& inputBase, Mesh& outMesh);
```

It has no option-accepting entry point and no ordered boundary-chain output.
The latest upstream `dcm3c/tangle` `master` (`66221c3`, v0.4.18) is the same:
`tangle_mesh.h` declares only that one function. The capability does not exist
upstream at any revision, so xfemm cannot simply move the pin.

## Decision

The required capability is added to the Tangle fork at
`https://github.com/crobarcro/tangle` on branch `xfemm-library-api` (current
commit `5ae72ff4dcfe88e496b0e8098f4bb796abd7fee4`):

- `MeshOptions` and an overload
  `tangle_mesh_fem(path, const MeshOptions&, Mesh&)` apply caller overrides
  (minimum angle, maximum element area, boundary Steiner suppression,
  unused-vertex jettison, verbosity) on top of the problem-derived defaults.
- `Mesh::boundary_matches` carries the ordered one-to-one correspondence
  between each declared pair of matched chains. The chains are extracted after
  refinement from the maintained `pbc_twin` map, so they reflect the final
  synchronously split discretisation and do not depend on `pbc_pairs`.
- `FemProblem` and an overload
  `tangle_mesh_fem(const FemProblem&, const MeshOptions&, Mesh&)` mesh a
  FEMM-like problem already held in memory, so a caller can mesh a tile it has
  parsed without a file round trip. `readFemFile` was factored into a stream
  reader (`readFemStream`) plus a thin file wrapper; the record entry point
  shares the exact file parser (arc discretization, LFS, PBC and AGE handling)
  and emits the same `boundary_matches`. The file overloads are unchanged.

xfemm temporarily pins `XFEMM_TANGLE_REPOSITORY` and the MEX
`gettanglesourcedir` helper to that fork commit. Both are documented as
temporary; when the changes are merged upstream, the repository returns to
`dcm3c/tangle` and the revision to the merged commit. No Tangle source is copied
or patched at configure time, and `TANGLE_PROVIDER=fetch` still exercises the
exact pinned revision.

On the xfemm side:

- `MeshingRequest` and `BoundaryMatch` describe topology-only matches using
  stable source-entity references; `createFieldConstraint` selects whether a
  match also becomes a periodic or antiperiodic field constraint.
- `MeshResult::boundaryMatches` exposes the ordered correspondence as
  backend-neutral seam data.
- `BoundaryMatchValidator` rejects stale, mixed line/arc, duplicated, missing,
  and non-periodic references before any backend runs.
- The Tangle backend emits periodic constraints only when requested; the
  Triangle backend reports `MeshStatus::Unsupported` for topology-only matches
  rather than emulating them with periodic physics.

## Rejected alternatives

### Wait for an upstream API before completing A2.3/B3

No upstream revision provides it, and the xfemm-side instancing work depends on
the capability now.

### Patch the fetched Tangle source at configure time

This is effectively an untracked fork, is fragile across the CMake and MEX
providers, and hides the upstream revision.

### Keep the file-based facade and ignore options

The task explicitly requires mapping supported options and diagnosing the rest
instead of silently ignoring them.

## Consequences

The dependency metadata, install license, and offline/`TANGLE_SOURCE_DIR`
workflows are unchanged; only the repository URL and revision move. Once the
fork changes are merged upstream, revert the repository/revision in
`cfemm/fmesher/CMakeLists.txt` and `mfemm/private/gettanglesourcedir.m` and
update this record.
