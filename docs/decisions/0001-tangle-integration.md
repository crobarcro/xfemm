# Decision 0001: consume Tangle as a pinned library dependency

- Status: accepted
- Date: 2026-08-15
- Tasks: A1.1–A1.3

## Context

Instanced-mesh development uses [Tangle](https://github.com/dcm3c/tangle) as its
production meshing engine.  xfemm needs reproducible online builds, an offline
development/release path, an in-process API, and correct redistribution of
Tangle's MIT license.

At the reviewed revision, Tangle is a C++17 project consisting of
`tangle.cpp`, `float256.cpp`, and public `tangle_mesh.h`.  It does not yet provide
a CMake project or package configuration.  Defining `TANGLE_AS_LIBRARY` removes
its CLI entry point, and `tangle_mesh_fem` returns a `Mesh` in memory.  That entry
point still reads a FEMM file path, so it is not sufficient for the final xfemm
backend contract; the follow-up backend task must add/use an upstream in-memory
problem/PSLG entry point rather than introducing temporary files.

## Decision

xfemm defines a thin CMake static-library target around the exact upstream sources
at revision `a808c624ec0584569e43593662f54890b602c6af`.

- `TANGLE_PROVIDER=auto` uses `TANGLE_SOURCE_DIR` when supplied and otherwise
  fetches the pinned Git revision.
- `TANGLE_PROVIDER=local` requires `TANGLE_SOURCE_DIR`, providing a deterministic
  offline and upstream-development workflow.
- `TANGLE_PROVIDER=fetch` always fetches the pinned revision, ignoring a local
  source hint so this mode explicitly exercises dependency retrieval.
- Configuration fails for invalid providers, missing local paths, or incomplete
  source trees. It never substitutes Triangle while reporting Tangle.
- Tangle is compiled with `TANGLE_AS_LIBRARY` and linked privately by `fmesher`.
- Tangle's `LICENSE` is installed as `doc/xfemm/LICENSE-tangle.txt`.

The target is intentionally owned by xfemm until upstream supplies a suitable
CMake package.  When it does, a future change may prefer `find_package(Tangle)`
while retaining the exact-revision and offline-provider guarantees.

## Rejected alternatives

### Copy Tangle sources into xfemm

This obscures the upstream revision, encourages snapshots to drift, and makes
security/license updates harder to audit.

### Invoke the Tangle executable

This requires temporary problem and mesh files, weakens typed error handling, and
violates the in-memory `MesherBackend` contract.

### Silently use Triangle if Tangle is unavailable

This makes backend diagnostics and differential tests unreliable.  A requested
Tangle configuration must either execute Tangle or fail clearly.

## Consequences

The xfemm CMake minimum is now 3.14 so `FetchContent` can provide the pinned online
source.  Normal offline builds must set `TANGLE_SOURCE_DIR` (usually together with
`TANGLE_PROVIDER=local`).  Merely linking the library does not complete the real
backend: A2.1/A2.2 remain responsible for value conversion and replacing the
current Triangle-delegating facade.
