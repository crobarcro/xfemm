# Decision 0003: the loaded problem is the rotational tile

- Status: accepted (initial milestone D scope)
- Date: 2026-09-15
- Tasks: D1–D6

## Context

Milestone D repeats a small stator or rotor tile by rotation. The plan's D2
describes extracting "one tile PSLG" from a larger model and preserving stable
mappings back to source entities.

Tangle's library entry point still accepts a FEMM file path only. It has no
in-memory `FemmProblem`/PSLG input at the pin or at upstream `master`
(decision 0002), and xfemm deliberately avoids temporary problem files on the
backend path. Reimplementing Tangle's FEMM reader (arc discretization, LFS,
PBC/AGE tagging) in xfemm, or adding group-filtered tile selection upstream,
is a larger change than the rest of milestone D.

## Decision

For the initial implementation, `TemplateRequest` treats the loaded problem as
the tile:

- the problem's own PSLG, markers, block-region attributes, holes, mesh
  controls, and units are handed to Tangle unchanged;
- its declared periodic boundaries are meshed with topology-only boundary
  matching (milestone B) and become the template seams;
- `TemplateRequest` carries the rotation centre, instance count, total angle,
  and optional seam boundary-property selection.

Consequently D2's "isolated tile PSLG" is the problem PSLG and per-source-entity
mappings are not produced. Tile selection from a larger model is deferred.

## Rejected alternatives

### Write a temporary tile `.fem` for Tangle

Violates the no-file-round-trip backend contract and reintroduces the failure
mode milestone B removed.

### Reimplement the FEMM-to-PSLG conversion in xfemm

Duplicates Tangle's reader and its arc/LFS/PBC/AGE handling, which would drift.

### Select tile geometry by group id and mesh the full problem

Does not reduce meshing work and would not prove the "triangulate the template
once" property.

## Consequences

A user meshes one tile model (for example one pole pitch) and xfemm repeats it.
`TemplateRequest` already has the shape needed for selection, so adding
`groupIds`/entity references later is additive. The D6 control test demonstrates
that a tile meshed once and repeated eight times agrees with a conventionally
meshed full ring, which is the milestone's exit criterion.
