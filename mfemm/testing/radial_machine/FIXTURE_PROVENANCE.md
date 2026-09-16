# Radial-machine fixture provenance

This file records how the checked-in `data/` fixtures were generated. It is the
reference for milestone F (RNFoundry machine validation) of the instanced-mesh
plan.

## Generator and revision

- Generator: `generate_radial_machine_fixtures.m` (this directory).
- RNFoundry revision: `f4d428055cc1b9f4e121a0c2761c819c8dd5fe7f`
  (`crobarcro/rnfoundry`). This deliberately updates the reconnaissance pin
  `9585b26767175dc53a7451acb12c38068d45b256` recorded in the implementation
  plan; the older fixtures did not reproduce the current design and produced
  non-finite circuit flux linkage in the redraw case.
- Machine model: RNFoundry `slottedfemmprob_radial`, 12-pole / 36-slot design,
  from `radial_machine_static_case.m` with `SolveMethod = xfemm_legacy` and
  `RotationMethod = MagnetRedraw`.
- Design/options: `NWindingLayers = design.CoilLayers`,
  `SplitSlot = (CoilLayers == 2) && (yd == 1)`, `CoilCurrent = zeros(Phases,1)`,
  and the `MagFEASim` region mesh sizes passed by the generator.

## Generated files

- `data/radial_machine_redraw_01.fem` .. `_10.fem` — independently drawn and
  meshed (MagnetRedraw) at `positions = linspace(0, 1, 10)` pole pitches.
- `data/radial_machine_sliding.fem` — one AGE/sliding-mesh machine at position 0.
- `data/positions.txt` — the ten pole-pitch fractions.

## Result schema and recorded quantities

`radial_machine_fixture_case.m` writes `schemaVersion = 2` with:

- `positions`, `fluxLinkage` (gauge-invariant winding flux linkage from
  opposing coil sides), `circuitFluxLinkage` (direct FEMM circuit flux
  linkage), and `coilFluxDensity` (flux-density magnitude at the coil sample
  points).

Tolerances are justified from the existing redraw-versus-sliding comparison and
are stored beside the fixture comparison helpers. They cannot be finalised until
the sliding-mesh session path produces finite results (see below).

## Regeneration

RNFoundry is a fixture-generation dependency only; ordinary CI consumes the
checked-in fixtures and never needs RNFoundry. To regenerate:

```matlab
addpath('/path/to/xfemm/mfemm');
addpath('/path/to/xfemm/mfemm/testing/radial_machine');
generate_radial_machine_fixtures('/path/to/rnfoundry', ...
                                 '/path/to/xfemm/mfemm/testing/radial_machine/data');
```

The generator records no external revision in its output; update the RNFoundry
revision above whenever it is rerun.

## Known blockers

- Independent stator/rotor templates (F3) require a multi-template
  AGE-coupling extension to the materializer: the current instancing path
  supports one rotational template and remaps each template's own AGE, but
  does not couple two templates through the air-gap rings.

## Investigation log (redraw versus sliding)

Two post-processing bugs were found and fixed while establishing this
comparison:

1. The in-memory post-processor used the solver's post-expansion series
   circuits, so a zero-current series circuit had no member labels and
   `fpproc::GetFluxLinkage` returned NaN. Fixed in
   `cfemm/fpproc/InMemorySolution.cpp`.
2. `loadfemmsolution` kept the surrounding quotes that `textscan` returns for
   `MagDirFctn`, so a `loadfemmfile`/`writefemmfile` round-trip wrote
   `""theta""`. The legacy redraw solve then evaluated a malformed
   magnetization-direction expression and produced a field roughly 2x wrong,
   which the winding-flux-linkage cancellation amplified to ~1000x. Fixed in
   `mfemm/loadfemmsolution.m`.

After both fixes the legacy redraw, the session redraw (Triangle and Tangle),
and the sliding session agree on the magnet, air-gap, and stator flux density
to within a few percent, and the `coil-region flux density` comparison passes.

The `coil flux linkage` comparison still exceeds its tolerance
(`max |delta| ~ 9.5e-7` against a `~1.5e-7` limit). The quantity is a small
difference of large cancelling `Turns * intA / area` terms (~2e-6 against
~0.14 terms), so it amplifies the remaining `MagnetRedraw`-versus-`SlidingMesh`
modelling difference. Running both sides through the same session/Tangle
backend still leaves a `~5e-7` difference, so the residual is the modelling
difference rather than the mesher or the post-processor. A tolerance decision
(or a more robust winding-flux-linkage metric) is needed before F4/F5.
