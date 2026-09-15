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

- The checked-in fixtures were regenerated with the revision above and the
  `redraw` comparison now produces finite results. The `sliding` case still
  produces non-finite circuit flux linkage through `xfemm.femmsession`, with
  both the Triangle and Tangle meshers and with the AGE angle left at its
  initial value. This is independent of the mesher and of the fixture revision
  and must be resolved before F2–F5 can be completed.
- Independent stator/rotor templates (F3) additionally require a
  multi-template AGE-coupling extension to the materializer: the current
  instancing path supports one rotational template and remaps each template's
  own AGE, but does not couple two templates through the air-gap rings.
