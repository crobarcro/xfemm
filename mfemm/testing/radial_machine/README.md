# Radial-machine static integration tests

These tests use the machine from
`example_radial_flux_permanent_magnet_machine_sim_using_femmsession.m`, but
stop after the magnetostatic position sweep. They do not perform SLM fitting,
plotting, inductance simulations, or an ODE simulation.

The recorded comparison quantities are:

- gauge-invariant coil flux linkage formed from opposing coil sides;
- direct FEMM circuit flux linkage (for identical-geometry release tests);
  and
- flux density at the coil-region sample points.

The checked-in `.fem` fixtures make both tests self-contained. RNFoundry is
needed only to run `generate_radial_machine_fixtures` after the source machine
design changes.

## Test 1: sliding mesh versus redraw

Build/setup the current mfemm MEX interfaces, then run:

```matlab
addpath('/path/to/xfemm/mfemm');
mfemm_setup('RunTests', false);
Test_radial_machine_static_rotation_methods();
```

The reference side redraws and remeshes the rotor for every position with
the current xfemm. The candidate side meshes once and changes the AGE angle
inside one `xfemm.femmsession`.

For a quicker development run using three fixtures across the pole pitch:

```matlab
Test_radial_machine_static_rotation_methods([1 5 10]);
```

The sliding fixture is generated directly by RNFoundry. Its drawing code must
use one AGE property for the complete annular interface and express
`MaxSegDegrees` in degrees; the test intentionally contains no fixture-side
normalization of those requirements.

## Test 2: current xfemm versus a release

Run:

```bash
mfemm/testing/radial_machine/run_radial_machine_release_comparison.sh
```

The default reference is tag `v4.0`, currently the latest published xfemm
release. If the tag exists locally, the runner extracts it without network
access; otherwise it shallow-clones that tag. Each version runs in a separate
Octave process to avoid loading two MEX implementations with identical names.
Both sides use `MagnetRedraw`, which isolates solver-version differences from
the sliding-mesh/redraw difference covered by Test 1.

Useful environment variables are:

- `REFERENCE_XFEMM_ROOT`: use an already unpacked reference tree;
- `REFERENCE_XFEMM_TAG`: select another tag (default `v4.0`);
- `N_POSITIONS`: number of evenly distributed fixture positions (default
  `10`); and
- `KEEP_TEST_OUTPUTS=1`: retain the result MAT files and temporary release
  tree for investigation.

For GitHub Actions, install `octave`, `octave-dev`, `cmake`, and a C++
compiler, and invoke the two commands above. The release job needs
`fetch-depth: 0` if it should use the local tag rather than downloading it.
