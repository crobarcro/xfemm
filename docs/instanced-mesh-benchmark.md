# Instanced solver performance benchmark

`cfemm/fsolver/test/instanced_solver_benchmark.cpp` (target
`instanced_solver_benchmark`) measures the cost of the three magnetostatic mesh
paths on the checked-in RNFoundry radial-machine fixtures:

| method | mesh | solve |
| --- | --- | --- |
| `conventional` | mesh the full redrawn `.fem` with Tangle | materialized `FSolver` |
| `materialized` | mesh each tile once, expand `InstancedMesh` into a `SolverMesh` | materialized `FSolver` |
| `native` | mesh each tile once, build a compact `LogicalMeshView` | `Static2DNative` |

Each row reports node/element counts, stored and expanded mesh bytes, mesh/view
and solve wall time, matrix bandwidth, a field checksum (`sum A^2`), and the
process peak resident set size. Run one method per process so the RSS is not
contaminated by an earlier method.

## Running

```bash
cmake --build <build-dir> --target instanced_solver_benchmark
./test/rmbench/benchmark_instanced.sh --repeats 3 --pin 2          # 60-degree sector
./test/rmbench/benchmark_instanced.sh --full --with-conventional   # whole machine
```

The shell runner defaults to `--instance-divisor 6`, which divides every tile's
repeat count and switches the closure to periodic, turning the 360-degree
machine into a 60-degree sector. The sector keeps the same materials, magnets,
circuits, AGE coupling, and nonlinear steel as the full fixture while running in
roughly one sixth of the time. `--full` uses the complete checked-in
`radial_machine_tiled.json`.

Direct invocation exposes the individual methods and a rotor-position sweep:

```bash
cfemm/bin/instanced_solver_benchmark --method native \
    --tiled mfemm/testing/radial_machine/data/radial_machine_tiled.json \
    --instance-divisor 6 --repeats 2 --position 15 0
```

## Sample results

Linux x86-64, pinned to one core, `--repeats 2` (the median is a warm solve),
60-degree sector of `radial_machine_tiled.json` (40 413 nodes, 80 038 elements,
nonlinear 1117 steel and NdFeB magnets):

| method | nodes | elements | store MB | expand MB | mesh s | solve s | band | sum A^2 | rss MB |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| materialized | 40413 | 80038 | 6.14 | 6.14 | 1.36 | 76.49 | 1132 | 1.243 | 68.4 |
| native | 40413 | 80038 | 1.57 | 6.14 | 1.24 | 74.17 | 1132 | 1.243 | 47.7 |

The identical `sum A^2` confirms the compressed native path reproduces the
materialized field. The native path stores about a quarter of the mesh bytes and
uses about 70% of the peak RSS, because it never builds the expanded `SolverMesh`
(only the solver's per-element working arrays are sized by the logical element
count). Solve time is comparable for this excited nonlinear problem; on a
zero-excitation fixture the assembly overhead dominates and the materialized
path is faster.

## Notes

- `mesh_s` is the Tangle meshing time for the tiles (instanced) or the whole
  model (conventional); `view_s` is currently not reported separately.
- `store_MB` is `LogicalMeshView::storedByteCount()` for the native path and the
  `SolverMesh` byte size otherwise; `expand_MB` is the equivalent expanded size.
- `--position INNER OUTER` sets the AGE rotor position before solving, which is
  useful for a sweep without remeshing.
- The full 360-degree machine is intended for release/extended runs; the sector
  is the CI-friendly default.
