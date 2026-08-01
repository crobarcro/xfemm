# xfemm Python interface

The Python package provides one stateful class, `xfemm.FemmSession`. It loads a
magnetic `.fem` model, meshes and solves it through xfemm's native C++ session,
and exposes NumPy results without requiring MATLAB or Octave.

```python
from xfemm import FemmSession

with FemmSession("model.fem") as session:
    session.setCircuit("phase-a", "current", 10)
    status = session.solve()
    nodal_solution = session.result()["A"]
```

The session constructs `FPProc` directly from the solved native state, without
writing or parsing an `.ans` file. Point values, smoothing, contours and line
integrals, block selection/integrals, problem and circuit information, and mesh
counts are thin wrappers around `FPProc`. Additional mesh-table, group-table,
and air-gap convenience bindings remain to be completed.
