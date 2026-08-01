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

This is an initial implementation. Session lifecycle and compact/full native
solver results are available. Post-processing methods will be exposed only
after `AnalysisSession` can hand an in-memory solution snapshot to the existing
`FPProc` C++ implementation; they will not be independently reimplemented in
the Python binding. See `docs/python-session-interface.md` for that work plan.
