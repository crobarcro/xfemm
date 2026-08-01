# xfemm Python interface

The Python package provides one stateful class, `xfemm.FemmSession`. It loads a
magnetic `.fem` model, meshes and solves it through xfemm's native C++ session,
and exposes NumPy results without requiring MATLAB or Octave.

```python
from xfemm import FemmSession

with FemmSession("model.fem") as session:
    session.setCircuit("phase-a", "current", 10)
    status = session.solve()
    potential = session.geta([[0, 0], [1, 0]])
```

This is an initial implementation. Session lifecycle, compact/full results,
mesh queries, and unsmoothed `geta`/`getb` point queries are available. The
remaining post-processing methods documented in
`docs/python-session-interface.md` will be added as the native solution snapshot
is expanded.

