"""Public single-class xfemm session API."""

from __future__ import annotations

from os import PathLike
from typing import Any

import numpy as np

from ._xfemm import NativeSession


class XfemmError(RuntimeError):
    """Base class for Python-facing xfemm errors."""


class SolverError(XfemmError):
    """Raised when a native meshing or solving operation fails."""


class SessionStateError(XfemmError):
    """Raised when an operation is invalid in the current session state."""


class FemmSession:
    """Stateful magnetic xfemm analysis and post-processing session.

    Method names intentionally match ``xfemm.femmsession`` in MATLAB. Point
    results are Python-oriented: samples occupy rows.
    """

    def __init__(self, filename: str | PathLike[str]):
        try:
            self._native: NativeSession | None = NativeSession(str(filename))
        except (RuntimeError, ValueError) as exc:
            raise XfemmError(str(exc)) from exc

    def __enter__(self) -> "FemmSession":
        self._require_open()
        return self

    def __exit__(self, *_: Any) -> None:
        self.close()

    def close(self) -> None:
        self._native = None

    def _require_open(self) -> NativeSession:
        if self._native is None:
            raise SessionStateError("session is closed")
        return self._native

    def setBackend(self, name: str) -> None:
        self._require_open().set_backend(name.lower())

    def mesh(self) -> int:
        native = self._require_open()
        try:
            return native.mesh()
        except RuntimeError as exc:
            raise SolverError(str(exc)) from exc

    def setCircuit(self, name: str, constraint: str, value: complex = 0) -> None:
        self._require_open().set_circuit(name, constraint.lower(), complex(value))

    def setAGEPosition(self, name: str, innerAngle: float, outerAngle: float) -> None:
        self._require_open().set_age_position(name, innerAngle, outerAngle)

    def setFrequency(self, frequency: float) -> None:
        self._require_open().set_frequency(frequency)

    def setTime(self, time: float) -> None:
        self._require_open().set_time(time)

    def solve(self) -> dict[str, Any]:
        native = self._require_open()
        try:
            return native.solve()
        except (RuntimeError, ValueError) as exc:
            raise SolverError(str(exc)) from exc

    def result(self) -> dict[str, Any]:
        native = self._require_open()
        try:
            return native.result()
        except RuntimeError as exc:
            raise SessionStateError(str(exc)) from exc

    def accept(self) -> None:
        native = self._require_open()
        try:
            native.accept()
        except RuntimeError as exc:
            raise SessionStateError(str(exc)) from exc

    def reject(self) -> None:
        self._require_open().reject()

    def getpointvalues(self, x: Any, y: Any = None) -> np.ndarray:
        return self._require_open().getpointvalues(_points(x, y))

    def geta(self, x: Any, y: Any = None) -> np.ndarray:
        return self.getpointvalues(x, y)[:, 0]

    def getb(self, x: Any, y: Any = None) -> np.ndarray:
        return self.getpointvalues(x, y)[:, 1:3]

    def geth(self, x: Any, y: Any = None) -> np.ndarray:
        return self.getpointvalues(x, y)[:, 5:7]

    def smoothon(self) -> None:
        self._require_open().smooth(True)

    def smoothoff(self) -> None:
        self._require_open().smooth(False)

    def smooth(self, flag: str) -> None:
        if flag not in ("on", "off"):
            raise ValueError("flag must be 'on' or 'off'")
        self._require_open().smooth(flag == "on")

    def clearcontour(self) -> None:
        self._require_open().clearcontour()

    def addcontour(self, x: Any, y: Any = None) -> None:
        self._require_open().addcontour(_points(x, y))

    def newcontour(self, x: Any, y: Any = None) -> None:
        self.clearcontour()
        self.addcontour(x, y)

    def lineintegral(self, type: int) -> np.ndarray:
        return self._require_open().lineintegral(type)

    def selectblock(self, x: float, y: float, clearselected: bool = False) -> None:
        if clearselected:
            self.clearblock()
        self._require_open().selectblock(x, y)

    def groupselectblock(self, groupno: Any, clearselected: bool = False) -> None:
        if clearselected:
            self.clearblock()
        groups = np.asarray(groupno, dtype=np.int64).reshape(-1)
        if groups.size == 0:
            self._require_open().selectallblocks()
            return
        for group in groups:
            self._require_open().groupselectblock(int(group))

    def selectallblocks(self) -> None:
        self._require_open().selectallblocks()

    def clearblock(self) -> None:
        self._require_open().clearblock()

    def blockintegral(self, type: int, x: Any = None, y: Any = None) -> complex:
        if x is not None:
            points = _points(x, y)
            self.clearblock()
            for px, py in points:
                self.selectblock(px, py)
        return self._require_open().blockintegral(type)

    def totalfieldenergy(self) -> complex:
        self.selectallblocks()
        return self.blockintegral(2)

    def totalfieldcoenergy(self) -> complex:
        self.selectallblocks()
        return self.blockintegral(17)

    def getprobleminfo(self) -> np.ndarray:
        return self._require_open().getprobleminfo()

    def getcircuitprops(self, circuitname: str) -> np.ndarray:
        return self._require_open().getcircuitprops(circuitname)

    def circuitRL(self, circuitname: str) -> tuple[complex, complex]:
        values = self.getcircuitprops(circuitname)
        return values[1] / values[0], values[2] / values[0]

    def nummeshnodes(self) -> int:
        return self._require_open().nummeshnodes()

    def numelements(self) -> int:
        return self._require_open().numelements()


def _points(x: Any, y: Any) -> np.ndarray:
    if y is None:
        points = np.asarray(x, dtype=np.float64)
        if points.ndim == 1 and points.size == 2:
            points = points.reshape(1, 2)
        if points.ndim < 2 or points.shape[-1] != 2:
            raise ValueError("a single coordinate argument must have final dimension 2")
        return np.ascontiguousarray(points.reshape(-1, 2))
    xa, ya = np.broadcast_arrays(np.asarray(x, dtype=np.float64), np.asarray(y, dtype=np.float64))
    return np.ascontiguousarray(np.column_stack((xa.reshape(-1), ya.reshape(-1))))
