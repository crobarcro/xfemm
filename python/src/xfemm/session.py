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
        except RuntimeError as exc:
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

    def nummeshnodes(self) -> int:
        return self._require_open().nummeshnodes()

    def numelements(self) -> int:
        return self._require_open().numelements()

    def getvertices(self, n: Any = None) -> np.ndarray:
        return self._require_open().getvertices(_indices(n))

    def getelements(self, n: Any = None) -> np.ndarray:
        return self._require_open().getelements(_indices(n))

    def getcentroids(self, n: Any = None) -> np.ndarray:
        return self._require_open().getcentroids(_indices(n))

    def getareas(self, n: Any = None) -> np.ndarray:
        return self._require_open().getareas(_indices(n))

    def geta(self, x: Any, y: Any = None) -> np.ndarray:
        points = _points(x, y)
        native = self._require_open()
        try:
            return native.geta(points)
        except RuntimeError as exc:
            raise XfemmError(str(exc)) from exc

    def getb(self, x: Any, y: Any = None) -> np.ndarray:
        points = _points(x, y)
        native = self._require_open()
        try:
            return native.getb(points)
        except RuntimeError as exc:
            raise XfemmError(str(exc)) from exc


def _indices(value: Any) -> np.ndarray:
    if value is None:
        return np.empty(0, dtype=np.int64)
    original = np.asarray(value)
    if original.dtype == np.bool_:
        raise TypeError("element indices cannot be boolean")
    result = np.asarray(value, dtype=np.int64).reshape(-1)
    if not np.all(np.asarray(value).reshape(-1) == result):
        raise ValueError("element indices must be integers")
    return np.ascontiguousarray(result)


def _points(x: Any, y: Any) -> np.ndarray:
    if y is None:
        result = np.asarray(x, dtype=np.float64)
        if result.ndim == 1 and result.size == 2:
            result = result.reshape(1, 2)
        if result.ndim < 2 or result.shape[-1] != 2:
            raise ValueError("a single coordinate argument must have final dimension 2")
        return np.ascontiguousarray(result.reshape(-1, 2))
    xa, ya = np.broadcast_arrays(np.asarray(x, dtype=np.float64),
                                 np.asarray(y, dtype=np.float64))
    return np.ascontiguousarray(np.column_stack((xa.reshape(-1), ya.reshape(-1))))
