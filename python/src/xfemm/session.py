"""Public single-class xfemm session API."""

from __future__ import annotations

from os import PathLike
from typing import Any

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
