from pathlib import Path

import numpy as np
import pytest

from xfemm import FemmSession, SessionStateError


MODEL = Path(__file__).parents[2] / "cfemm" / "fsolver" / "test" / "Temp.fem"
UNIFORM_MODEL = Path(__file__).parent / "data" / "uniform.fem"


def test_session_lifecycle_and_numpy_results():
    with FemmSession(MODEL) as session:
        count = session.mesh()
        assert count > 0
        status = session.solve()
        assert status["success"] is True
        assert status["elementCount"] == count

        result = session.result()
        assert result["A"].shape == (status["nodeCount"],)
        assert result["x"].dtype == np.float64
        session.accept()

    with pytest.raises(SessionStateError, match="closed"):
        session.mesh()


def test_result_requires_solve():
    with FemmSession(MODEL) as session:
        with pytest.raises(SessionStateError, match="solve first"):
            session.result()


def test_uniform_field_solve_is_in_memory(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    with FemmSession(UNIFORM_MODEL) as session:
        status = session.solve()
        assert status["nodeCount"] == session.result()["A"].size
        x = np.array([0.25, 0.5, 0.75])
        np.testing.assert_allclose(session.geta(x, 0.5).real, 0.25 * x, atol=1e-7)
        np.testing.assert_allclose(session.getb([[0.5, 0.5]]).real,
                                   [[0.0, -0.25]], atol=1e-7)

    # Solving must not use a hidden legacy solution round trip.
    assert list(tmp_path.iterdir()) == []
