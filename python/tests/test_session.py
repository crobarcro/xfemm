from pathlib import Path

import numpy as np
import pytest

from xfemm import FemmSession, SessionStateError, XfemmError


MODEL = Path(__file__).parents[2] / "cfemm" / "fsolver" / "test" / "Temp.fem"
UNIFORM_MODEL = Path(__file__).parent / "data" / "uniform.fem"
SERIES_MODEL = Path(__file__).parent / "data" / "series_circuit.fem"


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


def test_invalid_material_curve_is_reported_at_load(tmp_path):
    model = tmp_path / "invalid-curve.fem"
    text = UNIFORM_MODEL.read_text()
    text = text.replace("<BHPoints> = 0", "<BHPoints> = 1\n0 0", 1)
    model.write_text(text)

    with pytest.raises(XfemmError, match="zero points or at least two"):
        FemmSession(model)


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


def test_block_selection_change_rebuilds_stress_tensor_mask():
    with FemmSession(MODEL) as session:
        session.solve()
        session.groupselectblock(3, clearselected=True)
        selected = np.array([session.blockintegral(i) for i in (18, 19, 22)])
        session.clearblock()
        cleared = np.array([session.blockintegral(i) for i in (18, 19, 22)])

    # Temp.fem group 3 is the force-integral selection used by fpproc/main.cpp.
    # Its net weighted-stress result is near equilibrium but still nonzero, which
    # is enough to detect stale-mask reuse after the selection is cleared.
    selected_scale = np.max(np.abs(selected))
    assert selected_scale > 1e-14
    assert np.max(np.abs(cleared)) < selected_scale * 1e-3


def test_series_circuit_properties_use_user_circuit(tmp_path):
    def solve(model):
        with FemmSession(model) as session:
            # The file deliberately stores zero current. This checks that the
            # post-processor uses the solve-time circuit current as well as the
            # model-level series-circuit identity.
            session.setCircuit("a", "current", 2.0)
            session.solve()
            return session.getcircuitprops("a")

    props = solve(SERIES_MODEL)
    np.testing.assert_allclose(props[0], 2.0, rtol=0, atol=1e-12)
    assert abs(props[2]) > 1e-12

    # Reversing every winding sign reverses the solved field but not the
    # circuit's self-linkage. Before issue #67 was fixed, getcircuitprops()
    # instead returned the last expanded block's +/-200 A-turn current and a
    # linkage whose sign depended on winding-label order.
    reversed_model = tmp_path / "series_circuit_reversed.fem"
    text = SERIES_MODEL.read_text()
    text = text.replace(
        "0.5\t0.5\t1\t0.2\t1\t0\t0\t100\t0",
        "0.5\t0.5\t1\t0.2\t1\t0\t0\t-100\t0",
    )
    text = text.replace(
        "2.5\t0.5\t1\t0.2\t1\t0\t0\t-100\t0",
        "2.5\t0.5\t1\t0.2\t1\t0\t0\t100\t0",
    )
    reversed_model.write_text(text)

    reversed_props = solve(reversed_model)
    np.testing.assert_allclose(reversed_props[0], 2.0, rtol=0, atol=1e-12)
    np.testing.assert_allclose(reversed_props[2], props[2], rtol=1e-6, atol=1e-12)
