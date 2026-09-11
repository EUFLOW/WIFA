import os
from pathlib import Path

import numpy as np
import pytest

from wifa.windio_contract import TrackedDict, report_unread

test_path = Path(os.path.dirname(__file__))


def test_tracked_dict_records_reads():
    data = TrackedDict(
        {
            "a": 1,
            "nested": {"used": 2, "ignored": 3},
            "farms": [{"x": 1}, {"y": 2}],
            "untouched": {"deep": 4},
        }
    )

    _ = data["a"]
    _ = data["nested"]["used"]
    _ = data["farms"][0]["x"]

    unread = sorted(data.unread_paths())
    assert unread == [
        "wind_energy_system.farms[1].y",
        "wind_energy_system.nested.ignored",
        "wind_energy_system.untouched",
    ]


def test_tracked_dict_is_a_dict_and_iteration_counts_as_read():
    data = TrackedDict({"resource": {"ws": [1, 2], "wd": [3, 4]}})
    assert isinstance(data, dict)
    assert isinstance(data["resource"], dict)

    # Consuming a mapping wholesale (as windIO's dict_to_netcdf does)
    # counts as reading all of its keys
    dict(data["resource"].items())
    assert data.unread_paths() == []


def test_tracked_dict_get_and_pop():
    data = TrackedDict({"a": 1, "b": 2})
    assert data.get("a") == 1
    assert data.get("missing", "default") == "default"
    assert data.pop("b") == 2
    assert data.unread_paths() == []


def test_report_unread_warns():
    data = TrackedDict({"used": 1, "ignored": 2})
    _ = data["used"]
    with pytest.warns(UserWarning, match="ignored"):
        unread = report_unread(data, "testmodel")
    assert unread == ["wind_energy_system.ignored"]

    data = TrackedDict({"used": 1})
    _ = data["used"]
    assert report_unread(data, "testmodel") == []


def test_run_api_reports_unused_keys(tmp_path, monkeypatch):
    """End to end: an input key no runner consumes must be named in the
    warning, and the simulation result must be unaffected."""
    pytest.importorskip("py_wake")
    from wifa.main_api import run_api

    yaml_input = (
        test_path
        / "../examples/cases/multiple_wind_farms/wind_energy_system/system.yaml"
    )
    monkeypatch.chdir(tmp_path)

    with pytest.warns(UserWarning, match="electrical_substations"):
        results = run_api(str(yaml_input))
    assert isinstance(results, list) and len(results) == 3
    assert all(np.isfinite(a) and a > 0 for a in results)
