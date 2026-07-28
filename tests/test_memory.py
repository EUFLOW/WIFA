"""Regression guard for site-data load memory.

Keeping an included ``wind_resource.nc`` as numpy arrays (windIO's
``nc_data="array"``) must stay dramatically cheaper than the dict-of-lists
default for a large resource.  See tests/mem_bench.py for the full breakdown.
"""

import gc
import tracemalloc

import numpy as np
import pytest
import windIO
import xarray as xr

pytest.importorskip("netCDF4")


def _write_resource(tmp_path, n_times, n_turbines):
    rng = np.random.default_rng(0)
    ds = xr.Dataset(
        {
            v: (("time", "wind_turbine"), rng.random((n_times, n_turbines)))
            for v in ("wind_speed", "wind_direction", "turbulence_intensity", "density")
        },
        coords={"time": np.arange(n_times), "wind_turbine": np.arange(n_turbines)},
    )
    ds["operating"] = (
        ("time", "wind_turbine"),
        np.ones((n_times, n_turbines), dtype=int),
    )
    ds.to_netcdf(tmp_path / "wind_resource.nc")
    sysf = tmp_path / "system.yaml"
    sysf.write_text("wind_resource: !include wind_resource.nc\n")
    return sysf


def _peak_load(sysf, mode):
    gc.collect()
    tracemalloc.start()
    tracemalloc.clear_traces()
    obj = windIO.load_yaml(sysf, nc_data=mode)
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    del obj
    gc.collect()
    return peak


def test_array_loader_reduces_peak_memory(tmp_path):
    """Array-backed load must peak well below the dict-of-lists load."""
    sysf = _write_resource(tmp_path, 2000, 60)

    list_peak = _peak_load(sysf, "list")
    array_peak = _peak_load(sysf, "array")

    # Measured ~4x; require a comfortable >=2x to avoid flakiness.
    assert array_peak < 0.5 * list_peak, (
        f"array load peak {array_peak / 1e6:.1f}MB not < 0.5 x list "
        f"{list_peak / 1e6:.1f}MB"
    )


def test_array_loader_keeps_ndarrays(tmp_path):
    """The array loader must actually keep numpy arrays (not lists)."""
    sysf = _write_resource(tmp_path, 100, 4)
    wr = windIO.load_yaml(sysf, nc_data="array")["wind_resource"]
    assert isinstance(wr["wind_speed"]["data"], np.ndarray)
    # default stays list-backed (backwards compatible)
    wr_list = windIO.load_yaml(sysf)["wind_resource"]
    assert isinstance(wr_list["wind_speed"]["data"], list)
