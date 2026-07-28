"""Memory benchmark for WIFA site-data construction.

Measures the peak Python allocation of turning the windIO ``wind_resource``
dict into a pyWake site, isolating the copies this optimization targets:
  * the dict-of-lists footprint (what windIO's _ds2yml currently produces)
  * the extra peak allocated inside ``construct_site`` (np.array / [cases_idx]
    / dict_to_netcdf copies)
  * the Weibull path's duplicate ``dict_to_netcdf`` build

Run:  pixi run python tests/mem_bench.py
"""

import gc
import tempfile
import tracemalloc
from pathlib import Path

import numpy as np
import xarray as xr

from wifa.pywake_api import construct_site, create_turbines

_TURBINE = {
    "name": "generic",
    "hub_height": 100.0,
    "rotor_diameter": 120.0,
    "performance": {
        "power_curve": {
            "power_wind_speeds": [3.0, 8.0, 12.0, 25.0],
            "power_values": [0.0, 1.0e6, 3.0e6, 3.0e6],
        },
        "Ct_curve": {
            "Ct_wind_speeds": [3.0, 8.0, 12.0, 25.0],
            "Ct_values": [0.8, 0.8, 0.4, 0.2],
        },
    },
}

_ANALYSIS = {
    "wind_deficit_model": {
        "name": "Jensen",
        "wake_expansion_coefficient": {"k_a": 0.04, "k_b": 0.0},
    },
    "deflection_model": {"name": "None"},
    "turbulence_model": {"name": "STF2005", "c1": 1.0, "c2": 1.0},
    "superposition_model": {
        "ws_superposition": "Linear",
        "ti_superposition": "Squared",
    },
    "rotor_averaging": {"name": "Center"},
    "blockage_model": {"name": "None"},
}


def _attrs():
    return {
        "flow_model": {"name": "pywake"},
        "analysis": _ANALYSIS,
        "model_outputs_specification": {
            "turbine_outputs": {
                "turbine_nc_filename": "turbine_data.nc",
                "output_variables": ["power"],
            }
        },
    }


def _layout(n_turbines, spacing=600.0):
    return {
        "coordinates": {
            "x": [i * spacing for i in range(n_turbines)],
            "y": [0.0] * n_turbines,
        }
    }


def _bounds(n_turbines, spacing=600.0):
    hi = n_turbines * spacing + 100
    return {"polygons": [{"x": [-100, hi, hi, -100], "y": [100, 100, -100, -100]}]}


def timeseries_system(n_times, n_turbines, arrays=False):
    """Per-turbine time-series system.

    ``arrays=False`` emits a dict-of-Python-lists (what windIO's _ds2yml
    produces today).  ``arrays=True`` emits a dict-of-ndarrays (what the
    Phase 2 ``to_dict(data="array")`` loader would produce) to preview the win.
    """
    rng = np.random.default_rng(0)
    _c = (lambda a: a) if arrays else (lambda a: a.tolist())
    ws = _c(8.0 + 4.0 * rng.random((n_times, n_turbines)))
    wd = _c(270.0 + 10.0 * rng.random((n_times, n_turbines)))
    ti = _c(0.06 + 0.04 * rng.random((n_times, n_turbines)))
    rho = _c(1.18 + 0.06 * rng.random((n_times, n_turbines)))
    op = _c(np.ones((n_times, n_turbines), dtype=int))
    dims = ["time", "wind_turbine"]
    return {
        "name": "bench timeseries",
        "site": {
            "name": "s",
            "boundaries": _bounds(n_turbines),
            "energy_resource": {
                "name": "r",
                "wind_resource": {
                    "time": list(range(n_times)),
                    "wind_turbine": list(range(n_turbines)),
                    "wind_speed": {"data": ws, "dims": dims},
                    "wind_direction": {"data": wd, "dims": dims},
                    "turbulence_intensity": {"data": ti, "dims": dims},
                    "density": {"data": rho, "dims": dims},
                    "operating": {"data": op, "dims": dims},
                },
            },
        },
        "wind_farm": {
            "name": "f",
            "layouts": [_layout(n_turbines)],
            "turbines": _TURBINE,
        },
        "attributes": _attrs(),
    }


def weibull_system(n_dirs, n_turbines):
    """Per-turbine Weibull system as a dict-of-lists (exercises the Weibull path)."""
    rng = np.random.default_rng(1)
    a = (8.0 + 2.0 * rng.random((n_turbines, n_dirs))).tolist()
    k = (2.0 + 0.3 * rng.random((n_turbines, n_dirs))).tolist()
    p = rng.random((n_turbines, n_dirs))
    p = (p / p.sum(axis=1, keepdims=True)).tolist()
    wd = np.linspace(0, 360, n_dirs, endpoint=False).tolist()
    dims = ["wind_turbine", "wind_direction"]
    return {
        "name": "bench weibull",
        "site": {
            "name": "s",
            "boundaries": _bounds(n_turbines),
            "energy_resource": {
                "name": "r",
                "wind_resource": {
                    "wind_direction": wd,
                    "wind_turbine": list(range(n_turbines)),
                    "weibull_a": {"data": a, "dims": dims},
                    "weibull_k": {"data": k, "dims": dims},
                    "sector_probability": {"data": p, "dims": dims},
                    "turbulence_intensity": {
                        "data": (0.07 * np.ones((n_turbines, n_dirs))).tolist(),
                        "dims": dims,
                    },
                },
            },
        },
        "wind_farm": {
            "name": "f",
            "layouts": [_layout(n_turbines)],
            "turbines": _TURBINE,
        },
        "attributes": _attrs(),
    }


def _peak_of(fn):
    tracemalloc.start()
    tracemalloc.clear_traces()
    obj = fn()
    cur, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return obj, cur, peak


def _construct_peak(system):
    farm = system["wind_farm"]
    _, _, hub_heights = create_turbines(farm)
    x = farm["layouts"][0]["coordinates"]["x"]
    resource = system["site"]["energy_resource"]
    tracemalloc.start()
    tracemalloc.reset_peak()
    base, _ = tracemalloc.get_traced_memory()
    construct_site(system, resource, hub_heights, x)
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return peak - base


def nc_load_peaks(n_times, n_turbines):
    """Peak memory of windIO.load_yaml on a real wind_resource.nc, list vs array.

    Requires a windIO build with the ``nc_data`` loader option; returns None if
    unavailable so the benchmark still runs the in-process cases.
    """
    import windIO

    d = Path(tempfile.mkdtemp())
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
    ds.to_netcdf(d / "wind_resource.nc")
    sysf = d / "system.yaml"
    sysf.write_text("wind_resource: !include wind_resource.nc\n")

    def peak(mode):
        gc.collect()
        tracemalloc.start()
        tracemalloc.clear_traces()
        obj = windIO.load_yaml(sysf, nc_data=mode)
        _, pk = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        del obj
        gc.collect()
        return pk

    try:
        return peak("list"), peak("array")
    except TypeError:
        return None  # windIO without the nc_data option


def main():
    MB = 1e6
    print(f"{'case':<28}{'input dict':>14}{'construct_site':>18}")
    print("-" * 60)

    n_times, n_wt = 4000, 100
    nbytes = n_times * n_wt * 8 * 5 / MB  # 5 float arrays
    sysd, _, dpeak = _peak_of(lambda: timeseries_system(n_times, n_wt))
    cs = _construct_peak(sysd)
    print(
        f"{'timeseries 4000x100 (lists)':<28}{dpeak / MB:>11.1f}MB{cs / MB:>15.1f}MB"
        f"   (raw ndarray ~{nbytes:.1f}MB)"
    )

    sysa, _, apeak = _peak_of(lambda: timeseries_system(n_times, n_wt, arrays=True))
    csa = _construct_peak(sysa)
    print(
        f"{'timeseries 4000x100 (arrays)':<28}{apeak / MB:>11.1f}MB{csa / MB:>15.1f}MB"
        f"   <- Phase 2 preview"
    )

    n_dirs, n_wt_w = 12, 2000
    sysw, _, wpeak = _peak_of(lambda: weibull_system(n_dirs, n_wt_w))
    csw = _construct_peak(sysw)
    print(f"{'weibull 2000wt x 12dir':<28}{wpeak / MB:>11.1f}MB{csw / MB:>15.1f}MB")

    print()
    peaks = nc_load_peaks(n_times, n_wt)
    if peaks is None:
        print("load_yaml(.nc): windIO has no nc_data option (skipped)")
    else:
        pl, pa = peaks
        print(
            f"load_yaml(.nc) 4000x100   list={pl / MB:.1f}MB  "
            f"array={pa / MB:.1f}MB  ({pl / pa:.1f}x lower peak)"
        )


if __name__ == "__main__":
    main()
