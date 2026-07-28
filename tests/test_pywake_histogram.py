"""Tests for the histogram (wind_direction x wind_speed probability) site path."""

import os
import shutil
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip(
    "py_wake", reason="py_wake not installed, install with: pip install wifa[pywake]"
)

import xarray as xr
from py_wake.deficit_models.gaussian import BastankhahGaussian
from py_wake.examples.data.dtu10mw._dtu10mw import DTU10MW
from py_wake.tests import npt

from wifa.pywake_api import (
    _subsector_wind_directions,
    construct_site,
    dict_to_site,
    run_pywake,
)

test_path = Path(os.path.dirname(__file__))

N_SECTORS = 12
WS_CENTERS = np.arange(0.5, 26.0, 1.0)  # 26 one-m/s bins, 0.5 .. 25.5
WD_SECTORS = np.arange(0.0, 360.0, 360.0 / N_SECTORS)


def _conditional_ws_pdf(mean_ws):
    """Discretized Rayleigh-like pdf on WS_CENTERS, normalized to sum to 1."""
    p = WS_CENTERS * np.exp(-((WS_CENTERS / mean_ws) ** 2))
    return p / p.sum()


def _histogram_resource(n_turbines=2, with_ti=True, uniform_wd=False):
    """Build a windIO wind_resource dict with a per-turbine 2D histogram.

    ``probability`` is the per-sector conditional over wind speed (windkit
    convention), ``sector_probability`` the direction marginal.
    """
    rng = np.random.default_rng(42)
    if uniform_wd:
        sector_p = np.full((n_turbines, N_SECTORS), 1.0 / N_SECTORS)
    else:
        raw = 1.0 + rng.random((n_turbines, N_SECTORS))
        sector_p = raw / raw.sum(axis=1, keepdims=True)

    cond = np.empty((n_turbines, N_SECTORS, WS_CENTERS.size))
    for t in range(n_turbines):
        for s in range(N_SECTORS):
            mean_ws = 8.0 + 0.5 * t if uniform_wd else 7.0 + t + 0.2 * s
            cond[t, s] = _conditional_ws_pdf(mean_ws)

    resource = {
        "wind_direction": WD_SECTORS.tolist(),
        "wind_speed": WS_CENTERS.tolist(),
        "wind_turbine": list(range(n_turbines)),
        "probability": {
            "dims": ["wind_turbine", "wind_direction", "wind_speed"],
            "data": cond.tolist(),
        },
        "sector_probability": {
            "dims": ["wind_turbine", "wind_direction"],
            "data": sector_p.tolist(),
        },
    }
    if with_ti:
        resource["turbulence_intensity"] = {
            "dims": ["wind_turbine", "wind_direction"],
            "data": np.full((n_turbines, N_SECTORS), 0.08).tolist(),
        }
    return resource


# --- dict_to_site: joint P construction --------------------------------------


def test_dict_to_site_builds_joint_p():
    resource = _histogram_resource(n_turbines=2)
    site = dict_to_site(resource)

    assert "P" in site.ds
    assert "probability" not in site.ds
    assert "Sector_frequency" not in site.ds

    # Joint = conditional x marginal, and sums to 1 per turbine
    cond = np.array(resource["probability"]["data"])
    sector_p = np.array(resource["sector_probability"]["data"])
    # XRSite appends a wd=360 wrap slice for circular interpolation; compare
    # on the original sectors.
    expected = cond * sector_p[:, :, np.newaxis]
    p = site.ds["P"].sel(wd=WD_SECTORS).transpose("i", "wd", "ws").values
    npt.assert_array_almost_equal(p, expected, 12)
    npt.assert_array_almost_equal(p.sum(axis=(1, 2)), np.ones(2), 12)

    # TI carried into the site dataset (run_simulation prefers it there)
    assert "TI" in site.ds.data_vars


def test_dict_to_site_bare_probability_is_joint():
    resource = _histogram_resource(n_turbines=2)
    cond = np.array(resource["probability"]["data"])
    sector_p = np.array(resource["sector_probability"]["data"])
    joint = cond * sector_p[:, :, np.newaxis]
    resource["probability"]["data"] = joint.tolist()
    del resource["sector_probability"]

    site = dict_to_site(resource)
    assert "P" in site.ds
    p = site.ds["P"].sel(wd=WD_SECTORS).transpose("i", "wd", "ws").values
    npt.assert_array_almost_equal(p, joint, 12)


def test_dict_to_site_windkit_dim_order():
    # windkit exports dims as (wind_speed, wind_direction, wind_turbine);
    # the joint-P construction must be dim-order independent.
    resource = _histogram_resource(n_turbines=2)
    cond = np.array(resource["probability"]["data"])
    resource["probability"] = {
        "dims": ["wind_speed", "wind_direction", "wind_turbine"],
        "data": cond.transpose(2, 1, 0).tolist(),
    }
    site = dict_to_site(resource)

    sector_p = np.array(resource["sector_probability"]["data"])
    expected = cond * sector_p[:, :, np.newaxis]
    p = site.ds["P"].sel(wd=WD_SECTORS).transpose("i", "wd", "ws").values
    npt.assert_array_almost_equal(p, expected, 12)


# --- construct_site dispatch --------------------------------------------------


def _system_dat():
    return {
        "attributes": {},
        "site": {"boundaries": {"polygons": [{"x": [0, 1000], "y": [0, 1000]}]}},
    }


def test_construct_site_dispatches_histogram():
    resource_dat = {"wind_resource": _histogram_resource(n_turbines=2)}
    site_data = construct_site(_system_dat(), resource_dat, {"0": 119.0}, [0.0, 500.0])

    assert site_data["timeseries"] is False
    # ws verbatim: pyWake cannot interpolate a ws-dependent P
    npt.assert_array_equal(site_data["ws"], WS_CENTERS)
    # wd refined to sub-sectors
    assert len(site_data["wd"]) == N_SECTORS * 5
    # operating must be (n_turbines, 1) to pass pyWake's arg2ilk validation
    assert site_data["operating"].shape == (2, 1)
    assert "P" in site_data["site"].ds


def test_histogram_ti_default_when_absent():
    resource_dat = {"wind_resource": _histogram_resource(n_turbines=2, with_ti=False)}
    site_data = construct_site(_system_dat(), resource_dat, {"0": 119.0}, [0, 500])

    assert site_data["TI"] == 0.06
    assert "TI" not in site_data["site"].ds.data_vars


# --- sub-sector helper --------------------------------------------------------


def test_subsector_wind_directions():
    wd = _subsector_wind_directions(WD_SECTORS, 5)
    assert len(wd) == N_SECTORS * 5
    assert np.all(np.diff(wd) > 0)
    npt.assert_array_almost_equal(np.diff(wd), np.full(N_SECTORS * 5 - 1, 6.0), 10)
    # Sector-0 sub-sectors wrap across north onto a 0-anchored 6° grid
    assert wd[0] == pytest.approx(0.0)
    assert wd[-1] == pytest.approx(354.0)

    # Trailing 360 wrap value is stripped before refinement
    wd_wrapped = _subsector_wind_directions(np.append(WD_SECTORS, 360.0), 5)
    npt.assert_array_almost_equal(wd_wrapped, wd, 10)

    # No refinement for n_subsector=1 or fewer than 4 sectors
    npt.assert_array_equal(_subsector_wind_directions(WD_SECTORS, 1), WD_SECTORS)
    npt.assert_array_equal(
        _subsector_wind_directions([0.0, 120.0, 240.0], 5), [0.0, 120.0, 240.0]
    )


# --- end-to-end through pyWake ------------------------------------------------


def _run_histogram_simulation(resource, x, y):
    """Mirror run_simulation()'s kwargs for a histogram site."""
    resource_dat = {"wind_resource": resource}
    site_data = construct_site(_system_dat(), resource_dat, {"0": 119.0}, x)
    # No ``operating`` here: that kwarg only exists for WIFA-built turbines
    # (PowerCtFunctionList wrapper); the run_pywake test covers it end-to-end.
    wfm = BastankhahGaussian(site_data["site"], DTU10MW())
    sim_res = wfm(
        x,
        y,
        type=0,
        time=site_data["timeseries"],
        ws=site_data["ws"],
        wd=site_data["wd"],
        yaw=0,
        tilt=0,
    )
    return sim_res


def test_histogram_single_turbine_aep_matches_analytic():
    # wd-uniform resource: the sub-sector interpolation of P is exact, so the
    # no-wake AEP must equal sum_ws pdf(ws) * power(ws) * 8760 h / 1e9.
    resource = _histogram_resource(n_turbines=1, uniform_wd=True)
    sim_res = _run_histogram_simulation(resource, [0.0], [0.0])
    aep = float(sim_res.aep(normalize_probabilities=True).sum())

    pdf = np.array(resource["probability"]["data"])[0, 0]  # same for every sector
    power_w = DTU10MW().power(WS_CENTERS)
    expected = float((pdf * power_w).sum() * 24 * 365 / 1e9)
    npt.assert_array_almost_equal(aep, expected, 6)


def test_histogram_two_turbine_wake_loss():
    # Two turbines aligned east-west, 3 diameters apart: with per-turbine P
    # (i dim) the simulation runs end-to-end and produces a wake loss.
    resource = _histogram_resource(n_turbines=2)
    sim_res = _run_histogram_simulation(resource, [0.0, 3 * 178.3], [0.0, 0.0])
    aep_per_turbine = (
        sim_res.aep(normalize_probabilities=True).sum(["ws", "wd"]).to_numpy()
    )

    assert np.all(np.isfinite(aep_per_turbine))
    assert np.all(aep_per_turbine > 0)
    total = float(sim_res.aep(normalize_probabilities=True).sum())
    npt.assert_array_almost_equal(total, aep_per_turbine.sum(), 8)
    # Wake interaction must cost something relative to two isolated turbines
    solo = [
        float(
            _run_histogram_simulation(
                {
                    **resource,
                    "wind_turbine": [0],
                    "probability": {
                        "dims": resource["probability"]["dims"],
                        "data": [resource["probability"]["data"][t]],
                    },
                    "sector_probability": {
                        "dims": resource["sector_probability"]["dims"],
                        "data": [resource["sector_probability"]["data"][t]],
                    },
                    "turbulence_intensity": {
                        "dims": resource["turbulence_intensity"]["dims"],
                        "data": [resource["turbulence_intensity"]["data"][t]],
                    },
                },
                [0.0],
                [0.0],
            )
            .aep(normalize_probabilities=True)
            .sum()
        )
        for t in range(2)
    ]
    assert total < sum(solo)


def test_run_pywake_histogram_case(tmp_path):
    """Full run_pywake() on a copy of the heterogeneous case with a histogram
    resource replacing the per-turbine Weibull."""
    src = test_path / "../examples/cases/heterogeneous_wind_rose_at_turbines"
    case = tmp_path / "histogram_case"
    shutil.copytree(src, case)

    # Swap the flow model to pywake
    system_yaml = case / "wind_energy_system/system.yaml"
    system_yaml.write_text(
        system_yaml.read_text().replace("name: foxes", "name: PyWake")
    )

    # Replace the Weibull resource with a per-turbine histogram, keeping the
    # original turbine positions/heights
    nc_path = case / "plant_energy_resource/WTResource.nc"
    with xr.open_dataset(nc_path) as old:
        old_vars = {v: old[v].values.copy() for v in ["x", "y", "height"]}
        n_turbines = old.sizes["wind_turbine"]

    resource = _histogram_resource(n_turbines=n_turbines)
    ds = xr.Dataset(
        {
            "probability": (
                tuple(resource["probability"]["dims"]),
                np.array(resource["probability"]["data"]),
            ),
            "sector_probability": (
                tuple(resource["sector_probability"]["dims"]),
                np.array(resource["sector_probability"]["data"]),
            ),
            "turbulence_intensity": (
                tuple(resource["turbulence_intensity"]["dims"]),
                np.array(resource["turbulence_intensity"]["data"]),
            ),
            **{name: (("wind_turbine",), vals) for name, vals in old_vars.items()},
        },
        coords={
            "wind_turbine": np.arange(n_turbines),
            "wind_direction": WD_SECTORS,
            "wind_speed": WS_CENTERS,
        },
    )
    nc_path.unlink()
    ds.to_netcdf(nc_path)

    aep = run_pywake(system_yaml, output_dir=str(tmp_path / "output"))
    assert np.isfinite(float(aep))
    assert float(aep) > 0
