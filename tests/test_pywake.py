import os
import shutil
from pathlib import Path

import numpy as np
import pytest
import xarray as xr

pytest.importorskip(
    "py_wake", reason="py_wake not installed, install with: pip install wifa[pywake]"
)

from py_wake.deficit_models.gaussian import BastankhahGaussian
from py_wake.examples.data.dtu10mw._dtu10mw import DTU10MW
from py_wake.examples.data.hornsrev1 import Hornsrev1Site
from py_wake.rotor_avg_models import RotorCenter
from py_wake.site import XRSite
from py_wake.superposition_models import LinearSum
from py_wake.tests import npt
from py_wake.turbulence_models import CrespoHernandez
from py_wake.wind_turbines import WindTurbine
from py_wake.wind_turbines.power_ct_functions import PowerCtFunctionList, PowerCtTabular
from scipy.special import gamma
from windIO import __path__ as wiop
from windIO import validate as validate_yaml

from wifa.pywake_api import run_pywake

test_path = Path(os.path.dirname(__file__))
windIO_path = Path(wiop[0])
# sys.path.append(windIO.__path__[0])

# todo
# - set up KUL with constant thrust turbine / wake model
# - set up four turbines case with multiple turbine types


@pytest.fixture
def four_turbine_site(config_params):
    x = [0, 1248.1, 2496.2, 3744.3]
    y = [0, 0, 0, 0]
    # ws = [10.09, 8.798, 10.31]
    # wd = [271.8, 268.7, 271.1]
    config_name, ws, wd = config_params
    turbine = DTU10MW()
    site = Hornsrev1Site()
    # deficit = BastankhahGaussianDeficit()
    wfm = BastankhahGaussian(
        site,
        turbine,
        k=0.04,
        use_effective_ws=True,
        superpositionModel=LinearSum(),
        rotorAvgModel=RotorCenter(),
    )
    return wfm(x, y, ws=ws, wd=wd, time=True), config_name


def test_pywake_KUL():
    yaml_input = (
        test_path / "../examples/cases/KUL_LES/wind_energy_system/system_pywake.yaml"
    )

    # validate input
    validate_yaml(yaml_input, Path("plant/wind_energy_system"))

    # compute AEP (next step is to return a richer set of outputs)
    output_dir_name = "output_pywake_4wts"
    Path(output_dir_name).mkdir(parents=True, exist_ok=True)
    pywake_aep = run_pywake(yaml_input, output_dir=output_dir_name)
    # print(pywake_aep)

    # Check result
    pywake_aep_expected = 7515.2
    npt.assert_array_almost_equal(pywake_aep, pywake_aep_expected, 1)


@pytest.fixture(
    params=[
        # config_name, ws values, wd values
        ("windio_4turbines", [10.09, 8.798, 10.31], [271.8, 268.7, 271.1]),
        ("windio_4turbines_ABL", [10.09, 8.798, 10.31], [271.8, 268.7, 271.1]),
        ("windio_4turbines_ABL_stable", [10.09, 8.798, 10.31], [271.8, 268.7, 271.1]),
        (
            "windio_4turbines_profiles_stable",
            [9.708, 10.1, 11.25],
            [271.8, 268.7, 271.0],
        ),
    ]
)
def config_params(request):
    """Fixture that provides configuration parameters for each test case"""
    return request.param


def test_pywake_4wts(four_turbine_site):
    wfm, config_name = four_turbine_site

    yaml_input = (
        test_path / f"../examples/cases/{config_name}/wind_energy_system/system.yaml"
    )

    # validate input
    validate_yaml(yaml_input, Path("plant/wind_energy_system"))

    # compute AEP (next step is to return a richer set of outputs)
    output_dir_name = "output_pywake_4wts"
    Path(output_dir_name).mkdir(parents=True, exist_ok=True)
    pywake_aep = run_pywake(yaml_input, output_dir=output_dir_name)
    # print(pywake_aep)

    # Check result
    pywake_aep_expected = wfm.aep().sum()
    npt.assert_array_almost_equal(pywake_aep, pywake_aep_expected, 0)


def test_pywake_4wts_operating_flag():
    x = [0, 1248.1, 2496.2, 3744.3]
    y = [0, 0, 0, 0]
    ws = [10.0910225, 10.233016, 8.797999, 9.662098, 9.78371, 10.307792]
    wd = [271.82462, 266.20148, 268.6852, 273.61642, 263.45584, 271.05014]
    config_name = "timeseries_with_operating_flag"
    turbine = DTU10MW()
    turbine.powerCtFunction = PowerCtFunctionList(
        key="operating",
        powerCtFunction_lst=[
            PowerCtTabular(
                ws=[0, 100], power=[0, 0], power_unit="w", ct=[0, 0]
            ),  # 0=No power and ct
            turbine.powerCtFunction,
        ],  # 1=Normal operation
        default_value=1,
    )
    site = Hornsrev1Site()
    # deficit = BastankhahGaussianDeficit()
    wfm = BastankhahGaussian(
        site,
        turbine,
        k=0.04,
        use_effective_ws=True,
        superpositionModel=LinearSum(),
        rotorAvgModel=RotorCenter(),
    )

    operating = np.ones((len(ws), len(x)))
    operating[:-2, 0] = 0
    res = wfm(x, y, ws=ws, wd=wd, time=True, operating=operating.T)

    yaml_input = (
        test_path / f"../examples/cases/{config_name}/wind_energy_system/system.yaml"
    )

    # validate input
    validate_yaml(yaml_input, Path("plant/wind_energy_system"))

    # compute AEP (next step is to return a richer set of outputs)
    output_dir_name = "output_pywake_4wts"
    Path(output_dir_name).mkdir(parents=True, exist_ok=True)
    pywake_aep = run_pywake(yaml_input, output_dir=output_dir_name)
    # print(pywake_aep)

    # Check result
    pywake_aep_expected = res.aep().sum()
    npt.assert_array_almost_equal(pywake_aep, pywake_aep_expected, 0)


# fmt: off
POWER_CT_TABLE = PowerCtTabular(
    [
        3, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0,
        12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0,
        20.0, 21.0, 22.0, 23.0, 24.0, 25.0,
    ],
    [
        0, 263388.0, 751154.0, 1440738.0, 2355734.0, 3506858.0, 4993092.0,
        6849310.0, 9116402.0, 10000754.0, 10009590.0, 10000942.0, 10042678.0,
        10003480.0, 10001600.0, 10001506.0, 10013632.0, 10007428.0, 10005360.0,
        10002728.0, 10001130.0, 10004984.0, 9997558.0,
    ],
    "W",
    [
        0.923, 0.923, 0.919, 0.904, 0.858, 0.814, 0.814, 0.814, 0.814, 0.577,
        0.419, 0.323, 0.259, 0.211, 0.175, 0.148, 0.126, 0.109, 0.095, 0.084,
        0.074, 0.066, 0.059,
    ],
)
# fmt: on


def test_simple_wind_rose():
    _ = run_pywake(
        test_path / "../examples/cases/simple_wind_rose/wind_energy_system/system.yaml"
    )
    x = [0, 1248.1, 2496.2, 3744.3]
    y = [0, 0, 0, 0]
    site = Hornsrev1Site()
    turbine = WindTurbine(
        name="test",
        diameter=178.3,
        hub_height=119.0,
        powerCtFunction=POWER_CT_TABLE,
    )

    #  power_curve:
    #    power_values: [0, 263388., 751154., 1440738., 2355734., 3506858., 4993092., 6849310., 9116402., 10000754., 10009590., 10000942., 10042678., 10003480., 10001600., 10001506., 10013632., 10007428., 10005360., 10002728., 10001130., 10004984., 9997558.]
    #    power_wind_speeds: [3, 4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.,16.,17.,18.,19.,20.,21.,22.,23.,24.,25.]
    #  Ct_curve:
    #    Ct_values: [0.923, 0.923,0.919,0.904,0.858,0.814,0.814,0.814,0.814,0.577,0.419,0.323,0.259,0.211,0.175,0.148,0.126,0.109,0.095,0.084,0.074,0.066,0.059]
    #    Ct_wind_speeds: [3, 4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.,16.,17.,18.,19.,20.,21.,22.,23.,24.,25.]
    # hub_height: 119.0
    # rotor_diameter: 178.3

    # turbine = DTU10MW()
    wfm = BastankhahGaussian(
        site,
        turbine,
        k=0.04,
        ceps=0.2,
        use_effective_ws=True,
        superpositionModel=LinearSum(),
        rotorAvgModel=RotorCenter(),
    )
    res = wfm(x, y, wd=np.arange(0, 361, 30), TI=0.1)
    assert xr.load_dataset("output/PowerTable.nc").power.mean() == res.Power.mean()

    # def simple_yaml_to_pywake(ymlfile):
    #    dat = load_yaml(ymlfile)
    #    speeds = dat['cp_ws']
    #    pows = dat['power']
    #    cts = dat['ct']
    #
    #    hello
    #


def test_heterogeneous_wind_rose_grid():
    turbine = WindTurbine(
        name="test",
        diameter=178.3,
        hub_height=119.0,
        powerCtFunction=POWER_CT_TABLE,
    )
    dat = xr.load_dataset(
        test_path
        / "../examples/cases/heterogeneous_wind_rose_map/plant_energy_resource/Stochastic_atHubHeight.nc"
    )
    dat = dat.rename(
        {
            "wind_direction": "wd",
            "sector_probability": "Sector_frequency",
            "weibull_a": "Weibull_A",
            "weibull_k": "Weibull_k",
            "turbulence_intensity": "TI",
            "height": "h",
        }
    )
    mean_ws = dat["Weibull_A"].values * gamma(
        1 + 1.0 / dat["Weibull_k"].values
    )  # shape (x,y,h,wd)
    max_mean = np.max(mean_ws, axis=(0, 1))  # shape (h,wd)
    speedup = mean_ws / max_mean  # normalized speed-up (x,y,h,wd)
    dat["Speedup"] = (("x", "y", "h", "wd"), speedup)

    site = XRSite(dat)
    wfm = BastankhahGaussian(
        site,
        turbine,
        k=0.04,
        ceps=0.2,
        use_effective_ws=True,
        superpositionModel=LinearSum(),
        rotorAvgModel=RotorCenter(),
    )
    x = [0, 1248.1, 2496.2, 3744.3]
    y = [0, 0, 0, 0]

    # Compute Speedup-adjusted ws range (same logic as WIFA)
    A_vals = dat["Weibull_A"].values
    k_vals = dat["Weibull_k"].values
    ws_999 = A_vals * (-np.log(0.001)) ** (1.0 / k_vals)
    min_su = np.min(speedup)
    ws_max_ref = np.max(ws_999) / max(min_su, 0.1)
    ws_range = np.arange(0, np.ceil(ws_max_ref) + 0.5, 0.5)

    # Compute sub-sector wd (same logic as WIFA)
    wd_sectors = dat["wd"].values
    if len(wd_sectors) > 1 and np.isclose(wd_sectors[-1], 360.0):
        wd_sectors = wd_sectors[:-1]
    n_sub = 5
    sw = 360.0 / len(wd_sectors)
    ssw = sw / n_sub
    offsets = np.linspace(-sw / 2 + ssw / 2, sw / 2 - ssw / 2, n_sub)
    wd_fine = np.sort(
        (wd_sectors[:, np.newaxis] + offsets[np.newaxis, :]).ravel() % 360
    )

    # compute AEP with PyWake
    res_aep = wfm(x, y, ws=ws_range, wd=wd_fine).aep(normalize_probabilities=True).sum()

    # compute AEP with API
    wifa_res = run_pywake(
        test_path
        / "../examples/cases/heterogeneous_wind_rose_map/wind_energy_system/system.yaml"
    )

    # we need these to match
    assert wifa_res == res_aep


def test_heterogeneous_wind_rose_arbitrary_points():
    ds = xr.open_dataset(
        test_path
        / "../examples/cases/heterogeneous_wind_rose_at_turbines/plant_energy_resource/WTResource.nc"
    )
    ds = ds.rename(
        {
            "wind_direction": "wd",
            "wind_speed": "ws",
            "wind_turbine": "i",
            "sector_probability": "Sector_frequency",
            "weibull_a": "Weibull_A",
            "weibull_k": "Weibull_k",
            "turbulence_intensity": "TI",
            "height": "h",
        }
    )
    mean_ws = ds["Weibull_A"].values * gamma(
        1 + 1.0 / ds["Weibull_k"].values
    )  # shape (i,wd)
    max_mean = np.max(mean_ws, axis=0)  # shape (wd,)
    Speedup = mean_ws / max_mean  # normalized speed-up (i,wd)
    ds["Speedup"] = (("i", "wd"), Speedup)
    site = XRSite(ds)

    turbine = WindTurbine(
        name="test",
        diameter=178.3,
        hub_height=119.0,
        powerCtFunction=POWER_CT_TABLE,
    )
    wfm = BastankhahGaussian(
        site,
        turbine,
        k=0.04,
        ceps=0.2,
        use_effective_ws=True,
        superpositionModel=LinearSum(),
        rotorAvgModel=RotorCenter(),
    )
    x = ds["x"].values
    y = ds["y"].values

    res_aep = (
        wfm(x, y, ws=ds["ws"], wd=ds["wd"]).aep(normalize_probabilities=True).sum()
    )

    wifa_res = run_pywake(
        test_path
        / "../examples/cases/heterogeneous_wind_rose_at_turbines/wind_energy_system/system.yaml"
    )

    assert wifa_res == res_aep


def test_turbine_specific_speeds_timeseries():
    """
    Test case for time-series simulation where inflow wind speed/direction
    is specific to each turbine (dimensions: time x turbine).
    Validates that the API correctly reduces 2D inputs to 1D reference arrays
    for the simulation call while preserving local site data.
    """
    from windIO import load_yaml  # Import helper to read YAML

    case_name = "turbine_specific_speeds_timeseries"
    system_yaml = (
        test_path / f"../examples/cases/{case_name}/wind_energy_system/system.yaml"
    )
    resource_nc = (
        test_path
        / f"../examples/cases/{case_name}/plant_energy_resource/Stochastic_atHubHeight.nc"
    )

    # 1. Run via API
    wifa_res = run_pywake(system_yaml)

    # 2. Run manually to verify logic

    # Load coordinates from YAML
    sys_dat = load_yaml(system_yaml)
    farm_layout = sys_dat["wind_farm"]["layouts"]
    if isinstance(farm_layout, list):
        coords = farm_layout[0]["coordinates"]
    else:
        coords = farm_layout["coordinates"]
    x = coords["x"]
    y = coords["y"]

    # Load and prepare Site Data
    ds = xr.open_dataset(resource_nc)
    ds = ds.rename(
        {
            "wind_direction": "WD",
            "wind_speed": "WS",
            "wind_turbine": "i",
            "turbulence_intensity": "TI",
        }
    )

    # FIX 1: Re-index 'i' to be 0-based to match PyWake's internal numbering
    # The file has [1, 2, 3, 4], PyWake expects [0, 1, 2, 3]
    ds = ds.assign_coords(i=np.arange(len(ds.i)))

    # FIX 2: Transpose to (i, time) for XRSite linear interpolator
    ds = ds.transpose("i", "time")

    # FIX 3: Add uniform probability 'P'
    n_time = len(ds.time)
    ds["P"] = (("time"), np.ones(n_time) / n_time)

    # Initialize Site
    site = XRSite(ds, interp_method="linear")

    # Define Turbine
    turbine = WindTurbine(
        name="test",
        diameter=178.3,
        hub_height=119.0,
        powerCtFunction=POWER_CT_TABLE,
    )

    # Define Wake Model
    wfm = BastankhahGaussian(
        site,
        turbine,
        k=0.04,
        ceps=0.2,
        superpositionModel=LinearSum(),
        use_effective_ws=True,
        turbulenceModel=CrespoHernandez(),
    )

    # Calculate reference arrays (API Logic)
    if "i" in ds.WS.dims:
        ws_ref = ds.WS.mean(dim="i").values
    else:
        ws_ref = ds.WS.values

    if "i" in ds.WD.dims:
        rads = np.deg2rad(ds.WD)
        mean_sin = np.sin(rads).mean(dim="i")
        mean_cos = np.cos(rads).mean(dim="i")
        wd_ref = np.rad2deg(np.arctan2(mean_sin, mean_cos)) % 360
        wd_ref = wd_ref.values
    else:
        wd_ref = ds.WD.values

    # Manual Simulation
    res_manual = wfm(x, y, time=ds.time, ws=ws_ref, wd=wd_ref)

    manual_aep = res_manual.aep(normalize_probabilities=False).sum()

    # 3. Assert match
    npt.assert_allclose(wifa_res, manual_aep, rtol=1e-6)


def test_pywake_dict_timeseries_per_turbine_with_density(tmp_path):
    from conftest import make_timeseries_per_turbine_system_dict

    # Run with density
    system_dict = make_timeseries_per_turbine_system_dict("pywake")
    output_dir = tmp_path / "output_pywake_ts"
    aep_with = run_pywake(system_dict, output_dir=str(output_dir))
    assert np.isfinite(aep_with) and aep_with > 0

    # Run without density — same config but density removed
    system_dict_no = make_timeseries_per_turbine_system_dict("pywake")
    del system_dict_no["site"]["energy_resource"]["wind_resource"]["density"]
    output_dir_no = tmp_path / "output_pywake_ts_no_density"
    aep_without = run_pywake(system_dict_no, output_dir=str(output_dir_no))

    # Density correction should change AEP (test data varies around 1.225)
    assert aep_with != aep_without


def test_weibull_speedup_dim_ordering(tmp_path):
    """Regression test: per-turbine Weibull Speedup with both dim orderings.

    flow_model_chain (via windkit) writes wind_resource.nc with dims
    (wind_direction, wind_turbine), while WIFA's own test fixtures use
    (wind_turbine, wind_direction).  A bug in _construct_weibull_site()
    previously hardcoded axis=0 for the Speedup normalisation, which only
    worked for the turbine-first ordering.  With direction-first data the
    Speedup dims were silently swapped and PyWake ignored the variable,
    removing all terrain-induced wind speed inhomogeneity from the wake
    simulation and inflating wake losses from ~10 % to ~39 %.

    This test runs the same per-turbine Weibull case with BOTH dim
    orderings and asserts identical AEP.
    """
    from conftest import _ANALYSIS, _TURBINE

    n_wd = 4
    n_wt = 4
    wd_vals = [0.0, 90.0, 180.0, 270.0]
    ws_vals = list(np.arange(4.0, 26.0, 1.0).tolist())

    # Per-turbine, per-sector Weibull A — turbine 3 is windiest
    # Shape: (wind_direction, wind_turbine) = (4, 4)
    A_data = [
        [7.0, 8.0, 9.0, 10.0],  # sector 0°
        [6.5, 7.5, 8.5, 9.5],  # sector 90°
        [8.0, 9.0, 10.0, 11.0],  # sector 180°
        [6.0, 7.0, 8.0, 9.0],  # sector 270°
    ]
    k_data = [[2.0] * n_wt] * n_wd
    freq_data = [[1.0 / n_wd] * n_wt] * n_wd
    ti_data = [[0.06] * n_wt] * n_wd

    common_site = {
        "name": "Test site",
        "boundaries": {
            "polygons": [{"x": [-90, 5000, 5000, -90], "y": [90, 90, -90, -90]}]
        },
    }
    common_farm = {
        "name": "Test farm",
        "layouts": [
            {"coordinates": {"x": [0, 1248.1, 2496.2, 3744.3], "y": [0, 0, 0, 0]}}
        ],
        "turbines": _TURBINE,
    }
    common_attrs = {
        "flow_model": {"name": "pywake"},
        "analysis": _ANALYSIS,
        "model_outputs_specification": {
            "turbine_outputs": {
                "turbine_nc_filename": "PowerTable.nc",
                "output_variables": ["power"],
            },
        },
    }

    def _make_system(data, dims, name):
        return {
            "name": name,
            "site": {
                **common_site,
                "energy_resource": {
                    "name": "Test resource",
                    "wind_resource": {
                        "wind_direction": wd_vals,
                        "wind_speed": ws_vals,
                        "wind_turbine": list(range(n_wt)),
                        "reference_height": 119.0,
                        "weibull_a": {"data": data["A"], "dims": dims},
                        "weibull_k": {"data": data["k"], "dims": dims},
                        "sector_probability": {"data": data["f"], "dims": dims},
                        "turbulence_intensity": {"data": data["ti"], "dims": dims},
                    },
                },
            },
            "wind_farm": common_farm,
            "attributes": common_attrs,
        }

    # --- 1. Direction-first ordering (flow_model_chain convention) --------
    wd_first = _make_system(
        {"A": A_data, "k": k_data, "f": freq_data, "ti": ti_data},
        ["wind_direction", "wind_turbine"],
        "Direction-first",
    )
    aep_wd_first = run_pywake(wd_first, output_dir=str(tmp_path / "wd_first"))
    assert np.isfinite(aep_wd_first) and aep_wd_first > 0

    # --- 2. Turbine-first ordering (WIFA test-fixture convention) ---------
    A_T = np.array(A_data).T.tolist()
    k_T = np.array(k_data).T.tolist()
    freq_T = np.array(freq_data).T.tolist()
    ti_T = np.array(ti_data).T.tolist()

    wt_first = _make_system(
        {"A": A_T, "k": k_T, "f": freq_T, "ti": ti_T},
        ["wind_turbine", "wind_direction"],
        "Turbine-first",
    )
    aep_wt_first = run_pywake(wt_first, output_dir=str(tmp_path / "wt_first"))

    # Both orderings must produce identical AEP
    npt.assert_allclose(aep_wd_first, aep_wt_first, rtol=1e-6)


# if __name__ == "__main__":
#    test_heterogeneous_wind_rose()
#     simple_yaml_to_pywake('../examples/cases/windio_4turbines_multipleTurbines/plant_energy_turbine/IEA_10MW_turbine.yaml')
#    test_simple_wind_rose()
#    test_pywake_4wts_operating_flag()
#    test_pywake_4wts()
#    test_pywake_KUL()
