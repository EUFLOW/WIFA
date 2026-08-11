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

    # compute AEP with PyWake
    res_aep = (
        wfm(x, y, ws=np.arange(2, 30, 1), wd=dat["wd"])
        .aep(normalize_probabilities=True)
        .sum()
    )

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


def _two_farm_system_dict():
    """Tiny two-farm config sharing one turbine type and one wind resource."""
    from conftest import _ANALYSIS, _TURBINE

    farm_a = {
        "name": "Farm A (upwind)",
        "layouts": [{"coordinates": {"x": [0.0, 700.0], "y": [0.0, 0.0]}}],
        "turbines": _TURBINE,
    }
    farm_b = {
        "name": "Farm B (downwind)",
        "layouts": [{"coordinates": {"x": [3000.0, 3700.0], "y": [0.0, 0.0]}}],
        "turbines": _TURBINE,
    }
    return {
        "name": "multi-farm test",
        "site": {
            "name": "Test site",
            "boundaries": {
                "polygons": [
                    {"x": [-100, 4000, 4000, -100], "y": [100, 100, -100, -100]}
                ]
            },
            "energy_resource": {
                "name": "Test resource",
                "wind_resource": {
                    "wind_direction": [270.0],
                    "wind_speed": list(range(4, 26)),
                    "weibull_a": {"data": [9.0], "dims": ["wind_direction"]},
                    "weibull_k": {"data": [2.2], "dims": ["wind_direction"]},
                    "sector_probability": {"data": [1.0], "dims": ["wind_direction"]},
                    "turbulence_intensity": {
                        "data": [0.07],
                        "dims": ["wind_direction"],
                    },
                },
            },
        },
        "wind_farm": [farm_a, farm_b],
        "attributes": {
            "flow_model": {"name": "pywake"},
            "analysis": _ANALYSIS,
            "model_outputs_specification": {},
        },
    }


def test_pywake_multifarm_neighbor_wakes(tmp_path):
    system = _two_farm_system_dict()

    # Single-farm reference: each farm in isolation (no neighbor wakes)
    solo_a = dict(system, wind_farm=system["wind_farm"][0])
    solo_b = dict(system, wind_farm=system["wind_farm"][1])

    aep_multi = run_pywake(system, output_dir=str(tmp_path / "multi"))
    aep_a_solo = run_pywake(solo_a, output_dir=str(tmp_path / "a"))
    aep_b_solo = run_pywake(solo_b, output_dir=str(tmp_path / "b"))

    # Returns one AEP per farm and writes outputs
    assert isinstance(aep_multi, list) and len(aep_multi) == 2
    assert (tmp_path / "multi" / "output.yaml").exists()

    # Neighbor effect: farm B (downwind) loses energy when A is present
    aep_a_multi, aep_b_multi = aep_multi
    assert aep_b_multi < aep_b_solo
    # Farm A is upwind of B, so its AEP should be ~unchanged
    np.testing.assert_allclose(aep_a_multi, aep_a_solo, rtol=1e-3)


def test_pywake_multifarm_conflicting_turbine_specs(tmp_path):
    system = _two_farm_system_dict()

    # Same turbine name, different spec, must be rejected
    import copy

    conflicting = copy.deepcopy(system["wind_farm"][1]["turbines"])
    conflicting["hub_height"] = conflicting["hub_height"] + 10.0
    system["wind_farm"][1]["turbines"] = conflicting

    with pytest.raises(ValueError, match="defined differently"):
        run_pywake(system, output_dir=str(tmp_path))


def test_pywake_multifarm_rowp_example(tmp_path):
    """Run the IEA 22MW reference offshore wind plant (three farms, shared site)."""
    yaml_input = (
        test_path
        / "../examples/cases/multiple_wind_farms/wind_energy_system/system.yaml"
    )
    aep = run_pywake(str(yaml_input), output_dir=str(tmp_path))

    assert isinstance(aep, list) and len(aep) == 3
    assert all(farm_aep > 0 for farm_aep in aep)
    assert (tmp_path / "output.yaml").exists()


def test_pywake_multifarm_sum_matches_total(tmp_path):
    """Per-farm AEPs must sum to the AEP of the merged single-farm run."""
    system = _two_farm_system_dict()
    per_farm = run_pywake(system, output_dir=str(tmp_path / "multi"))

    merged = _two_farm_system_dict()
    farm_a, farm_b = merged["wind_farm"]
    coords_a = farm_a["layouts"][0]["coordinates"]
    coords_b = farm_b["layouts"][0]["coordinates"]
    merged["wind_farm"] = {
        "name": "merged",
        "layouts": [
            {
                "coordinates": {
                    "x": coords_a["x"] + coords_b["x"],
                    "y": coords_a["y"] + coords_b["y"],
                }
            }
        ],
        "turbines": farm_a["turbines"],
    }
    total = run_pywake(merged, output_dir=str(tmp_path / "single"))

    npt.assert_allclose(sum(per_farm), total, rtol=1e-6)


def test_pywake_multifarm_timeseries_sum_matches_total(tmp_path):
    """Same sum-to-total invariant on the time-series resource path."""
    from conftest import make_timeseries_per_turbine_system_dict

    base = make_timeseries_per_turbine_system_dict("pywake")
    total = run_pywake(base, output_dir=str(tmp_path / "single"))

    multi = make_timeseries_per_turbine_system_dict("pywake")
    farm = multi["wind_farm"]
    coords = farm["layouts"][0]["coordinates"]
    turbine = farm["turbines"]
    multi["wind_farm"] = [
        {
            "name": "Farm A",
            "layouts": [{"coordinates": {"x": coords["x"][:2], "y": coords["y"][:2]}}],
            "turbines": turbine,
        },
        {
            "name": "Farm B",
            "layouts": [{"coordinates": {"x": coords["x"][2:], "y": coords["y"][2:]}}],
            "turbines": turbine,
        },
    ]
    per_farm = run_pywake(multi, output_dir=str(tmp_path / "multi"))

    assert isinstance(per_farm, list) and len(per_farm) == 2
    npt.assert_allclose(sum(per_farm), total, rtol=1e-6)


def test_pywake_singlefarm_type_map_without_layout_types(tmp_path):
    """A single-spec turbine_types mapping without per-position layout types
    must fall back to that type (regression: KeyError on schema-valid input)."""
    system = _two_farm_system_dict()
    farm = system["wind_farm"][0]
    farm["turbine_types"] = {1: farm.pop("turbines")}
    system["wind_farm"] = farm

    aep = run_pywake(system, output_dir=str(tmp_path))
    assert np.isfinite(aep) and aep > 0


def test_pywake_multifarm_dict_form_layouts(tmp_path):
    """windIO allows `layouts` as a single mapping instead of a list."""
    system = _two_farm_system_dict()
    for farm in system["wind_farm"]:
        farm["layouts"] = farm["layouts"][0]

    per_farm = run_pywake(system, output_dir=str(tmp_path))
    assert isinstance(per_farm, list) and len(per_farm) == 2
    assert all(np.isfinite(a) and a > 0 for a in per_farm)


def test_pywake_layout_types_length_mismatch(tmp_path):
    """Fewer layout turbine_types entries than coordinates must raise, not
    silently misattribute turbines between farms."""
    system = _two_farm_system_dict()
    farm = system["wind_farm"][0]
    farm["turbine_types"] = {1: farm.pop("turbines")}
    farm["layouts"][0]["turbine_types"] = [1]  # 1 entry, 2 coordinates

    with pytest.raises(ValueError, match="turbine positions but"):
        run_pywake(system, output_dir=str(tmp_path))


def test_pywake_multitype_one_based_keys(tmp_path):
    """windIO turbine_types mappings may use arbitrary (e.g. 1-based) keys;
    layout entries are keys, not positional indices."""
    yaml_input = (
        test_path
        / "../examples/cases/windio_4turbines_multipleTurbines/wind_energy_system/system.yaml"
    )
    aep = run_pywake(str(yaml_input), output_dir=str(tmp_path))
    assert np.isfinite(aep) and aep > 0


def test_pywake_timeseries_two_hub_heights_without_ti(tmp_path):
    """Multi-hub-height time-series without turbulence_intensity must fall back
    to a default TI instead of crashing with site=None."""
    import copy

    from conftest import make_timeseries_per_turbine_system_dict

    system = make_timeseries_per_turbine_system_dict("pywake")
    del system["site"]["energy_resource"]["wind_resource"]["turbulence_intensity"]

    farm = system["wind_farm"]
    short = farm.pop("turbines")
    tall = copy.deepcopy(short)
    tall["name"] = short["name"] + " tall"
    tall["hub_height"] = short["hub_height"] + 20.0
    farm["turbine_types"] = {1: short, 2: tall}
    farm["layouts"][0]["turbine_types"] = [1, 1, 2]

    aep = run_pywake(system, output_dir=str(tmp_path))
    assert np.isfinite(aep) and aep > 0


def test_pywake_layout_types_positional_indices(tmp_path):
    """Layout turbine_types entries that are not mapping keys but are valid
    0-based indices must be interpreted positionally (windIO schema wording)."""
    system = _two_farm_system_dict()
    farm = system["wind_farm"][0]
    turbine = farm.pop("turbines")
    import copy

    other = copy.deepcopy(turbine)
    other["name"] = turbine["name"] + " b"
    farm["turbine_types"] = {"typeA": turbine, "typeB": other}
    farm["layouts"][0]["turbine_types"] = [0, 1]
    system["wind_farm"] = farm

    with pytest.warns(UserWarning, match="positional indices"):
        aep = run_pywake(system, output_dir=str(tmp_path))
    assert np.isfinite(aep) and aep > 0


def test_pywake_multifarm_shared_spec_with_nan(tmp_path):
    """Identical turbine specs containing NaN must merge, not be rejected as
    conflicting."""
    system = _two_farm_system_dict()
    for farm in system["wind_farm"]:
        farm["turbines"]["performance"]["cutout_wind_speed"] = float("nan")

    per_farm = run_pywake(system, output_dir=str(tmp_path))
    assert isinstance(per_farm, list) and len(per_farm) == 2


def test_specs_equal_edge_cases():
    """_specs_equal must handle NaN, numpy arrays, and nested list-of-dict
    specs without raising, and still detect genuine differences."""
    from wifa.pywake_api import _specs_equal

    nan_spec = {"a": float("nan"), "curve": [1.0, float("nan")]}
    assert _specs_equal(nan_spec, {"a": float("nan"), "curve": [1.0, float("nan")]})

    np_spec = {"modes": [{"power": np.array([1.0, 2.0])}]}
    assert _specs_equal(np_spec, {"modes": [{"power": np.array([1.0, 2.0])}]})
    assert not _specs_equal(np_spec, {"modes": [{"power": np.array([1.0, 3.0])}]})

    assert not _specs_equal({"a": 1}, {"a": 2})
    assert not _specs_equal({"a": [1, 2]}, {"a": [1, 2, 3]})
    assert _specs_equal({"a": np.array([1, 2])}, {"a": [1, 2]})


def test_pywake_timeseries_height_coord_per_turbine_ti(tmp_path):
    """Height-coordinate wind data combined with per-turbine TI on the
    multi-hub-height path must run (regression: uncaught IndexError)."""
    import copy

    from conftest import make_timeseries_per_turbine_system_dict

    system = make_timeseries_per_turbine_system_dict("pywake")
    resource = system["site"]["energy_resource"]["wind_resource"]
    n_times = len(resource["time"])

    resource["height"] = [80.0, 140.0]
    resource["wind_speed"] = {
        "data": [[8.0 + 0.1 * t, 9.0 + 0.1 * t] for t in range(n_times)],
        "dims": ["time", "height"],
    }
    resource["wind_direction"] = {
        "data": [[270.0, 272.0] for _ in range(n_times)],
        "dims": ["time", "height"],
    }
    # turbulence_intensity keeps dims ["time", "wind_turbine"] with no height

    farm = system["wind_farm"]
    short = farm.pop("turbines")
    tall = copy.deepcopy(short)
    tall["name"] = short["name"] + " tall"
    tall["hub_height"] = short["hub_height"] + 20.0
    farm["turbine_types"] = {1: short, 2: tall}
    farm["layouts"][0]["turbine_types"] = [1, 1, 2]

    with pytest.warns(UserWarning, match="averaged across"):
        aep = run_pywake(system, output_dir=str(tmp_path))
    assert np.isfinite(aep) and aep > 0


def test_pywake_timeseries_subset_multi_hub_height(tmp_path):
    """A times_run subset without an operating array must size the default
    operating array to the subset, not the full mask length."""
    import copy

    from conftest import make_timeseries_per_turbine_system_dict

    system = make_timeseries_per_turbine_system_dict("pywake")
    resource = system["site"]["energy_resource"]["wind_resource"]
    del resource["operating"]
    n_times = len(resource["time"])
    subset = [True, False, True, False, True, False][:n_times]
    system["attributes"]["model_outputs_specification"]["run_configuration"] = {
        "times_run": {"all_occurences": False, "subset": subset}
    }

    farm = system["wind_farm"]
    short = farm.pop("turbines")
    tall = copy.deepcopy(short)
    tall["name"] = short["name"] + " tall"
    tall["hub_height"] = short["hub_height"] + 20.0
    farm["turbine_types"] = {1: short, 2: tall}
    farm["layouts"][0]["turbine_types"] = [1, 1, 2]

    aep = run_pywake(system, output_dir=str(tmp_path))
    assert np.isfinite(aep) and aep > 0


def test_pywake_timeseries_turbine_first_dims(tmp_path):
    """Resource variables declared with dims ['wind_turbine', 'time'] must be
    averaged/subset along the declared axes, not hard-coded ones."""
    import copy

    from conftest import make_timeseries_per_turbine_system_dict

    system = make_timeseries_per_turbine_system_dict("pywake")
    resource = system["site"]["energy_resource"]["wind_resource"]
    for var in ("wind_speed", "wind_direction", "turbulence_intensity", "density"):
        data = np.array(resource[var]["data"])
        resource[var] = {"data": data.T.tolist(), "dims": ["wind_turbine", "time"]}
    del resource["operating"]

    farm = system["wind_farm"]
    short = farm.pop("turbines")
    tall = copy.deepcopy(short)
    tall["name"] = short["name"] + " tall"
    tall["hub_height"] = short["hub_height"] + 20.0
    farm["turbine_types"] = {1: short, 2: tall}
    farm["layouts"][0]["turbine_types"] = [1, 1, 2]

    with pytest.warns(UserWarning, match="averaged across"):
        aep = run_pywake(system, output_dir=str(tmp_path))
    assert np.isfinite(aep) and aep > 0


def test_farm_turbine_specs_resolution_rules():
    """Key resolution: exact keys always win over string coercion; mixed
    key/positional resolution raises; booleans are not positional indices."""
    from wifa.pywake_api import _farm_turbine_specs

    layout = {"coordinates": {"x": [0.0, 1.0], "y": [0.0, 0.0]}}
    spec_a, spec_b = {"name": "A"}, {"name": "B"}

    # keys 1 and "1" are distinct; each entry must resolve to its exact key
    farm = {
        "name": "f",
        "layouts": [dict(layout, turbine_types=["1", 1])],
        "turbine_types": {1: spec_a, "1": spec_b},
    }
    specs, per_pos = _farm_turbine_specs(farm)
    assert [specs[i]["name"] for i in per_pos] == ["B", "A"]

    # not all entries match keys, but all are valid 0-based indices:
    # positional interpretation with a loud per-entry warning (main-branch
    # semantics; supports 1-based keys used with 0-based indices)
    farm = {
        "name": "f",
        "layouts": [dict(layout, turbine_types=[0, 1])],
        "turbine_types": {2: spec_a, 0: spec_b},
    }
    with pytest.warns(UserWarning, match="positional indices"):
        specs, per_pos = _farm_turbine_specs(farm)
    assert [specs[i]["name"] for i in per_pos] == ["A", "B"]

    # 1-based integer keys with 0-based positional entries resolve positionally
    farm = {
        "name": "f",
        "layouts": [dict(layout, turbine_types=[0, 1])],
        "turbine_types": {1: spec_a, 2: spec_b},
    }
    with pytest.warns(UserWarning, match="positional indices"):
        specs, per_pos = _farm_turbine_specs(farm)
    assert [specs[i]["name"] for i in per_pos] == ["A", "B"]

    # integer entries resolve to quoted string keys via coercion
    farm = {
        "name": "f",
        "layouts": [dict(layout, turbine_types=[1, 2])],
        "turbine_types": {"1": spec_a, "2": spec_b},
    }
    specs, per_pos = _farm_turbine_specs(farm)
    assert [specs[i]["name"] for i in per_pos] == ["A", "B"]

    # booleans are rejected outright, even when they would alias integer
    # keys 0/1 through dict hashing
    for type_map in ({"a": spec_a, "b": spec_b}, {0: spec_a, 1: spec_b}):
        farm = {
            "name": "f",
            "layouts": [dict(layout, turbine_types=[True, False])],
            "turbine_types": type_map,
        }
        with pytest.raises(ValueError, match="boolean"):
            _farm_turbine_specs(farm)

    # entries matching nothing and not valid indices raise
    farm = {
        "name": "f",
        "layouts": [dict(layout, turbine_types=[0, 5])],
        "turbine_types": {"a": spec_a, "b": spec_b},
    }
    with pytest.raises(ValueError, match="neither match"):
        _farm_turbine_specs(farm)


def test_specs_equal_zero_dim_arrays():
    """0-d numpy arrays (e.g. values read from netCDF) must compare as
    scalars instead of raising 'iteration over a 0-d array'."""
    from wifa.pywake_api import _specs_equal

    assert _specs_equal(np.array(1.0), np.array(1.0))
    assert _specs_equal(np.array(1.0), 1.0)
    assert _specs_equal({"hh": np.array(119.0)}, {"hh": 119.0})
    assert not _specs_equal(np.array(1.0), 2.0)


def _with_height_profile_resource(system, n_times):
    """Replace ws/wd in a conftest system dict with height-profile data."""
    resource = system["site"]["energy_resource"]["wind_resource"]
    resource["height"] = [80.0, 140.0]
    resource["wind_speed"] = {
        "data": [[8.0 + 0.1 * t, 9.0 + 0.1 * t] for t in range(n_times)],
        "dims": ["time", "height"],
    }
    resource["wind_direction"] = {
        "data": [[270.0, 272.0] for _ in range(n_times)],
        "dims": ["time", "height"],
    }
    return resource


def _with_two_types(system, delta_hh=20.0):
    """Split the conftest single-type farm into two turbine types."""
    import copy

    farm = system["wind_farm"]
    short = farm.pop("turbines")
    tall = copy.deepcopy(short)
    tall["name"] = short["name"] + " tall"
    tall["hub_height"] = short["hub_height"] + delta_hh
    farm["turbine_types"] = {1: short, 2: tall}
    farm["layouts"][0]["turbine_types"] = [1, 1, 2]
    return system


def test_pywake_timeseries_ti_2d_without_dims(tmp_path):
    """2-D (time, height) TI data with no declared dims must be recognized as
    height-resolved when the resource declares a height coordinate."""
    from conftest import make_timeseries_per_turbine_system_dict

    system = make_timeseries_per_turbine_system_dict("pywake")
    n_times = len(system["site"]["energy_resource"]["wind_resource"]["time"])
    resource = _with_height_profile_resource(system, n_times)
    resource["turbulence_intensity"] = {
        "data": [[0.06, 0.08] for _ in range(n_times)]  # no "dims" key
    }
    del resource["operating"]
    _with_two_types(system)

    aep = run_pywake(system, output_dir=str(tmp_path))
    assert np.isfinite(aep) and aep > 0


def test_interp_helpers_floor_and_wraparound():
    """The extrapolation floor must not clamp in-range values, and direction
    interpolation must respect the 0/360 wraparound."""
    from wifa.pywake_api import _interp_along_height, _interp_direction_along_height

    heights = [80.0, 140.0]
    ti = np.array([[0.01, 0.01]])  # (time, height), legitimately below 0.02

    # in-range target: no clamping
    out = _interp_along_height(ti, ["time", "height"], heights, 110.0, min_val=0.02)
    npt.assert_allclose(out, 0.01)

    # out-of-range target: floor applies
    ws = np.array([[5.0, 8.0]])  # strong shear; extrapolates negative below
    out = _interp_along_height(ws, ["time", "height"], heights, 10.0, min_val=0.0)
    assert np.all(out >= 0.0)

    # 350 deg at 80 m and 10 deg at 140 m must interpolate near 0/360,
    # not to 180
    wd = np.array([[350.0, 10.0]])
    out = _interp_direction_along_height(wd, ["time", "height"], heights, 110.0)
    assert np.all((out >= 350.0) | (out <= 10.0))


def test_pywake_timeseries_density_with_height_dim(tmp_path):
    """Height-resolved density must be interpolated to hub height, not crash
    XRSite with a 2-D ('time',) assignment."""
    from conftest import make_timeseries_per_turbine_system_dict

    system = make_timeseries_per_turbine_system_dict("pywake")
    n_times = len(system["site"]["energy_resource"]["wind_resource"]["time"])
    resource = _with_height_profile_resource(system, n_times)
    resource["density"] = {
        "data": [[1.25, 1.22] for _ in range(n_times)],
        "dims": ["time", "height"],
    }
    del resource["operating"]
    _with_two_types(system)

    aep = run_pywake(system, output_dir=str(tmp_path))
    assert np.isfinite(aep) and aep > 0


def test_pywake_timeseries_ambiguous_undeclared_dims(tmp_path):
    """2-D data without declared dims whose axis lengths cannot distinguish
    height from turbine must raise, asking for explicit dims."""
    from conftest import make_timeseries_per_turbine_system_dict

    system = make_timeseries_per_turbine_system_dict("pywake")
    resource = system["site"]["energy_resource"]["wind_resource"]
    n_times = len(resource["time"])
    resource["height"] = [80.0, 140.0, 200.0]  # 3 heights == 3 turbines
    resource["wind_speed"] = {
        "data": [[8.0, 8.5, 9.0] for _ in range(n_times)]  # no dims key
    }
    del resource["operating"]
    _with_two_types(system)

    with pytest.raises(ValueError, match="declare its 'dims'"):
        run_pywake(system, output_dir=str(tmp_path))


def test_pywake_timeseries_height_dim_without_heights_coord(tmp_path):
    """A variable declaring a height dim while the resource has no height
    coordinate must raise a clear ValueError, not an opaque scipy error."""
    from conftest import make_timeseries_per_turbine_system_dict

    system = make_timeseries_per_turbine_system_dict("pywake")
    resource = system["site"]["energy_resource"]["wind_resource"]
    n_times = len(resource["time"])
    resource["turbulence_intensity"] = {
        "data": [[0.06, 0.08] for _ in range(n_times)],
        "dims": ["time", "height"],
    }
    del resource["operating"]
    _with_two_types(system)

    with pytest.raises(ValueError, match="height"):
        run_pywake(system, output_dir=str(tmp_path))


def test_pywake_timeseries_height_first_dims_equivalent(tmp_path):
    """dims ['height','time'] must give the same result as ['time','height'],
    including when a times_run subset selects exactly as many timesteps as
    there are heights (the axis-guessing poison case)."""
    from conftest import make_timeseries_per_turbine_system_dict

    def build(height_first):
        system = make_timeseries_per_turbine_system_dict("pywake")
        resource = _with_height_profile_resource(
            system, len(system["site"]["energy_resource"]["wind_resource"]["time"])
        )
        n_times = len(resource["time"])
        if height_first:
            for var in ("wind_speed", "wind_direction"):
                data = np.array(resource[var]["data"])
                resource[var] = {
                    "data": data.T.tolist(),
                    "dims": ["height", "time"],
                }
        del resource["operating"]
        # subset of exactly len(height)=2 timesteps
        subset = [True, True] + [False] * (n_times - 2)
        system["attributes"]["model_outputs_specification"]["run_configuration"] = {
            "times_run": {"all_occurences": False, "subset": subset}
        }
        return _with_two_types(system)

    aep_th = run_pywake(build(False), output_dir=str(tmp_path / "th"))
    aep_ht = run_pywake(build(True), output_dir=str(tmp_path / "ht"))
    npt.assert_allclose(aep_ht, aep_th, rtol=1e-9)


def test_pywake_timeseries_hub_height_outside_resource_heights(tmp_path):
    """A hub height outside the resource height range must extrapolate
    (consistently with the multi-type branch) instead of crashing."""
    from conftest import make_timeseries_per_turbine_system_dict

    system = make_timeseries_per_turbine_system_dict("pywake")
    resource = _with_height_profile_resource(
        system, len(system["site"]["energy_resource"]["wind_resource"]["time"])
    )
    resource["height"] = [150.0, 200.0]  # hub height 119 m is below the range
    n_times = len(resource["time"])
    resource["turbulence_intensity"] = {
        "data": [[0.06, 0.08] for _ in range(n_times)],
        "dims": ["time", "height"],
    }
    del resource["operating"]

    aep = run_pywake(system, output_dir=str(tmp_path))
    assert np.isfinite(aep) and aep > 0


def test_pywake_timeseries_two_types_same_hub_height(tmp_path):
    """Two turbine types sharing one hub height on the time-series path must
    not crash on the height-deduplication (regression: xarray dim conflict)."""
    import copy

    from conftest import make_timeseries_per_turbine_system_dict

    system = make_timeseries_per_turbine_system_dict("pywake")
    farm = system["wind_farm"]
    type_a = farm.pop("turbines")
    type_b = copy.deepcopy(type_a)
    type_b["name"] = type_a["name"] + " b"
    farm["turbine_types"] = {1: type_a, 2: type_b}
    farm["layouts"][0]["turbine_types"] = [1, 1, 2]

    aep = run_pywake(system, output_dir=str(tmp_path))
    assert np.isfinite(aep) and aep > 0


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


# if __name__ == "__main__":
#    test_heterogeneous_wind_rose()
#     simple_yaml_to_pywake('../examples/cases/windio_4turbines_multipleTurbines/plant_energy_turbine/IEA_10MW_turbine.yaml')
#    test_simple_wind_rose()
#    test_pywake_4wts_operating_flag()
#    test_pywake_4wts()
#    test_pywake_KUL()
