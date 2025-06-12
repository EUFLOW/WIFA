import os
from pathlib import Path

import numpy as np
import pytest
import xarray as xr
from py_wake.deficit_models.gaussian import BastankhahGaussian
from py_wake.examples.data.dtu10mw._dtu10mw import DTU10MW
from py_wake.examples.data.hornsrev1 import Hornsrev1Site
from py_wake.rotor_avg_models import RotorCenter
from py_wake.site import XRSite
from py_wake.superposition_models import LinearSum
from py_wake.tests import npt
from py_wake.wind_turbines import WindTurbine
from py_wake.wind_turbines.power_ct_functions import (
    PowerCtFunctionList,
    PowerCtTabular,
)
from scipy.special import gamma
from windIO import __path__ as wiop
from windIO.utils.yml_utils import validate_yaml

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
    validate_yaml(yaml_input, windIO_path / Path("plant/wind_energy_system.yaml"))

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
    validate_yaml(yaml_input, windIO_path / Path("plant/wind_energy_system.yaml"))

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
    validate_yaml(yaml_input, windIO_path / Path("plant/wind_energy_system.yaml"))

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

def test_het_timeseries():
    resource_ds = xr.load_dataset('examples/cases/timeseries_with_turbine_specific_speeds/plant_energy_resource/Stochastic_atHubHeight.nc')
    reference_point_for_speed_up_calc = resource_ds.median(dim=('x', 'y'))
    resource_ds["Speedup"] = resource_ds['wind_speed'] / reference_point_for_speed_up_calc['wind_speed']
    
    resource_ds['P'] = 1
    resource_ds = resource_ds.rename({'wind_speed': 'WS', 'wind_direction': 'WD', 'turbulence_intensity': 'TI'})
    site = XRSite(resource_ds)
    
    x = [0, 1248.1, 2496.2, 3744.3]
    y = [0, 0, 0, 0]
    turbine = DTU10MW()
    # site = Hornsrev1Site()
    wfm = BastankhahGaussian(
        site,
        turbine,
        k=0.04,
        use_effective_ws=True,
        superpositionModel=LinearSum(),
        rotorAvgModel=RotorCenter(),
    )
    
    _f = np.zeros_like(resource_ds.time)
    
    sim_res = wfm(x, y, time=_f, ws=_f, wd=_f)
 
    wifa_res = run_pywake('examples/cases/timeseries_with_turbine_specific_speeds/wind_energy_system/system.yaml')

    assert(sim_res.aep() == wifa_res)

test_het_timeseries()
# if __name__ == "__main__":
#    test_heterogeneous_wind_rose()
#     simple_yaml_to_pywake('../examples/cases/windio_4turbines_multipleTurbines/plant_energy_turbine/IEA_10MW_turbine.yaml')
#    test_simple_wind_rose()
#    test_pywake_4wts_operating_flag()
#    test_pywake_4wts()
#    test_pywake_KUL()
