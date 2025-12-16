import argparse
import os
import warnings
from collections import OrderedDict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
import yaml
from scipy.interpolate import interp1d
from scipy.special import gamma
from windIO import dict_to_netcdf, load_yaml
from windIO import validate as validate_yaml

# Define default values for wind_deficit_model parameters
DEFAULTS = {
    "wind_deficit_model": {
        #'wind_deficit_model': {
        "name": "Jensen",
        #'k': 0.04,  # Default wake expansion coefficient for Jensen
        #'k2': 0.0,  # TI wake expansion modifier
    },
    "deflection_model": {
        "name": "Jimenez",
        "beta": 0.1,  # Default Jimenez deflection coefficient
    },
    "turbulence_model": {
        "name": "STF2005",
        "c1": 1.0,  # Default STF C1 value
        "c2": 1.0,  # Default STF C2 value
    },
    "superposition_model": {
        "ws_superposition": "Linear",
    },
    "rotor_averaging": {
        "name": "Center",
    },
    "blockage_model": {"name": None, "ss_alpha": 0.8888888888888888},
}


def get_with_default(data, key, defaults):
    """
    Retrieve a value from a dictionary, using a default if the key is not present.
    If the value is a dictionary, apply the same process recursively.
    """
    if key not in data:
        print("WARNING: Using default value for ", key)
        return defaults[key]
    elif isinstance(data[key], dict):
        # For nested dictionaries, ensure all subkeys are checked for defaults
        return {
            sub_key: get_with_default(data[key], sub_key, defaults[key])
            for sub_key in defaults[key]
        }
    else:
        return data[key]


def weighted_quantile(
    values, quantiles, sample_weight=None, values_sorted=False, old_style=False
):
    """Very close to numpy.percentile, but supports weights.
    NOTE: quantiles should be in [0, 1]!
    :param values: numpy.array with data
    :param quantiles: array-like with many quantiles needed
    :param sample_weight: array-like of the same length as `array`
    :param values_sorted: bool, if True, then will avoid sorting of
        initial array
    :param old_style: if True, will correct output to be consistent
        with numpy.percentile.
    :return: numpy.array with computed quantiles.
    """
    values = np.array(values)
    quantiles = np.array(quantiles)
    if sample_weight is None:
        sample_weight = np.ones(len(values))
    sample_weight = np.array(sample_weight)
    assert np.all(quantiles >= 0) and np.all(
        quantiles <= 1
    ), "quantiles should be in [0, 1]"

    if not values_sorted:
        sorter = np.argsort(values)
        values = values[sorter]
        sample_weight = sample_weight[sorter]

    weighted_quantiles = np.cumsum(sample_weight) - 0.5 * sample_weight
    if old_style:
        # To be convenient with numpy.percentile
        weighted_quantiles -= weighted_quantiles[0]
        weighted_quantiles /= weighted_quantiles[-1]
    else:
        weighted_quantiles /= np.suyamlFilem(sample_weight)
    return np.interp(quantiles, weighted_quantiles, values)


def run_pywake(yamlFile, output_dir="output"):
    from py_wake import NOJ, BastankhahGaussian, HorizontalGrid
    from py_wake.deficit_models import SelfSimilarityDeficit2020
    from py_wake.deficit_models.fuga import FugaDeficit
    from py_wake.deficit_models.gaussian import (
        BastankhahGaussianDeficit,
        BlondelSuperGaussianDeficit2020,
        TurboGaussianDeficit,
    )
    from py_wake.deficit_models.noj import NOJLocalDeficit
    from py_wake.deflection_models import JimenezWakeDeflection
    from py_wake.examples.data.hornsrev1 import Hornsrev1Site
    from py_wake.rotor_avg_models import (
        CGIRotorAvg,
        EqGridRotorAvg,
        GaussianOverlapAvgModel,
        GQGridRotorAvg,
        GridRotorAvg,
        PolarGridRotorAvg,
        RotorCenter,
        polar_gauss_quadrature,
    )
    from py_wake.site import XRSite
    from py_wake.superposition_models import LinearSum, SquaredSum
    from py_wake.turbulence_models import (
        CrespoHernandez,
        STF2005TurbulenceModel,
        STF2017TurbulenceModel,
    )
    from py_wake.wind_farm_models import All2AllIterative, PropagateDownwind
    from py_wake.wind_turbines import WindTurbine, WindTurbines
    from py_wake.wind_turbines.power_ct_functions import (
        PowerCtFunctionList,
        PowerCtFunctions,
        PowerCtTabular,
    )

    def dict_to_site(resource_dict):
        resource_ds = dict_to_netcdf(resource_dict)
        rename_map = {
            "height": "h",
            "weibull_a": "Weibull_A",
            "weibull_k": "Weibull_k",
            "sector_probability": "Sector_frequency",
            "turbulence_intensity": "TI",
            "wind_turbine": "i",
        }

        # Smart rename for wind_direction and wind_speed
        for key, coord_name, var_name in [
            ("wind_direction", "wd", "WD"),
            ("wind_speed", "ws", "WS"),
        ]:
            if key in resource_ds:
                # If it's a coordinate (dimension), use lowercase (wd, ws)
                # If it's a data variable (time series/map), use uppercase (WD, WS)
                rename_map[key] = coord_name if key in resource_ds.coords else var_name

        for name in rename_map:
            if name in resource_ds:
                resource_ds = resource_ds.rename({name: rename_map[name]})

        if "P" not in resource_ds and "time" in resource_ds.dims:
            n_time = len(resource_ds.time)
            # Create uniform probability array (1/N)
            resource_ds["P"] = (("time",), np.ones(n_time) / n_time)
        if "i" in resource_ds.dims:
            other_dims = [d for d in resource_ds.dims if d != "i"]
            # The transpose operation ensures that 'i' (turbine index) is the first dimension.
            # This is required for XRSite's linear interpolation, which expects the turbine index
            # as the leading dimension. See test at line 392 for details.
            resource_ds = resource_ds.transpose("i", *other_dims)
        print("making site with ", resource_ds)
        return XRSite(resource_ds)

    # allow yamlFile to be an already parsed input dict
    if not isinstance(yamlFile, dict):
        validate_yaml(yamlFile, "plant/wind_energy_system")
        system_dat = load_yaml(Path(yamlFile))
    else:
        system_dat = yamlFile

    # output_dir priority: 1) yaml file, 2) function argument, 3) default
    output_dir = str(
        system_dat["attributes"]
        .get("model_outputs_specification", {})
        .get("output_folder", output_dir)
    )

    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # check for multiple turbines?

    # define turbine
    farm_dat = system_dat["wind_farm"]
    if "turbines" in farm_dat:
        turbine_dats = [farm_dat["turbines"]]
        type_names = "0"
    else:
        turbine_dats = [
            farm_dat["turbine_types"][key] for key in farm_dat["turbine_types"]
        ]
        type_names = list(farm_dat["turbine_types"].keys())

    turbines = []
    hub_heights = {}
    for turbine_dat, key in zip(turbine_dats, type_names):
        hh = turbine_dat["hub_height"]
        rd = turbine_dat["rotor_diameter"]
        hub_heights[key] = hh
        if "Cp_curve" in turbine_dat["performance"]:
            cp = turbine_dat["performance"]["Cp_curve"]["Cp_values"]
            cp_ws = turbine_dat["performance"]["Cp_curve"]["Cp_wind_speeds"]
            power_curve_type = "cp"
        elif "power_curve" in turbine_dat["performance"]:
            cp_ws = turbine_dat["performance"]["power_curve"]["power_wind_speeds"]
            pows = turbine_dat["performance"]["power_curve"]["power_values"]
            power_curve_type = "power"
        else:
            raise Exception("Bad Power Curve")
        ct = turbine_dat["performance"]["Ct_curve"]["Ct_values"]
        ct_ws = turbine_dat["performance"]["Ct_curve"]["Ct_wind_speeds"]
        speeds = np.arange(np.min([cp_ws, ct_ws]), np.max([cp_ws, ct_ws]) + 1, 1)
        cts_int = np.interp(speeds, ct_ws, ct)
        if power_curve_type == "power":
            powers = np.interp(speeds, cp_ws, pows)
        else:
            cps_int = np.interp(speeds, cp_ws, cp)
            powers = 0.5 * cps_int * speeds**3 * 1.225 * (rd / 2) ** 2 * np.pi

        if "cutin_wind_speed" in turbine_dat["performance"]:
            cutin = turbine_dat["performance"]["cutin_wind_speed"]
        else:
            cutin = 0
        if "cutout_wind_speed" in turbine_dat["performance"]:
            cutout = turbine_dat["performance"]["cutout_wind_speed"]
        else:
            cutout = None
        this_turbine = WindTurbine(
            name=turbine_dat["name"],
            diameter=rd,
            hub_height=hh,
            powerCtFunction=PowerCtTabular(speeds, powers, power_unit="W", ct=cts_int),
            ws_cutin=cutin,
            ws_cutout=cutout,
        )
        this_turbine.powerCtFunction = PowerCtFunctionList(
            key="operating",
            powerCtFunction_lst=[
                PowerCtTabular(
                    ws=[0, 100], power=[0, 0], power_unit="w", ct=[0, 0]
                ),  # 0=No power and ct
                this_turbine.powerCtFunction,
            ],  # 1=Normal operation
            default_value=1,
        )
        turbines.append(this_turbine)

    if len(turbines) == 1:
        turbine = this_turbine
        turbine_types = 0
    else:
        turbine = WindTurbines.from_WindTurbine_lst(turbines)
        turbine_types = farm_dat["layouts"][0]["turbine_types"]

    # farm = 'examples/plant/plant_wind_farm/IEA37_case_study_3_wind_farm.yaml'
    # with open(farm, "r") as stream:
    #    try:
    #        farm_dat = yaml.safe_load(stream)
    #    except yaml.YAMLError as exc:
    #        print(exc)

    # resource = 'examples/plant/plant_energy_resource/IEA37_case_study_4_energy_resource.yaml'
    # with open(resource, "r") as stream:
    #    try:
    #        resource_dat = yaml.safe_load(stream)
    #    except yaml.YAMLError as exc:
    #        print(exc)
    resource_dat = system_dat["site"]["energy_resource"]
    WFXLB = np.min(system_dat["site"]["boundaries"]["polygons"][0]["x"])
    WFXUB = np.max(system_dat["site"]["boundaries"]["polygons"][0]["x"])
    WFYLB = np.min(system_dat["site"]["boundaries"]["polygons"][0]["y"])
    WFYUB = np.max(system_dat["site"]["boundaries"]["polygons"][0]["y"])
    checkk = "flow_field" in system_dat["attributes"]["model_outputs_specification"]
    if checkk:
        checkk = (
            "z_planes"
            in system_dat["attributes"]["model_outputs_specification"]["flow_field"]
        )
    if (
        checkk
        and "xlb"
        in system_dat["attributes"]["model_outputs_specification"]["flow_field"][
            "z_planes"
        ]
    ):
        WFXLB = system_dat["attributes"]["model_outputs_specification"]["flow_field"][
            "z_planes"
        ]["xlb"]
    if (
        checkk
        and "xub"
        in system_dat["attributes"]["model_outputs_specification"]["flow_field"][
            "z_planes"
        ]
    ):
        WFXUB = system_dat["attributes"]["model_outputs_specification"]["flow_field"][
            "z_planes"
        ]["xub"]
    if (
        checkk
        and "ylb"
        in system_dat["attributes"]["model_outputs_specification"]["flow_field"][
            "z_planes"
        ]
    ):
        WFYLB = system_dat["attributes"]["model_outputs_specification"]["flow_field"][
            "z_planes"
        ]["ylb"]
    if (
        checkk
        and "yub"
        in system_dat["attributes"]["model_outputs_specification"]["flow_field"][
            "z_planes"
        ]
    ):
        WFYUB = system_dat["attributes"]["model_outputs_specification"]["flow_field"][
            "z_planes"
        ]["yub"]

    if (
        checkk
        and "dx"
        in system_dat["attributes"]["model_outputs_specification"]["flow_field"][
            "z_planes"
        ]
    ):
        WFDX = system_dat["attributes"]["model_outputs_specification"]["flow_field"][
            "z_planes"
        ]["dx"]
    else:
        WFDX = (WFXUB - WFXLB) / 100

    if (
        checkk
        and "dy"
        in system_dat["attributes"]["model_outputs_specification"]["flow_field"][
            "z_planes"
        ]
    ):
        WFDY = system_dat["attributes"]["model_outputs_specification"]["flow_field"][
            "z_planes"
        ]["dy"]
    else:
        WFDY = (WFYUB - WFYLB) / 100

    # get x and y positions
    if type(farm_dat["layouts"]) == list:
        x = farm_dat["layouts"][0]["coordinates"]["x"]
        y = farm_dat["layouts"][0]["coordinates"]["y"]
    else:
        x = farm_dat["layouts"]["coordinates"]["x"]
        y = farm_dat["layouts"]["coordinates"]["y"]

    ##################
    # construct site
    ##################
    operating = 1
    additional_heights = []
    if "time" in resource_dat["wind_resource"].keys():
        timeseries = True
        wind_resource_timeseries = resource_dat["wind_resource"]["time"]
        times = wind_resource_timeseries
        cases_idx = np.ones(len(times)).astype(bool)
        if (
            "model_outputs_specification" in system_dat["attributes"]
            and "run_configuration"
            in system_dat["attributes"]["model_outputs_specification"]
        ):
            run_config = system_dat["attributes"]["model_outputs_specification"][
                "run_configuration"
            ]
            if "times_run" in run_config and not run_config["times_run"].get(
                "all_occurences", True
            ):
                if "subset" in run_config["times_run"]:
                    cases_idx = run_config["times_run"]["subset"]

        if "height" in resource_dat["wind_resource"].keys():
            heights = resource_dat["wind_resource"]["height"]
        else:
            heights = None

        # Helper to get data and dimensions safely
        def get_resource_data(var_name):
            data_obj = resource_dat["wind_resource"][var_name]
            # Handle windIO schema structure (data + dims)
            vals = np.array(data_obj["data"])
            # Default to ["time"] if dims missing but data exists
            dims = data_obj.get("dims", ["time"])
            return vals, dims

        # Extract Raw Data
        ws_vals, ws_dims = get_resource_data("wind_speed")
        wd_vals, wd_dims = get_resource_data("wind_direction")

        # Apply subsetting (cases_idx)
        # Assuming 'time' is always the first dimension
        ws_vals = ws_vals[cases_idx]
        wd_vals = wd_vals[cases_idx]

        # PREPARE REFERENCE ARRAYS FOR SIMULATION CALL
        # If data is turbine specific (2D: time x turbine), flatten to 1D (time)
        if "wind_turbine" in ws_dims:
            ws = np.mean(ws_vals, axis=1)
        else:
            ws = ws_vals

        if "wind_turbine" in wd_dims:
            # Vector mean for direction to handle 360/0 boundary
            rads = np.deg2rad(wd_vals)
            mean_sin = np.mean(np.sin(rads), axis=1)
            mean_cos = np.mean(np.cos(rads), axis=1)
            wd = np.mod(np.rad2deg(np.arctan2(mean_sin, mean_cos)), 360)
        else:
            wd = wd_vals

        if "operating" in resource_dat["wind_resource"].keys():
            operating = np.array(resource_dat["wind_resource"]["operating"]["data"])[
                cases_idx
            ].T
            assert operating.shape[0] == len(x)
        else:
            operating = np.ones((len(x), len(cases_idx)))

        if len(hub_heights) > 1:
            speeds = []
            dirs = []
            seen = []
            if (
                "z_planes"
                in system_dat["attributes"]["model_outputs_specification"]["flow_field"]
            ):
                z_planes = system_dat["attributes"]["model_outputs_specification"][
                    "flow_field"
                ]
                if z_planes != "hub_heights":
                    additional_heights = system_dat["attributes"][
                        "model_outputs_specification"
                    ]["flow_field"]["z_planes"]["z_list"]
            else:
                additional_heights = []
            for hh in sorted(np.append(list(hub_heights.values()), additional_heights)):
                if hh in seen:
                    continue
                seen.append(hh)
                if heights:
                    try:
                        ws_int = interp1d(
                            heights, ws, axis=1, fill_value="extrapolate"
                        )(hh)
                        wd_int = interp1d(
                            heights, wd, axis=1, fill_value="extrapolate"
                        )(hh)
                    except ValueError:
                        ws_int = interp1d(
                            heights, np.array(ws).T, axis=1, fill_value="extrapolate"
                        )(hh)
                        wd_int = interp1d(
                            heights, np.array(wd).T, axis=1, fill_value="extrapolate"
                        )(hh)
                else:
                    ws_int = ws
                    wd_int = wd
                speeds.append(ws_int)
                dirs.append(wd_int)
            ws = ws_int
            wd = wd_int
        else:
            print(np.array(ws).shape, np.array(heights).shape)
            if heights:
                try:
                    ws = interp1d(heights, ws, axis=1, fill_value="extrapolate")(hh)
                    wd = interp1d(heights, wd, axis=1, fill_value="extrapolate")(hh)
                except ValueError:
                    ws = interp1d(
                        heights, np.array(ws).T, axis=1, fill_value="extrapolate"
                    )(hh)
                    wd = interp1d(
                        heights, np.array(wd).T, axis=1, fill_value="extrapolate"
                    )(hh)
            assert len(np.array(times)[cases_idx]) == len(ws)
            assert len(wd) == len(ws)
            if "wind_turbine" in ws_dims or "wind_turbine" in wd_dims:
                # Turbine-specific data -> Use XRSite via dict_to_site
                site = dict_to_site(resource_dat["wind_resource"])
            else:
                # Uniform data -> Use Hornsrev1Site (Old behavior)
                # This preserves behavior for KUL and 4wts tests which use uniform inflow
                site = Hornsrev1Site()
        if "turbulence_intensity" not in resource_dat["wind_resource"]:
            TI = 0.02
        else:
            TI = np.array(
                resource_dat["wind_resource"]["turbulence_intensity"]["data"]
            )[cases_idx]
            if len(hub_heights) > 1:
                TIs = []
                seen = []
                for hh in sorted(
                    np.append(
                        list(hub_heights.values()),
                        system_dat["attributes"]["model_outputs_specification"][
                            "flow_field"
                        ]["z_planes"]["z_list"],
                    )
                ):
                    # hh = hub_heights[tt]
                    if heights:
                        if hh in seen:
                            continue
                        seen.append(hh)
                        try:
                            ti_int = np.maximum(
                                interp1d(heights, TI, axis=1, fill_value="extrapolate")(
                                    hh
                                ),
                                2e-2,
                            )
                        except ValueError:
                            ti_int = np.maximum(
                                interp1d(
                                    heights,
                                    np.array(TI).T,
                                    axis=1,
                                    fill_value="extrapolate",
                                )(hh),
                                2e-2,
                            )
                    else:
                        ti_int = TI
                    TIs.append(ti_int[cases_idx])
                TI = ti_int
                site = XRSite(
                    xr.Dataset(
                        data_vars={
                            "WS": (["h", "time"], np.array(speeds)),
                            "WD": (["h", "time"], np.array(dirs)),
                            "TI": (["h", "time"], np.array(TIs)),
                            "P": 1,
                        },
                        coords={"h": seen, "time": np.arange(len(times))},
                    )
                )
            else:
                if heights:
                    TI = interp1d(heights, TI, axis=1)(hh)

        # ite = XRSite(xr.Dataset(
        # data_vars={'P': (('time'), np.ones(len(ws)) / len(speeds)), },
        # coords={'time': range(len(times)),
        #         'ws': speeds,
        #         'wd': wd}))

    elif "weibull_k" in resource_dat["wind_resource"].keys():
        A = resource_dat["wind_resource"]["weibull_a"]
        k = resource_dat["wind_resource"]["weibull_k"]
        freq = resource_dat["wind_resource"]["sector_probability"]
        wd = resource_dat["wind_resource"]["wind_direction"]
        if "wind_speed" in resource_dat["wind_resource"]:
            ws = resource_dat["wind_resource"]["wind_speed"]
        else:
            ws = np.arange(2, 30, 1)

        if (
            "wind_turbine"
            in resource_dat["wind_resource"]["sector_probability"]["dims"]
        ):
            mean_ws = np.array(A["data"]) * gamma(
                1 + 1.0 / np.array(k["data"])
            )  # shape (i,wd)
            max_mean = np.max(mean_ws, axis=0)  # shape (wd,)
            Speedup = mean_ws / max_mean  # normalized speed-up (i,wd)
            resource_dat["wind_resource"]["Speedup"] = {
                "dims": ["wind_turbine", "wd"],
                "data": Speedup,
            }

        if all(
            [
                k in resource_dat["wind_resource"]["sector_probability"]["dims"]
                for k in ["x", "y"]
            ]
        ):
            mean_ws = np.array(A["data"]) * gamma(
                1 + 1.0 / np.array(k["data"])
            )  # shape (x,y,h,wd)
            max_mean = np.max(mean_ws, axis=(0, 1))  # shape (h,wd)
            Speedup = mean_ws / max_mean  # normalized speed-up (x,y,h,wd)
            resource_dat["wind_resource"]["Speedup"] = {
                "dims": ["x", "y", "height", "wind_direction"],
                "data": Speedup,
            }

        site = dict_to_site(resource_dat["wind_resource"])

        timeseries = False
        site_ds = dict_to_netcdf(resource_dat["wind_resource"])
        if "x" in site_ds.turbulence_intensity.dims:
            interpolated_ti = site_ds.turbulence_intensity.interp(x=x, y=y)
            if "height" in interpolated_ti.dims:
                interpolated_ti = interpolated_ti.interp(height=hub_heights["0"])
            TI = np.array(
                [interpolated_ti.isel(x=i, y=i).values for i in range(len(x))]
            )
        else:
            TI = resource_dat["wind_resource"]["turbulence_intensity"]["data"]

    else:
        timeseries = False
        ws = resource_dat["wind_resource"]["wind_speed"]
        wd = resource_dat["wind_resource"]["wind_direction"]
        P = np.array(resource_dat["wind_resource"]["probability"]["data"])
        site = dict_to_site(resource_dat["wind_resource"])
        TI = resource_dat["wind_resource"]["turbulence_intensity"]["data"]

    # if 'name' in system_dat['attributes']['analysis']['model_outputs_specification']:
    #   output_dir = system_dat['attributes']['analysis']['model_outputs_specification']['name']
    #   if not os.path.exists(output_dir):
    #      os.makedirs(output_dir)

    wind_deficit_model_data = get_with_default(
        system_dat["attributes"]["analysis"], "wind_deficit_model", DEFAULTS
    )

    deficit_args = {}
    deficit_param_mapping = {}
    wake_deficit_key = None
    print("Running deficit ", wind_deficit_model_data)
    if wind_deficit_model_data["name"] == "Jensen":
        wakeModel = NOJLocalDeficit
        # deficit_param_mapping = {"k": "k", "k2": "k2"}
        if (
            "k_b"
            in system_dat["attributes"]["analysis"]["wind_deficit_model"][
                "wake_expansion_coefficient"
            ]
        ):
            # Handle k2 if present, default to 0.0 if not
            if (
                "k_a"
                in system_dat["attributes"]["analysis"]["wind_deficit_model"][
                    "wake_expansion_coefficient"
                ]
            ):
                k_a = system_dat["attributes"]["analysis"]["wind_deficit_model"][
                    "wake_expansion_coefficient"
                ]["k_a"]
            else:
                k_a = 0

            k_b = system_dat["attributes"]["analysis"]["wind_deficit_model"][
                "wake_expansion_coefficient"
            ]["k_b"]
            deficit_args["a"] = [k_a, k_b]

    elif wind_deficit_model_data["name"].lower() == "bastankhah2014":
        wakeModel = BastankhahGaussianDeficit
        if (
            "k_b"
            in system_dat["attributes"]["analysis"]["wind_deficit_model"][
                "wake_expansion_coefficient"
            ]
        ):
            deficit_args["k"] = system_dat["attributes"]["analysis"][
                "wind_deficit_model"
            ]["wake_expansion_coefficient"]["k_b"]
        elif (
            "k"
            in system_dat["attributes"]["analysis"]["wind_deficit_model"][
                "wake_expansion_coefficient"
            ]
        ):
            deficit_args["k"] = system_dat["attributes"]["analysis"][
                "wind_deficit_model"
            ]["wake_expansion_coefficient"]["k"]
        if "ceps" in system_dat["attributes"]["analysis"]["wind_deficit_model"]:
            deficit_args["ceps"] = system_dat["attributes"]["analysis"][
                "wind_deficit_model"
            ]["ceps"]
        # deficit_param_mapping = {'k': 'k', 'ceps': 'ceps'}
        # from py_wake.deficit_models.utils import ct2a_mom1d
        # deficit_args['ct2a'] = ct2a_mom1d
    elif wind_deficit_model_data["name"] == "SuperGaussian":
        wakeModel = BlondelSuperGaussianDeficit2020
        # deficit_param_mapping = {'k': 'k', 'ceps': 'ceps'}
    elif wind_deficit_model_data["name"] == "TurboPark":
        wakeModel = TurboGaussianDeficit
        # wake_deficit_key = 'WS_jlk'
    elif wind_deficit_model_data["name"].upper() == "FUGA":
        wakeModel = FugaDeficit
        from pyfuga import get_luts

        lut = get_luts(
            folder="luts",  # Path where all files (intermediate and final) are stored
            zeta0=0,  # Stability parameter
            nkz0=8,
            nbeta=32,
            diameter=rd,
            zhub=hh,
            z0=0.00001,
            zi=500,
            zlow=70,
            zhigh=70,
            lut_vars=["UL"],
            nx=2048,
            ny=512,
            n_cpu=1,
        )
        deficit_args["LUT_path"] = (
            r"luts/LUTs_Zeta0=0.00e+00_8_32_D%.1f_zhub%.1f_zi500_z0=0.00001000_z69.2-72.8_UL_nx2048_ny512_dx44.575_dy11.14375.nc"
            % (rd, hh)
        )
        # deficit_args['LUT_path'] = 'luts/LUTs_Zeta0=0.00e+00_8_32_D%i_zhub%i_zi500_z0=0.00001000_z70.0_UL_nx2048_ny512_dx20.0_dy5.0.nc' % (rd, hh)
    else:
        raise Exception(
            "%s wake model not implemented in PyWake" % wind_deficit_model_data["name"]
        )
    for key in wind_deficit_model_data.keys():
        if key == "name":
            continue
        if key in deficit_param_mapping.keys():
            deficit_args[deficit_param_mapping[key]] = wind_deficit_model_data[key]

    if "k2" in deficit_args:
        k = deficit_args.pop("k")
        k2 = deficit_args.pop("k2")
        deficit_args["a"] = [k2, k]

    print("deficit args ", deficit_args)

    # Continuing from the previous example...

    deflection_model_data = get_with_default(
        system_dat["attributes"]["analysis"], "deflection_model", DEFAULTS
    )
    turbulence_model_data = get_with_default(
        system_dat["attributes"]["analysis"], "turbulence_model", DEFAULTS
    )
    superposition_model_data = get_with_default(
        system_dat["attributes"]["analysis"], "superposition_model", DEFAULTS
    )
    rotor_averaging_data = get_with_default(
        system_dat["attributes"]["analysis"], "rotor_averaging", DEFAULTS
    )
    blockage_data = get_with_default(
        system_dat["attributes"]["analysis"], "blockage_model", DEFAULTS
    )

    # Map the deflection model
    if deflection_model_data["name"].lower() == "none":
        deflectionModel = None
    elif deflection_model_data["name"].lower() == "jimenez":
        deflectionModel = JimenezWakeDeflection(
            beta=deflection_model_data["beta"]
        )  # Assuming Jimenez takes 'beta' as an argument
    elif deflection_model_data["name"] == "None":
        deflectionModel = None
    else:
        raise Exception(
            "%s deflection model not implemented" % deflection_model_data["name"]
        )
    deficit_args["use_effective_ws"] = True

    # Map the turbulence model
    if turbulence_model_data["name"].lower() == "none":
        turbulenceModel = None
    elif turbulence_model_data["name"].upper() == "STF2005":
        turbulenceModel = STF2005TurbulenceModel(
            c=[turbulence_model_data["c1"], turbulence_model_data["c2"]]
        )
    elif turbulence_model_data["name"].upper() == "STF2017":
        turbulenceModel = STF2017TurbulenceModel(
            c=[turbulence_model_data["c1"], turbulence_model_data["c2"]]
        )
        # turbulenceModel = STF(turbulence_model_data['c1'], turbulence_model_data['c2'])  # Assuming STF takes 'c1' and 'c2' as arguments
    elif turbulence_model_data["name"].upper() == "CRESPOHERNANDEZ":
        turbulenceModel = CrespoHernandez()
    else:
        raise Exception(
            "%s turbulence model not implemented" % turbulence_model_data["name"]
        )

    # Map the superposition model
    if superposition_model_data["ws_superposition"].lower() == "linear":
        superpositionModel = LinearSum()
    elif superposition_model_data["ws_superposition"].lower() == "squared":
        superpositionModel = SquaredSum()
    else:
        raise Exception(
            "%s superposition model not implemented" % superposition_model_data["name"]
        )

    print("using superposition ", superposition_model_data)
    # Map the rotor averaging model
    if rotor_averaging_data["name"].lower() == "center":
        print("Using Center Average")
        rotorAveraging = RotorCenter()
    elif rotor_averaging_data["name"].lower() == "avg_deficit":
        rotorAveraging = GridRotorAvg()
    else:
        raise Exception(
            "%s rotor averaging model not implemented" % rotor_averaging_data["name"]
        )

    if blockage_data["name"] == "None" or blockage_data["name"] is None:
        blockage = None
    elif blockage_data["name"] == "SelfSimilarityDeficit2020":
        blockage = SelfSimilarityDeficit2020(ss_alpha=blockage_data["ss_alpha"])
    elif blockage_data["name"].upper() == "FUGA":
        blockage = FugaDeficit(deficit_args["LUT_path"])
    else:
        raise Exception("Bad Blockage Specified")

    solver_args = {}
    if blockage is not None:
        solver = All2AllIterative
        solver_args["blockage_deficitModel"] = blockage
    else:
        solver = PropagateDownwind

    print("Running ", wakeModel, deficit_args)
    deficit_model = wakeModel(
        rotorAvgModel=rotorAveraging, groundModel=None, **deficit_args
    )
    if wake_deficit_key:
        deficit_model.WS_key = wake_deficit_key

    windFarmModel = solver(
        site,
        turbine,
        wake_deficitModel=deficit_model,
        superpositionModel=superpositionModel,
        deflectionModel=deflectionModel,
        turbulenceModel=turbulenceModel,
        **solver_args,
    )
    # noj = NOJ(site, turbine, turbulenceModel=None)
    #    sim_res = noj(x, y)
    # sim_res = windFarmModel(x, y, type=turbine_types, time=timeseries, ws=ws, wd=wd, TI=0, yaw=0, tilt=0)
    sim_res = windFarmModel(
        x,
        y,
        type=turbine_types,
        time=timeseries,
        ws=ws,
        wd=wd,
        # TI=TI, # TI is defined in XRSite
        yaw=0,
        tilt=0,
        operating=operating,
    )
    aep = sim_res.aep(normalize_probabilities=not timeseries).sum()
    print("aep is ", aep, "GWh")
    # print('aep is ', sim_res.aep().sum(), 'GWh')
    # print('(%.2f capcacity factor)' % ( aep / (len(x) * turbine.power(10000) * 8760 / 1e9)))

    ######################
    # Construct Outputs
    ####################
    # turbine specific AEP
    if timeseries:
        aep_per_turbine = (
            sim_res.aep(normalize_probabilities=True).sum(["time"]).to_numpy()
        )
    else:
        aep_per_turbine = (
            sim_res.aep(normalize_probabilities=True).sum(["ws", "wd"]).to_numpy()
        )

    # names -- to be updated according to team consensus
    # (also, TODO: read these in from the input WES file!)
    name = "FLOW tool output"
    wind_deficit_model = "Bastankhah’s Gaussian wake model"
    statistical_description = "PDF"  # "time_series"
    statistical_dimensions = ["wind_direction", "wind_velocity"]

    # flow field NetCDF file
    # -----------------
    wind_output_file = "FarmFlow.nc"
    wind_output_variables = [
        "velocity_u",
        "turbulence_intensity",
        "turbulence_k",
        "pressure",
        "wind_loss_to_inflow",
    ]

    power_output_file = "FarmPower.nc"
    power_output_variables = [
        "power_per_turbine",
        "load_per_turbine",
        "Remaining_Useful_life_per_turbine",
        "power_loss_to_inflow_per_turbine",
    ]

    # Create a dictionary to represent the YAML file structure
    # Create an OrderedDict to represent the YAML file structure
    data = OrderedDict(
        [
            ("name", name),
            (
                "FLOW_simulation_config",
                dict(
                    [
                        ("tool", "PyWake"),
                        ("wind_deficit_model", wind_deficit_model),
                        # ('wind_energy_system', '!include recorded_inputs.yaml')
                        ("wind_energy_system", "INCLUDE_YAML_PLACEHOLDER"),
                    ]
                ),
            ),
            (
                "FLOW_simulation_outputs",
                dict(
                    [
                        ("statistical_description", statistical_description),
                        # ('statistical_dimensions', statistical_dimensions),
                        # ('power_percentiles', power_percentiles),
                        # ('AEP', aep),
                        # ('AEP_per_turbine', aep_per_turbine),
                        # ('lifetime_per_turbine', lifetime_per_turbine),
                        # ('wind_output_file', wind_output_file),
                        # ('wind_output_variables', wind_output_variables),
                        # ('power_output_file', power_output_file),
                        # ('power_output_variables', power_output_variables)
                    ]
                ),
            ),
        ]
    )

    # if 'AEP' in system_dat['attributes']['outputs']:
    #   data['FLOW_simulation_outputs']['AEP'] = float(aep)

    # if 'power_percentiles' in system_dat['attributes']['analysis']['outputs']:
    # if system_dat['attributes']['analysis']['outputs']['power_percentiles']['report']:
    #   # compute power percentiles
    #   percentiles = np.array(system_dat['attributes']['analysis']['outputs']['power_percentiles']['percentiles']) / 100
    #   if not timeseries:

    #      # Flatten the power and P arrays
    #      power_values = sim_res.Power.sum('wt').values.flatten()
    #      probabilities = sim_res.P.values.flatten()

    #      # Compute the weighted percentiles
    #      power_percentiles = weighted_quantile(power_values, percentiles, sample_weight=probabilities)

    #   else:
    #       total_power = sim_res.Power.sum(dim='wt')
    #       power_percentiles = total_power.quantile(percentiles, dim=['time']).to_numpy()
    #   data['FLOW_simulation_outputs']['computed_percentiles'] = system_dat['attributes']['analysis']['outputs']['power_percentiles']['percentiles']
    #   data['FLOW_simulation_outputs']['power_percentiles'] = power_percentiles

    os.makedirs(output_dir, exist_ok=True)
    if "turbine_outputs" in system_dat["attributes"]["model_outputs_specification"]:
        # print('aep per turbine', list(aep_per_turbine)); hey
        # data['FLOW_simulation_outputs']['AEP_per_turbine'] = [float(value) for value in aep_per_turbine]
        sim_res_formatted = sim_res[["Power", "WS_eff"]].rename(
            {"Power": "power", "WS_eff": "effective_wind_speed", "wt": "turbine"}
        )
        turbine_nc_filename = str(
            system_dat["attributes"]
            .get("model_outputs_specification", {})
            .get("turbine_outputs", {})
            .get("turbine_nc_filename", "PowerTable.nc")
        )
        turbine_nc_filepath = output_dir + os.sep + turbine_nc_filename
        sim_res_formatted.to_netcdf(turbine_nc_filepath)

    print(sim_res)

    # flow field handling
    flow_map = None
    if (
        "flow_field" in system_dat["attributes"]["model_outputs_specification"]
        and not timeseries
    ):
        # z_planes = system_dat['attributes']['model_outputs_specification']['flow_field']
        # sorted(np.append(list(hub_heights.values()), additional_heights))
        # if 'x_bounds' in  z_planes
        # compute flow map for specified directions (wd) and speeds (ws)
        if timeseries:
            flow_map = sim_res.flow_box(
                x=np.arange(WFXLB, WFXUB + WFDX, WFDX),
                y=np.arange(WFYLB, WFYUB + WFDY, WFDY),
                h=additional_heights,
                time=sim_res.time,
                # operating=operating # TODO
            )

            # remove unwanted data
            flow_map = flow_map.drop_vars(["WD", "WS", "TI", "P"])
        else:
            flow_map = sim_res.flow_box(
                x=np.arange(WFXLB, WFXUB + WFDX, WFDX),
                y=np.arange(WFYLB, WFYUB + WFDY, WFDY),
                h=list(hub_heights.values())
                # operating=operating # TODO
            )

        # raise warning if user requests data we can not provide
        if any(
            element not in ["velocity_u", "turbulence_intensity"]
            for element in system_dat["attributes"]["model_outputs_specification"][
                "flow_field"
            ]["output_variables"]
        ):
            warnings.warn("PyWake can only output velocity_u and turbulence_intensity")

        # remove TI or WS if they are not requested
        if (
            "turbulence_intensity"
            not in system_dat["attributes"]["model_outputs_specification"][
                "flow_field"
            ]["output_variables"]
        ):
            # flow_map = flow_map.drop_vars(["TI_eff"])
            pass
        if (
            "velocity_u"
            not in system_dat["attributes"]["model_outputs_specification"][
                "flow_field"
            ]["output_variables"]
        ):
            pass
            # flow_map = flow_map.drop_vars(["WS_eff"])

    elif (
        "flow_field" in system_dat["attributes"]["model_outputs_specification"]
        and timeseries
    ):
        # time_to_plot = system_dat['attributes']['outputs']['flow_field']['time'][0]
        # print('TIMES ', time_to_plot)
        # print(system_dat['attributes']['model_outputs_specification']['flow_field']['z_planes']['z_list'])
        # print("Flow box bounds: ", WFXLB, WFXUB, WFYLB, WFYUB)
        if (
            system_dat["attributes"]["model_outputs_specification"]["flow_field"][
                "report"
            ]
            is not False
        ):
            z_planes = system_dat["attributes"]["model_outputs_specification"][
                "flow_field"
            ]
            if "z_list" in z_planes:
                additional_heights = z_planes["z_list"]
            else:
                additional_heights = sorted(list(hub_heights.values()))
            flow_map = sim_res.flow_box(
                x=np.arange(WFXLB, WFXUB + WFDX, WFDX),
                y=np.arange(WFYLB, WFYUB + WFDY, WFDY),
                h=additional_heights,
                time=sim_res.time.values,
                # operating=operating,
            )
        # flow_map = sim_res.flow_map(HorizontalGrid(x = np.linspace(WFXLB, WFXUB, 100),
    # y = np.linspace(WFYLB, WFYUB, 100)))
    # wd=wd,
    # ws=ws)
    # power table
    # if 'power_table' in system_dat['attributes']['analysis']['outputs']:
    #   # todo: more in depth stuff here, include loads
    #   #data['FLOW_simulation_outputs']['power_output_variables'] = tuple(system_dat['attributes']['analysis']['outputs']['flow_field']['output_variables'])
    #
    if flow_map:
        # save data
        flow_map = (
            flow_map[["WS_eff", "TI_eff"]]
            # .drop(["wd", "ws"])
            .rename(
                {
                    "h": "z",
                    "WS_eff": "wind_speed",
                    "TI_eff": "turbulence_intensity",
                }
            )
        )
        flow_map.to_netcdf(output_dir + os.sep + "FarmFlow.nc")

        # record data
        data["FLOW_simulation_outputs"]["wind_output_file"] = "FarmFlow.nc"
        data["FLOW_simulation_outputs"]["flow_field"] = system_dat["attributes"][
            "model_outputs_specification"
        ]["flow_field"]

    # Write out the YAML data
    output_yaml_nam = output_dir + os.sep + "output.yaml"
    #    with open(output_yaml_nam, 'w') as file:
    #        yaml.dump(data, file, default_flow_style=False, allow_unicode=True)
    #
    #    with open(output_yaml_nam, 'r') as f:
    #       yaml_content = f.read()
    #
    #    yaml_content = yaml_content.replace('INCLUDE_YAML_PLACEHOLDER', '!include recorded_inputs.yaml')
    #
    #    # Save the post-processed YAML
    #    with open(output_yaml_nam, 'w') as f:
    #       f.write(yaml_content)

    ######################
    # Construct Outputs
    ####################
    # Create a dictionary to represent the simplified YAML file structure
    data = {
        "wind_energy_system": "INCLUDE_YAML_PLACEHOLDER",
        "power_table": "INCLUDE_POWER_TABLE_PLACEHOLDER",
        "flow_field": "INCLUDE_FLOW_FIELD_PLACEHOLDER",
    }

    # Write out the YAML data to the specified output directory
    output_yaml_name = os.path.join(output_dir, "output.yaml")
    with open(output_yaml_name, "w") as file:
        yaml_content = yaml.dump(
            data, file, default_flow_style=False, allow_unicode=True
        )

    # Replace placeholders with actual include directives, ensuring no quotes are used
    with open(output_yaml_name, "r") as file:
        yaml_content = file.read()

    # yaml_content = yaml_content.replace('INCLUDE_YAML_PLACEHOLDER', '!include recorded_inputs.yaml')
    #
    ## Save the post-processed YAML
    # with open(output_yaml_nam, 'w') as f:
    #   f.write(yaml_content)
    data = {
        "wind_energy_system": "INCLUDE_YAML_PLACEHOLDER",
        "power_table": "INCLUDE_POWER_TABLE_PLACEHOLDER",
        "flow_field": "INCLUDE_FLOW_FIELD_PLACEHOLDER",
    }

    # Write out the YAML data to the specified output directory
    output_yaml_name = os.path.join(output_dir, "output.yaml")
    with open(output_yaml_name, "w") as file:
        yaml_content = yaml.dump(
            data, file, default_flow_style=False, allow_unicode=True
        )

    # Replace placeholders with actual include directives, ensuring no quotes are used
    with open(output_yaml_name, "r") as file:
        yaml_content = file.read()

    yaml_content = yaml_content.replace(
        "INCLUDE_YAML_PLACEHOLDER", "!include recorded_inputs.yaml"
    )
    yaml_content = yaml_content.replace(
        "INCLUDE_POWER_TABLE_PLACEHOLDER", "!include PowerTable.nc"
    )
    yaml_content = yaml_content.replace(
        "INCLUDE_FLOW_FIELD_PLACEHOLDER", "!include FarmFlow.nc"
    )

    # Save the post-processed YAML
    with open(output_yaml_name, "w") as file:
        file.write(yaml_content)

    return aep


def run():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_yaml", help="The input yaml file")
    args = parser.parse_args()

    run_pywake(args.input_yaml)


if __name__ == "__main__":
    run()
