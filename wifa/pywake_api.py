import argparse
import warnings
from pathlib import Path

import numpy as np
import xarray as xr
from scipy.interpolate import interp1d
from scipy.special import gamma

from wifa._optional import require

from wifa._optional import require

# Define default values for wind_deficit_model parameters


def _normalize_name(name):
    """Normalize model name for case-insensitive matching."""
    return name.strip().lower().replace("-", "").replace("_", "")


DEFAULTS = {
    "wind_deficit_model": {
        "name": "Jensen",
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
        warnings.warn(f"Using default value for {key}")
        return defaults[key]

    if isinstance(data[key], dict):
        # Merge defaults into user dict: fill missing keys from defaults,
        # but preserve all extra user-provided keys (e.g. n, n_x_grid_points).
        # Recurse when both user and default values are dicts.
        merged = dict(data[key])
        for sub_key in defaults[key]:
            if sub_key not in merged:
                warnings.warn(f"Using default value for {sub_key}")
                merged[sub_key] = defaults[key][sub_key]
            elif isinstance(merged[sub_key], dict) and isinstance(
                defaults[key][sub_key], dict
            ):
                merged[sub_key] = get_with_default(data[key], sub_key, defaults[key])
        return merged

    return data[key]


def load_and_validate_config(yaml_input, default_output_dir="output"):
    """Load and validate a wind energy system YAML configuration.

    Args:
        yaml_input: Path to YAML file (str) or pre-parsed dict
        default_output_dir: Default output directory if not specified in config

    Returns:
        tuple: (system_dat, output_dir) where system_dat is the parsed config dict
    """
    from windIO import load_yaml
    from windIO import validate as validate_yaml

    if not isinstance(yaml_input, dict):
        validate_yaml(yaml_input, "plant/wind_energy_system")
        system_dat = load_yaml(Path(yaml_input))
    else:
        system_dat = yaml_input

    # output_dir priority: 1) yaml file, 2) function argument, 3) default
    output_dir = str(
        system_dat["attributes"]
        .get("model_outputs_specification", {})
        .get("output_folder", default_output_dir)
    )

    Path(output_dir).mkdir(parents=True, exist_ok=True)

    return system_dat, output_dir


def create_turbines(farm_dat):
    """Create turbine objects from farm configuration.

    Args:
        farm_dat: Farm data dictionary containing turbine specifications

    Returns:
        tuple: (turbine, turbine_types, hub_heights) where:
            - turbine: WindTurbine or WindTurbines object
            - turbine_types: int (0) for single turbine or list of types
            - hub_heights: dict mapping type names to hub heights
    """
    from py_wake.wind_turbines import WindTurbine, WindTurbines
    from py_wake.wind_turbines.power_ct_functions import (
        PowerCtFunctionList,
        PowerCtTabular,
    )

    # Handle single vs multiple turbine types
    if "turbines" in farm_dat:
        turbine_dats = [farm_dat["turbines"]]
        type_names = "0"
    else:
        turbine_dats = list(farm_dat["turbine_types"].values())
        type_names = list(farm_dat["turbine_types"].keys())

    turbines = []
    hub_heights = {}

    for turbine_dat, key in zip(turbine_dats, type_names):
        hh = turbine_dat["hub_height"]
        rd = turbine_dat["rotor_diameter"]
        hub_heights[key] = hh

        # Parse power/Cp curves
        if "Cp_curve" in turbine_dat["performance"]:
            cp = turbine_dat["performance"]["Cp_curve"]["Cp_values"]
            cp_ws = turbine_dat["performance"]["Cp_curve"]["Cp_wind_speeds"]
            power_curve_type = "cp"
        elif "power_curve" in turbine_dat["performance"]:
            cp_ws = turbine_dat["performance"]["power_curve"]["power_wind_speeds"]
            pows = turbine_dat["performance"]["power_curve"]["power_values"]
            power_curve_type = "power"
        else:
            raise ValueError(
                "Missing Cp_curve or power_curve in turbine performance data"
            )

        ct = turbine_dat["performance"]["Ct_curve"]["Ct_values"]
        ct_ws = turbine_dat["performance"]["Ct_curve"]["Ct_wind_speeds"]
        speeds = np.arange(np.min([cp_ws, ct_ws]), np.max([cp_ws, ct_ws]) + 1, 1)
        cts_int = np.interp(speeds, ct_ws, ct)

        if power_curve_type == "power":
            powers = np.interp(speeds, cp_ws, pows)
        else:
            cps_int = np.interp(speeds, cp_ws, cp)
            powers = 0.5 * cps_int * speeds**3 * 1.225 * (rd / 2) ** 2 * np.pi

        cutin = turbine_dat["performance"].get("cutin_wind_speed", 0)
        cutout = turbine_dat["performance"].get("cutout_wind_speed")

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
        turbine = turbines[0]
        turbine_types = 0
    else:
        turbine = WindTurbines.from_WindTurbine_lst(turbines)
        turbine_types = farm_dat["layouts"][0]["turbine_types"]

    return turbine, turbine_types, hub_heights


def dict_to_site(resource_dict):
    """Convert a wind resource dictionary to a PyWake XRSite.

    Args:
        resource_dict: Wind resource dictionary from windIO

    Returns:
        XRSite object configured with the wind resource data
    """
    from py_wake.site import XRSite
    from windIO import dict_to_netcdf

    resource_ds = dict_to_netcdf(resource_dict)
    rename_map = {
        "height": "h",
        "weibull_a": "Weibull_A",
        "weibull_k": "Weibull_k",
        "sector_probability": "Sector_frequency",
        "turbulence_intensity": "TI",
        "wind_turbine": "i",
        "density": "Air_density",
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
        # as the leading dimension.
        resource_ds = resource_ds.transpose("i", *other_dims)
    return XRSite(resource_ds)


def get_flow_field_param(system_dat, param_name, default=None):
    """Extract flow field parameter with safe nested access.

    Args:
        system_dat: System data dictionary
        param_name: Name of parameter to extract (e.g., 'xlb', 'dx')
        default: Default value if parameter not found

    Returns:
        Parameter value or default
    """
    try:
        return system_dat["attributes"]["model_outputs_specification"]["flow_field"][
            "z_planes"
        ][param_name]
    except KeyError:
        return default


def construct_site(system_dat, resource_dat, hub_heights, x_positions):
    """Construct site object and wind conditions for simulation.

    Args:
        system_dat: System data dictionary
        resource_dat: Energy resource dictionary
        hub_heights: dict mapping turbine type names to hub heights
        x_positions: list of turbine x positions (for operating array sizing)

    Returns:
        dict with keys: site, ws, wd, TI, timeseries, operating, additional_heights,
                       cases_idx, flow_bounds
    """
    # Compute flow field bounds from site boundaries, with optional overrides
    flow_bounds = _compute_flow_bounds(system_dat)

    # Determine site type and construct accordingly
    wind_resource = resource_dat["wind_resource"]
    if "time" in wind_resource:
        result = _construct_timeseries_site(
            system_dat, resource_dat, hub_heights, x_positions
        )
    elif "weibull_k" in wind_resource:
        result = _construct_weibull_site(resource_dat, hub_heights, x_positions)
    else:
        result = {
            "site": dict_to_site(wind_resource),
            "ws": wind_resource["wind_speed"],
            "wd": wind_resource["wind_direction"],
            "TI": wind_resource["turbulence_intensity"]["data"],
            "timeseries": False,
            "operating": np.ones((len(x_positions), 1)),
            "additional_heights": [],
            "cases_idx": np.ones(1).astype(bool),
        }

    result["flow_bounds"] = flow_bounds
    return result


def _compute_flow_bounds(system_dat):
    """Compute flow field bounds from site boundaries with optional overrides."""
    boundaries = system_dat["site"]["boundaries"]["polygons"][0]
    xlb = get_flow_field_param(system_dat, "xlb", np.min(boundaries["x"]))
    xub = get_flow_field_param(system_dat, "xub", np.max(boundaries["x"]))
    ylb = get_flow_field_param(system_dat, "ylb", np.min(boundaries["y"]))
    yub = get_flow_field_param(system_dat, "yub", np.max(boundaries["y"]))
    return {
        "xlb": xlb,
        "xub": xub,
        "ylb": ylb,
        "yub": yub,
        "dx": get_flow_field_param(system_dat, "dx", (xub - xlb) / 100),
        "dy": get_flow_field_param(system_dat, "dy", (yub - ylb) / 100),
    }


def _construct_timeseries_site(system_dat, resource_dat, hub_heights, x_positions):
    """Construct site from timeseries data.

    Internal helper for construct_site().
    """
    from py_wake.examples.data.hornsrev1 import Hornsrev1Site
    from py_wake.site import XRSite

    wind_resource = resource_dat["wind_resource"]
    times = wind_resource["time"]
    cases_idx = np.ones(len(times)).astype(bool)

    # Check for subset configuration
    times_run = (
        system_dat["attributes"]
        .get("model_outputs_specification", {})
        .get("run_configuration", {})
        .get("times_run", {})
    )
    if not times_run.get("all_occurences", True) and "subset" in times_run:
        cases_idx = times_run["subset"]

    heights = wind_resource.get("height")

    # Helper to get data and dimensions safely
    def get_resource_data(var_name):
        data_obj = wind_resource[var_name]
        vals = np.array(data_obj["data"])
        dims = data_obj.get("dims", ["time"])
        return vals, dims

    # Extract raw data
    ws_vals, ws_dims = get_resource_data("wind_speed")
    wd_vals, wd_dims = get_resource_data("wind_direction")

    # Apply subsetting
    ws_vals = ws_vals[cases_idx]
    wd_vals = wd_vals[cases_idx]

    # Prepare reference arrays - average across turbines if turbine-specific
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

    # Handle operating status
    if "operating" in wind_resource:
        operating = np.array(wind_resource["operating"]["data"])[cases_idx].T
        assert operating.shape[0] == len(x_positions)
    else:
        operating = np.ones((len(x_positions), len(cases_idx)))

    # Handle multi-height interpolation
    additional_heights = []
    hh = list(hub_heights.values())[0]
    site = None

    if len(hub_heights) > 1:
        # Multiple turbine types - need height interpolation
        flow_field_spec = (
            system_dat["attributes"]
            .get("model_outputs_specification", {})
            .get("flow_field", {})
        )
        if (
            "z_planes" in flow_field_spec
            and flow_field_spec["z_planes"] != "hub_heights"
        ):
            additional_heights = flow_field_spec.get("z_planes", {}).get("z_list", [])

        speeds, dirs, TIs, seen = [], [], [], []
        for hh in sorted(np.append(list(hub_heights.values()), additional_heights)):
            if hh in seen:
                continue
            seen.append(hh)
            ws_int, wd_int = _interpolate_wind_data(heights, ws, wd, hh)
            speeds.append(ws_int)
            dirs.append(wd_int)

        ws, wd = ws_int, wd_int

        # Handle TI interpolation
        if "turbulence_intensity" not in wind_resource:
            TI = 0.02
        else:
            TI_data = np.array(wind_resource["turbulence_intensity"]["data"])[cases_idx]
            for hh in sorted(np.append(list(hub_heights.values()), additional_heights)):
                if hh in seen[len(speeds) :]:
                    continue
                if heights:
                    ti_int = _interpolate_with_min(heights, TI_data, hh, min_val=0.02)
                else:
                    ti_int = TI_data
                TIs.append(ti_int)
            TI = ti_int

            data_vars = {
                "WS": (["h", "time"], np.array(speeds)),
                "WD": (["h", "time"], np.array(dirs)),
                "TI": (["h", "time"], np.array(TIs)),
                "P": 1,
            }
            if "density" in wind_resource:
                density_vals = np.array(wind_resource["density"]["data"])[cases_idx]
                density_dims = wind_resource["density"].get("dims", ["time"])
                if "wind_turbine" in density_dims:
                    density_vals = np.mean(density_vals, axis=1)
                data_vars["Air_density"] = (["time"], density_vals)
            site = XRSite(
                xr.Dataset(
                    data_vars=data_vars,
                    coords={"h": seen, "time": np.arange(len(times))},
                )
            )
    else:
        # Single turbine type
        if heights:
            ws, wd = _interpolate_wind_data(heights, ws, wd, hh)

        assert len(np.array(times)[cases_idx]) == len(ws)
        assert len(wd) == len(ws)

        if "wind_turbine" in ws_dims or "wind_turbine" in wd_dims:
            site = dict_to_site(wind_resource)
        else:
            site = Hornsrev1Site()
            if "density" in wind_resource:
                density_vals = np.array(wind_resource["density"]["data"])[cases_idx]
                site.ds["Air_density"] = (("time",), density_vals)

        # Handle TI
        if "turbulence_intensity" not in wind_resource:
            TI = 0.02
        else:
            TI = np.array(wind_resource["turbulence_intensity"]["data"])[cases_idx]
            if heights:
                TI = interp1d(heights, TI, axis=1)(hh)

    return {
        "site": site,
        "ws": ws,
        "wd": wd,
        "TI": TI,
        "timeseries": True,
        "operating": operating,
        "additional_heights": additional_heights,
        "cases_idx": cases_idx,
    }


def _construct_weibull_site(resource_dat, hub_heights, x_positions):
    """Construct site from Weibull distribution data.

    Internal helper for construct_site().
    """
    from windIO import dict_to_netcdf

    wind_resource = resource_dat["wind_resource"]
    A = wind_resource["weibull_a"]
    k = wind_resource["weibull_k"]
    wd = wind_resource["wind_direction"]
    ws = wind_resource.get("wind_speed", np.arange(2, 30, 1))

    # Handle turbine-specific Weibull
    if "wind_turbine" in wind_resource["sector_probability"]["dims"]:
        mean_ws = np.array(A["data"]) * gamma(1 + 1.0 / np.array(k["data"]))
        max_mean = np.max(mean_ws, axis=0)
        Speedup = mean_ws / max_mean
        wind_resource["Speedup"] = {
            "dims": ["wind_turbine", "wd"],
            "data": Speedup,
        }

    # Handle spatial Weibull
    if all(key in wind_resource["sector_probability"]["dims"] for key in ["x", "y"]):
        mean_ws = np.array(A["data"]) * gamma(1 + 1.0 / np.array(k["data"]))
        max_mean = np.max(mean_ws, axis=(0, 1))
        Speedup = mean_ws / max_mean
        wind_resource["Speedup"] = {
            "dims": ["x", "y", "height", "wind_direction"],
            "data": Speedup,
        }

    site = dict_to_site(wind_resource)

    # Handle TI
    site_ds = dict_to_netcdf(wind_resource)
    if "x" in site_ds.turbulence_intensity.dims:
        interpolated_ti = site_ds.turbulence_intensity.interp(
            x=x_positions, y=x_positions
        )
        if "height" in interpolated_ti.dims:
            interpolated_ti = interpolated_ti.interp(height=hub_heights["0"])
        TI = np.array(
            [interpolated_ti.isel(x=i, y=i).values for i in range(len(x_positions))]
        )
    else:
        TI = wind_resource["turbulence_intensity"]["data"]

    return {
        "site": site,
        "ws": ws,
        "wd": wd,
        "TI": TI,
        "timeseries": False,
        "operating": np.ones((len(x_positions), 1)),
        "additional_heights": [],
        "cases_idx": np.ones(1).astype(bool),
    }


def _interpolate_wind_data(heights, ws, wd, target_height):
    """Interpolate wind speed and direction to target height.

    Handles automatic transpose for shape mismatches.
    """
    if heights is None:
        return ws, wd

    try:
        ws_int = interp1d(heights, ws, axis=1, fill_value="extrapolate")(target_height)
        wd_int = interp1d(heights, wd, axis=1, fill_value="extrapolate")(target_height)
    except ValueError:
        ws_int = interp1d(heights, np.array(ws).T, axis=1, fill_value="extrapolate")(
            target_height
        )
        wd_int = interp1d(heights, np.array(wd).T, axis=1, fill_value="extrapolate")(
            target_height
        )

    return ws_int, wd_int


def _interpolate_with_min(heights, values, target_height, min_val=0.02):
    """Interpolate values to target height with minimum value clipping."""
    try:
        return np.maximum(
            interp1d(heights, values, axis=1, fill_value="extrapolate")(target_height),
            min_val,
        )
    except ValueError:
        return np.maximum(
            interp1d(heights, np.array(values).T, axis=1, fill_value="extrapolate")(
                target_height
            ),
            min_val,
        )


def configure_wake_model(system_dat, rotor_diameter, hub_height):
    """Configure the wake model components based on system configuration.

    Args:
        system_dat: System data dictionary
        rotor_diameter: Rotor diameter for FUGA LUT generation
        hub_height: Hub height for FUGA LUT generation

    Returns:
        dict with keys: wake_model_class, deficit_args, deflection_model,
                       turbulence_model, superposition_model, rotor_averaging,
                       blockage_model, solver_class, solver_args
    """
    from py_wake.wind_farm_models import All2AllIterative, PropagateDownwind

    analysis = system_dat["attributes"]["analysis"]

    # Resolve each submodel config, filling missing keys from DEFAULTS
    wind_deficit_data = get_with_default(analysis, "wind_deficit_model", DEFAULTS)
    deflection_data = get_with_default(analysis, "deflection_model", DEFAULTS)
    turbulence_data = get_with_default(analysis, "turbulence_model", DEFAULTS)
    superposition_data = get_with_default(analysis, "superposition_model", DEFAULTS)
    rotor_avg_data = get_with_default(analysis, "rotor_averaging", DEFAULTS)
    blockage_data = get_with_default(analysis, "blockage_model", DEFAULTS)

    wake_model_class, deficit_args = _configure_deficit_model(
        wind_deficit_data, analysis, rotor_diameter, hub_height
    )
    deflection_model = _configure_deflection_model(deflection_data)
    turbulence_model = _configure_turbulence_model(turbulence_data)
    superposition_model = _configure_superposition_model(superposition_data)
    rotor_averaging = _configure_rotor_averaging(rotor_avg_data)
    blockage_model = _configure_blockage_model(blockage_data, deficit_args)

    # Blockage requires All2AllIterative solver
    solver_args = {}
    if blockage_model is not None:
        solver_class = All2AllIterative
        solver_args["blockage_deficitModel"] = blockage_model
    else:
        solver_class = PropagateDownwind

    return {
        "wake_model_class": wake_model_class,
        "deficit_args": deficit_args,
        "wake_deficit_key": None,  # Deprecated: kept for API compatibility
        "deflection_model": deflection_model,
        "turbulence_model": turbulence_model,
        "superposition_model": superposition_model,
        "rotor_averaging": rotor_averaging,
        "blockage_model": blockage_model,
        "solver_class": solver_class,
        "solver_args": solver_args,
    }


def _configure_deficit_model(wind_deficit_data, analysis, rotor_diameter, hub_height):
    """Configure the wind deficit model.

    Returns:
        tuple: (wake_model_class, deficit_args)
    """
    from py_wake.deficit_models.fuga import FugaDeficit
    from py_wake.deficit_models.gaussian import (
        BastankhahGaussianDeficit,
        BlondelSuperGaussianDeficit2020,
        BlondelSuperGaussianDeficit2023,
        CarbajofuertesGaussianDeficit,
        NiayifarGaussianDeficit,
        TurboGaussianDeficit,
        ZongGaussianDeficit,
    )
    from py_wake.deficit_models.gcl import GCLDeficit
    from py_wake.deficit_models.noj import NOJLocalDeficit, TurboNOJDeficit

    model_name = wind_deficit_data["name"]
    normalized = _normalize_name(model_name)
    deficit_args = {"use_effective_ws": True}

    wind_deficit_cfg = analysis.get("wind_deficit_model", {})
    wake_expansion = wind_deficit_cfg.get("wake_expansion_coefficient", {})

    GAUSSIAN_MODELS = {
        "bastankhah2014": BastankhahGaussianDeficit,
        "niayifar2016": NiayifarGaussianDeficit,
        "zong2020": ZongGaussianDeficit,
        "carbajofuertes2018": CarbajofuertesGaussianDeficit,
    }
    # Models that accept a=[k_a, k_b] instead of k (scalar)
    A_PARAM_MODELS = {"niayifar2016", "zong2020", "carbajofuertes2018"}

    if normalized in ("jensen", "nojlocaldeficit"):
        wake_model_class = NOJLocalDeficit
        if "k_b" in wake_expansion:
            deficit_args["a"] = [wake_expansion.get("k_a", 0), wake_expansion["k_b"]]

    elif normalized in GAUSSIAN_MODELS:
        wake_model_class = GAUSSIAN_MODELS[normalized]
        if normalized in A_PARAM_MODELS:
            # Niayifar, Zong, Carbajofuertes use a=[k_a, k_b]
            if "k" in wake_expansion:
                warnings.warn(
                    f"{model_name} uses a=[k_a, k_b] for wake expansion, not scalar k. "
                    f"Scalar 'k' is ignored; specify k_a/k_b instead."
                )
            if "k_b" in wake_expansion:
                if "k_a" not in wake_expansion:
                    warnings.warn(
                        f"k_a not specified for {model_name}, defaulting to 0"
                    )
                deficit_args["a"] = [
                    wake_expansion.get("k_a", 0),
                    wake_expansion["k_b"],
                ]
        else:
            # Bastankhah2014 uses k (scalar)
            if "k_b" in wake_expansion:
                deficit_args["k"] = wake_expansion["k_b"]
            elif "k" in wake_expansion:
                deficit_args["k"] = wake_expansion["k"]
        # ceps: valid for Bastankhah, Niayifar, Carbajofuertes (not Zong)
        if normalized != "zong2020" and "ceps" in wind_deficit_cfg:
            deficit_args["ceps"] = wind_deficit_cfg["ceps"]
        # use_effective_ti: valid for Niayifar, Zong, Carbajofuertes (not Bastankhah)
        if normalized != "bastankhah2014" and "use_effective_ti" in wind_deficit_cfg:
            deficit_args["use_effective_ti"] = wind_deficit_cfg["use_effective_ti"]

    elif normalized == "supergaussian":
        wake_model_class = BlondelSuperGaussianDeficit2020

    elif normalized == "supergaussian2023":
        wake_model_class = BlondelSuperGaussianDeficit2023

    elif normalized == "turbopark":
        wake_model_class = TurboGaussianDeficit

    elif normalized == "turbonoj":
        wake_model_class = TurboNOJDeficit
        if "A" in wind_deficit_cfg:
            deficit_args["A"] = wind_deficit_cfg["A"]

    elif normalized == "gcl":
        wake_model_class = GCLDeficit

    elif normalized == "bastankhah2016":
        raise NotImplementedError(
            "Bastankhah2016 is not available in PyWake. Use flow_model 'foxes', "
            "or choose Bastankhah2014/Zong2020 for PyWake."
        )

    elif normalized == "fuga":
        wake_model_class = FugaDeficit
        from pyfuga import get_luts

        get_luts(
            folder="luts",
            zeta0=0,
            nkz0=8,
            nbeta=32,
            diameter=rotor_diameter,
            zhub=hub_height,
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
            f"luts/LUTs_Zeta0=0.00e+00_8_32_D{rotor_diameter:.1f}_zhub{hub_height:.1f}"
            f"_zi500_z0=0.00001000_z69.2-72.8_UL_nx2048_ny512_dx44.575_dy11.14375.nc"
        )

    else:
        raise NotImplementedError(f"Wake model '{model_name}' is not supported")

    return wake_model_class, deficit_args


def _configure_deflection_model(deflection_data):
    """Configure the wake deflection model."""
    from py_wake.deflection_models import JimenezWakeDeflection
    from py_wake.deflection_models.gcl_hill_vortex import GCLHillDeflection

    name = deflection_data.get("name")
    if name is None:
        return None

    normalized = _normalize_name(name)
    if normalized == "none":
        return None
    if normalized == "jimenez":
        return JimenezWakeDeflection(beta=deflection_data["beta"])
    if normalized == "gclhill":
        return GCLHillDeflection()
    if normalized == "bastankhah2016":
        raise NotImplementedError(
            "Bastankhah2016 deflection is not available in PyWake. Use flow_model "
            "'foxes', or choose Jimenez/GCLHill for PyWake."
        )
    raise NotImplementedError(f"Deflection model '{name}' is not supported")


def _configure_turbulence_model(turbulence_data):
    """Configure the turbulence model."""
    from py_wake.turbulence_models import (
        CrespoHernandez,
        STF2005TurbulenceModel,
        STF2017TurbulenceModel,
    )
    from py_wake.turbulence_models.gcl_turb import GCLTurbulence

    name = turbulence_data.get("name")
    if name is None:
        return None

    normalized = _normalize_name(name)
    if normalized == "none":
        return None

    STF_MODELS = {
        "stf2005": STF2005TurbulenceModel,
        "stf2017": STF2017TurbulenceModel,
        "iecti2019": STF2017TurbulenceModel,
    }

    if normalized in STF_MODELS:
        c = [turbulence_data.get("c1", 1.0), turbulence_data.get("c2", 1.0)]
        return STF_MODELS[normalized](c=c)
    if normalized == "crespohernandez":
        return CrespoHernandez()
    if normalized == "gcl":
        return GCLTurbulence()
    raise NotImplementedError(f"Turbulence model '{name}' is not supported")


def _configure_superposition_model(superposition_data):
    """Configure the superposition model."""
    from py_wake.superposition_models import (
        CumulativeWakeSum,
        LinearSum,
        MaxSum,
        SquaredSum,
        WeightedSum,
    )

    name = superposition_data["ws_superposition"]
    normalized = _normalize_name(name)

    SUPERPOSITION_MODELS = {
        "linear": LinearSum,
        "squared": SquaredSum,
        "max": MaxSum,
        "weighted": WeightedSum,
        "cumulative": CumulativeWakeSum,
    }

    if normalized in SUPERPOSITION_MODELS:
        return SUPERPOSITION_MODELS[normalized]()
    if normalized == "product":
        raise NotImplementedError("Product superposition is not available in PyWake.")
    raise NotImplementedError(f"Superposition model '{name}' is not supported")


def _configure_rotor_averaging(rotor_avg_data):
    """Configure the rotor averaging model."""
    from py_wake.rotor_avg_models import (
        CGIRotorAvg,
        EqGridRotorAvg,
        GQGridRotorAvg,
        GridRotorAvg,
        PolarGridRotorAvg,
        RotorCenter,
    )

    name = rotor_avg_data["name"]
    normalized = _normalize_name(name)

    if normalized == "center":
        return RotorCenter()
    if normalized == "avgdeficit":
        return GridRotorAvg()
    if normalized == "eqgrid":
        return EqGridRotorAvg(n=rotor_avg_data.get("n", 4))
    if normalized == "gqgrid":
        return GQGridRotorAvg(
            n_x=rotor_avg_data.get("n_x_grid_points", 4),
            n_y=rotor_avg_data.get("n_y_grid_points", 4),
        )
    if normalized == "polargrid":
        return PolarGridRotorAvg()
    if normalized == "cgi":
        return CGIRotorAvg(n=rotor_avg_data.get("n", 4))
    raise NotImplementedError(f"Rotor averaging model '{name}' is not supported")


def _configure_blockage_model(blockage_data, deficit_args):
    """Configure the blockage model."""
    from py_wake.deficit_models import (
        HybridInduction,
        RankineHalfBody,
        SelfSimilarityDeficit,
        SelfSimilarityDeficit2020,
        VortexCylinder,
        VortexDipole,
    )
    from py_wake.deficit_models.fuga import FugaDeficit
    from py_wake.deficit_models.rathmann import Rathmann

    name = blockage_data["name"]
    if name is None:
        return None

    normalized = _normalize_name(name)
    if normalized == "none":
        return None

    # Models that take no constructor arguments
    SIMPLE_BLOCKAGE_MODELS = {
        "selfsimilaritydeficit": SelfSimilarityDeficit,
        "rankinehalfbody": RankineHalfBody,
        "rathmann": Rathmann,
        "vortexcylinder": VortexCylinder,
        "vortexdipole": VortexDipole,
        "hybridinduction": HybridInduction,
    }

    if normalized == "selfsimilaritydeficit2020":
        return SelfSimilarityDeficit2020(
            ss_alpha=blockage_data.get("ss_alpha", 0.8888888888888888)
        )
    if normalized in SIMPLE_BLOCKAGE_MODELS:
        return SIMPLE_BLOCKAGE_MODELS[normalized]()
    if normalized == "fuga":
        return FugaDeficit(deficit_args["LUT_path"])
    raise NotImplementedError(f"Blockage model '{name}' is not supported")


def run_simulation(site, turbine, wake_config, site_data, x, y, turbine_types):
    """Run the PyWake simulation.

    Args:
        site: Site object (XRSite or similar)
        turbine: WindTurbine or WindTurbines object
        wake_config: dict from configure_wake_model()
        site_data: dict from construct_site()
        x: Turbine x positions
        y: Turbine y positions
        turbine_types: int (0) for single type or list of types

    Returns:
        dict with keys: sim_res, aep, aep_per_turbine
    """
    # Build deficit model
    deficit_model = wake_config["wake_model_class"](
        rotorAvgModel=wake_config["rotor_averaging"],
        groundModel=None,
        **wake_config["deficit_args"],
    )

    # Build wind farm model
    wind_farm_model = wake_config["solver_class"](
        site,
        turbine,
        wake_deficitModel=deficit_model,
        superpositionModel=wake_config["superposition_model"],
        deflectionModel=wake_config["deflection_model"],
        turbulenceModel=wake_config["turbulence_model"],
        **wake_config["solver_args"],
    )

    # Prepare simulation kwargs
    sim_kwargs = {
        "x": x,
        "y": y,
        "type": turbine_types,
        "time": site_data["timeseries"],
        "ws": site_data["ws"],
        "wd": site_data["wd"],
        "yaw": 0,
        "tilt": 0,
        "operating": site_data["operating"],
    }

    # Pass TI if not in site's data variables
    if "TI" not in site.ds.data_vars:
        sim_kwargs["TI"] = site_data["TI"]

    # Run simulation
    sim_res = wind_farm_model(**sim_kwargs)
    is_timeseries = site_data["timeseries"]
    aep = sim_res.aep(normalize_probabilities=not is_timeseries).sum()
    sum_dims = ["time"] if is_timeseries else ["ws", "wd"]
    aep_per_turbine = sim_res.aep(normalize_probabilities=True).sum(sum_dims).to_numpy()

    return {"sim_res": sim_res, "aep": aep, "aep_per_turbine": aep_per_turbine}


def generate_outputs(sim_results, system_dat, site_data, hub_heights, output_dir):
    """Generate output files from simulation results.

    Args:
        sim_results: dict from run_simulation()
        system_dat: System data dictionary
        site_data: dict from construct_site()
        hub_heights: dict mapping type names to hub heights
        output_dir: Output directory path

    Returns:
        float: AEP value
    """
    sim_res = sim_results["sim_res"]
    flow_bounds = site_data["flow_bounds"]
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Write turbine outputs if requested
    output_spec = system_dat["attributes"].get("model_outputs_specification", {})
    if "turbine_outputs" in output_spec:
        sim_res_formatted = sim_res[["Power", "WS_eff"]].rename(
            {"Power": "power", "WS_eff": "effective_wind_speed", "wt": "turbine"}
        )
        turbine_nc_filename = output_spec["turbine_outputs"].get(
            "turbine_nc_filename", "PowerTable.nc"
        )
        sim_res_formatted.to_netcdf(output_path / turbine_nc_filename)

    # Flow field handling
    flow_map = _generate_flow_field(
        sim_res, system_dat, site_data, hub_heights, flow_bounds
    )

    if flow_map is not None:
        flow_map = flow_map[["WS_eff", "TI_eff"]].rename(
            {
                "h": "z",
                "WS_eff": "wind_speed",
                "TI_eff": "turbulence_intensity",
            }
        )
        flow_map.to_netcdf(output_path / "FarmFlow.nc")

    # Write YAML output
    _write_yaml_output(output_dir)

    return sim_results["aep"]


def _generate_flow_field(sim_res, system_dat, site_data, hub_heights, flow_bounds):
    """Generate flow field data if requested.

    Returns:
        Flow map xarray or None
    """
    output_spec = system_dat["attributes"].get("model_outputs_specification", {})
    if "flow_field" not in output_spec:
        return None

    x_range = np.arange(
        flow_bounds["xlb"], flow_bounds["xub"] + flow_bounds["dx"], flow_bounds["dx"]
    )
    y_range = np.arange(
        flow_bounds["ylb"], flow_bounds["yub"] + flow_bounds["dy"], flow_bounds["dy"]
    )

    if not site_data["timeseries"]:
        flow_map = sim_res.flow_box(
            x=x_range,
            y=y_range,
            h=list(hub_heights.values()),
        )
        # Warn if user requests unsupported outputs
        requested_vars = output_spec["flow_field"].get("output_variables", [])
        unsupported = {"velocity_u", "turbulence_intensity"}
        if any(var not in unsupported for var in requested_vars):
            warnings.warn("PyWake can only output velocity_u and turbulence_intensity")
        return flow_map

    # Timeseries flow field
    flow_field_spec = output_spec["flow_field"]
    if flow_field_spec.get("report") is False:
        return None

    z_list = flow_field_spec.get("z_list", sorted(hub_heights.values()))
    return sim_res.flow_box(
        x=x_range,
        y=y_range,
        h=z_list,
        time=sim_res.time.values,
    )


def _write_yaml_output(output_dir):
    """Write the output YAML file with include directives."""
    # Write directly with !include tags (avoids round-trip through yaml.dump)
    content = (
        "flow_field: !include FarmFlow.nc\n"
        "power_table: !include PowerTable.nc\n"
        "wind_energy_system: !include recorded_inputs.yaml\n"
    )
    (Path(output_dir) / "output.yaml").write_text(content)


def run_pywake(yaml_input, output_dir="output"):
    """Run a PyWake wind farm simulation.

    This is the main entry point that orchestrates the simulation workflow:
    1. Load and validate configuration
    2. Create turbine objects
    3. Construct site with wind resource data
    4. Configure wake models
    5. Run simulation
    6. Generate outputs

    Args:
        yaml_input: Path to YAML file (str) or pre-parsed dict
        output_dir: Output directory (can be overridden in YAML config)

    Returns:
        float: Total AEP in GWh
    """
    # Step 1: Load and validate configuration
    require("py_wake")

    system_dat, output_dir = load_and_validate_config(yaml_input, output_dir)

    # Step 2: Create turbine objects
    farm_dat = system_dat["wind_farm"]
    turbine, turbine_types, hub_heights = create_turbines(farm_dat)

    # Get turbine positions
    if isinstance(farm_dat["layouts"], list):
        x = farm_dat["layouts"][0]["coordinates"]["x"]
        y = farm_dat["layouts"][0]["coordinates"]["y"]
    else:
        x = farm_dat["layouts"]["coordinates"]["x"]
        y = farm_dat["layouts"]["coordinates"]["y"]

    # Step 3: Construct site
    resource_dat = system_dat["site"]["energy_resource"]
    site_data = construct_site(system_dat, resource_dat, hub_heights, x)
    site = site_data["site"]

    # Step 4: Configure wake model
    # Use first turbine's dimensions for FUGA LUT if needed
    first_hh = list(hub_heights.values())[0]
    # Get rotor diameter from farm data
    if "turbines" in farm_dat:
        rd = farm_dat["turbines"]["rotor_diameter"]
    else:
        first_key = list(farm_dat["turbine_types"].keys())[0]
        rd = farm_dat["turbine_types"][first_key]["rotor_diameter"]

    wake_config = configure_wake_model(system_dat, rd, first_hh)

    # Step 5: Run simulation
    sim_results = run_simulation(
        site, turbine, wake_config, site_data, x, y, turbine_types
    )

    # Step 6: Generate outputs
    aep = generate_outputs(sim_results, system_dat, site_data, hub_heights, output_dir)

    return aep


def run():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_yaml", help="The input yaml file")
    args = parser.parse_args()

    run_pywake(args.input_yaml)


if __name__ == "__main__":
    run()
