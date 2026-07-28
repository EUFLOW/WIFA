import argparse
import os
import warnings
from pathlib import Path

import numpy as np
import xarray as xr
from scipy.interpolate import interp1d
from scipy.special import gamma

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
    "blockage_model": {
        "name": None,
        "ss_alpha": 0.8888888888888888,
        "ground_mirror": False,
    },
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
        # Keep an included wind_resource.nc as numpy arrays instead of exploding
        # it into Python lists; validate structure-only so jsonschema does not
        # require/iterate the bulk data.
        validate_yaml(yaml_input, "plant/wind_energy_system", array_data=True)
        system_dat = load_yaml(Path(yaml_input), nc_data="array")
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
        DensityCompensation,
        PowerCtFunctionList,
        PowerCtTabular,
        SimpleYawModel,
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

        # Use DensityCompensation (wind speed correction before lookup) to
        # match foxes' air density handling: ws *= (rho/rho_ref)^(1/3)
        try:
            yaw_model = SimpleYawModel(exp=2)  # PyWake < 2.6
        except TypeError:
            yaw_model = SimpleYawModel()  # PyWake >= 2.6
        density_models = [yaw_model, DensityCompensation(1.225)]

        this_turbine = WindTurbine(
            name=turbine_dat["name"],
            diameter=rd,
            hub_height=hh,
            powerCtFunction=PowerCtTabular(
                speeds,
                powers,
                power_unit="W",
                ct=cts_int,
                additional_models=density_models,
            ),
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

    if "probability" in resource_ds:
        # windIO histogram resource. Producers that follow the windkit/WAsP
        # convention write ``probability`` as the per-sector *conditional*
        # distribution over wind speed (sums to 1 per wind_direction slice)
        # with the direction marginal in ``sector_probability``; pyWake's
        # ``P`` is the *joint* distribution over (wd, ws).  Multiply the two
        # when both are present; a bare ``probability`` is taken as already
        # joint (matching py_wake.utils.ieawind37_utils).
        if "sector_probability" in resource_ds:
            resource_ds["P"] = (
                resource_ds["probability"] * resource_ds["sector_probability"]
            )
            resource_ds = resource_ds.drop_vars(["probability", "sector_probability"])
        else:
            resource_ds = resource_ds.rename({"probability": "P"})

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

    if "time" in resource_ds.dims:
        # Convert time coordinate to integer indices for GridInterpolator compatibility
        # (string or datetime time coords cannot be interpolated numerically)
        resource_ds = resource_ds.assign_coords(time=np.arange(len(resource_ds.time)))
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
    elif "probability" in wind_resource:
        result = _construct_histogram_site(resource_dat, hub_heights, x_positions)
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
    do_subset = False
    if not times_run.get("all_occurences", True) and "subset" in times_run:
        cases_idx = times_run["subset"]
        do_subset = True

    heights = wind_resource.get("height")
    # ``heights`` may be a list or a numpy array (array-backed resource), so use
    # an explicit boolean rather than the ambiguous truth value of an array.
    has_heights = heights is not None and len(heights) > 0

    def _subset(arr):
        # Only copy when a real time subset is requested.  An all-True boolean
        # index still copies the whole array, so skip it in the common case.
        return arr[cases_idx] if do_subset else arr

    # Helper to get data and dimensions safely.  np.asarray avoids copying when
    # the data is already an ndarray (e.g. an array-backed windIO resource).
    def get_resource_data(var_name):
        data_obj = wind_resource[var_name]
        vals = np.asarray(data_obj["data"])
        dims = data_obj.get("dims", ["time"])
        return vals, dims

    # Extract raw data
    ws_vals, ws_dims = get_resource_data("wind_speed")
    wd_vals, wd_dims = get_resource_data("wind_direction")

    # Apply subsetting
    ws_vals = _subset(ws_vals)
    wd_vals = _subset(wd_vals)

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
        operating = _subset(np.asarray(wind_resource["operating"]["data"])).T
        assert operating.shape[0] == len(x_positions)
    else:
        operating = np.ones((len(x_positions), len(cases_idx)))

    # Handle multi-height interpolation
    additional_heights = []
    hh = list(hub_heights.values())[0]
    site = None

    if len(hub_heights) > 1:
        # Multiple turbine types with differing hub heights.
        if ("wind_turbine" in ws_dims or "wind_turbine" in wd_dims) and not has_heights:
            # Per-turbine ws/wd are given without a vertical profile (e.g. terrain
            # speedups produced by a microscale model, already at each turbine's
            # hub height). Preserve the per-turbine values via dict_to_site rather
            # than averaging across turbines or interpolating a vertical profile.
            site = dict_to_site(wind_resource)
            if "turbulence_intensity" not in wind_resource:
                TI = 0.02
            else:
                TI = _subset(np.asarray(wind_resource["turbulence_intensity"]["data"]))
        elif "wind_turbine" in ws_dims or "wind_turbine" in wd_dims:
            # Both a vertical profile and per-turbine data is ambiguous — we
            # cannot tell whether height or turbine indexes the inflow.
            raise NotImplementedError(
                "A wind resource with both a 'height' profile and a "
                "'wind_turbine' dimension is not supported. Provide one or the "
                "other for mixed hub heights."
            )
        else:
            # Vertical wind profile (a `height` dimension) with mixed hub
            # heights: interpolate the profile to each distinct hub height (plus
            # any requested flow-field z-planes) and build a height-indexed
            # XRSite.  pyWake then assigns every turbine its inflow at its own
            # hub height.
            flow_field_spec = (
                system_dat["attributes"]
                .get("model_outputs_specification", {})
                .get("flow_field", {})
            )
            if (
                "z_planes" in flow_field_spec
                and flow_field_spec["z_planes"] != "hub_heights"
            ):
                additional_heights = flow_field_spec.get("z_planes", {}).get(
                    "z_list", []
                )

            h_levels = sorted(set(hub_heights.values()) | set(additional_heights))

            # Interpolate WS (linear) and WD (vector / circular) to each level.
            speeds = [_interpolate_wind_data(heights, ws, wd, h)[0] for h in h_levels]
            dirs = [_interpolate_wind_dir(heights, wd, h) for h in h_levels]
            n_time = np.asarray(speeds[0]).shape[0]

            # Interpolate TI to each level (floored), else a constant default.
            if "turbulence_intensity" not in wind_resource:
                ti_levels = [np.full(n_time, 0.02) for _ in h_levels]
            else:
                ti_data = _subset(
                    np.asarray(wind_resource["turbulence_intensity"]["data"])
                )
                ti_levels = [
                    _interpolate_with_min(heights, ti_data, h, min_val=0.02)
                    for h in h_levels
                ]

            data_vars = {
                "WS": (["h", "time"], np.array(speeds)),
                "WD": (["h", "time"], np.array(dirs)),
                "TI": (["h", "time"], np.array(ti_levels)),
                "P": (["time"], np.ones(n_time) / n_time),
            }
            if "density" in wind_resource:
                density_vals = _subset(np.asarray(wind_resource["density"]["data"]))
                density_dims = wind_resource["density"].get("dims", ["time"])
                if "height" in density_dims:
                    data_vars["Air_density"] = (
                        ["h", "time"],
                        np.array(
                            [
                                _interpolate_with_min(heights, density_vals, h, 0.0)
                                for h in h_levels
                            ]
                        ),
                    )
                elif "wind_turbine" in density_dims:
                    data_vars["Air_density"] = (["time"], np.mean(density_vals, axis=1))
                else:
                    data_vars["Air_density"] = (["time"], density_vals)

            site = XRSite(
                xr.Dataset(
                    data_vars=data_vars,
                    coords={"h": h_levels, "time": np.arange(n_time)},
                )
            )

            # Reference inflow arrays (1-D over time).  The height-indexed site
            # is authoritative per turbine, so these only define the flow cases;
            # use the turbine-averaged inflow for a representative reference.
            turbine_types = system_dat["wind_farm"]["layouts"][0]["turbine_types"]
            ordered_hh = list(hub_heights.values())
            level_of = {h: i for i, h in enumerate(h_levels)}
            idx_per_turbine = [level_of[ordered_hh[t]] for t in turbine_types]
            ws = np.mean([speeds[i] for i in idx_per_turbine], axis=0)
            rads = np.deg2rad([dirs[i] for i in idx_per_turbine])
            wd = np.mod(
                np.rad2deg(
                    np.arctan2(
                        np.mean(np.sin(rads), axis=0), np.mean(np.cos(rads), axis=0)
                    )
                ),
                360.0,
            )
            TI = np.mean([ti_levels[i] for i in idx_per_turbine], axis=0)
    else:
        # Single turbine type
        if has_heights:
            ws, wd = _interpolate_wind_data(heights, ws, wd, hh)

        assert len(_subset(np.asarray(times))) == len(ws)
        assert len(wd) == len(ws)

        if "wind_turbine" in ws_dims or "wind_turbine" in wd_dims:
            site = dict_to_site(wind_resource)
        else:
            site = Hornsrev1Site()
            if "density" in wind_resource:
                density_vals = _subset(np.asarray(wind_resource["density"]["data"]))
                site.ds["Air_density"] = (("time",), density_vals)

        # Handle TI
        if "turbulence_intensity" not in wind_resource:
            TI = 0.02
        else:
            TI = _subset(np.asarray(wind_resource["turbulence_intensity"]["data"]))
            if has_heights:
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


def _coord_values(coord):
    """Extract a windIO coordinate as a float array.

    windIO permits coordinates either as plain arrays or as
    ``{"data": ..., "dims": ...}`` mappings.
    """
    if isinstance(coord, dict):
        coord = coord["data"]
    return np.asarray(coord, dtype=float)


def _subsector_wind_directions(wd_raw, n_subsector):
    """Split each wind-direction sector into ``n_subsector`` sub-directions.

    Strips a trailing 360° wrap-around value first.  Refinement requires at
    least 4 sectors (equidistant, full circle); otherwise the sector centers
    are returned unchanged.
    """
    wd_sectors = np.asarray(wd_raw, dtype=float)
    if len(wd_sectors) > 1 and np.isclose(wd_sectors[-1], 360.0):
        wd_sectors = wd_sectors[:-1]
    if n_subsector > 1 and len(wd_sectors) >= 4:
        n_sectors = len(wd_sectors)
        sector_width = 360.0 / n_sectors
        subsector_width = sector_width / n_subsector
        offsets = np.linspace(
            -sector_width / 2 + subsector_width / 2,
            sector_width / 2 - subsector_width / 2,
            n_subsector,
        )
        return np.sort(
            (wd_sectors[:, np.newaxis] + offsets[np.newaxis, :]).ravel() % 360
        )
    return wd_sectors


def _construct_histogram_site(resource_dat, hub_heights, x_positions, n_subsector=5):
    """Construct site from a (wind_direction × wind_speed) histogram resource.

    Internal helper for construct_site().

    The windIO resource carries ``probability`` on a shared wind-speed bin
    grid; dict_to_site() builds the joint pyWake ``P`` from it (multiplying
    in ``sector_probability`` when present).  Flow cases: pyWake cannot
    interpolate a ws-dependent ``P`` onto foreign wind speeds, so the
    histogram bin centers are used verbatim; wind directions are refined to
    sub-sectors (as in the Weibull path) — XRSite interpolates ``P``
    circularly over wd and rescales by the sub-sector width.

    Parameters
    ----------
    resource_dat : dict
        Energy resource dictionary from windIO.
    hub_heights : dict
        Mapping of turbine type names to hub heights.
    x_positions : list
        Turbine x positions (for operating array sizing).
    n_subsector : int
        Number of sub-directions per wind direction sector.  Default 5
        (matching the Weibull path).
    """
    wind_resource = resource_dat["wind_resource"]

    ws = _coord_values(wind_resource["wind_speed"])
    wd = _subsector_wind_directions(
        _coord_values(wind_resource["wind_direction"]), n_subsector
    )

    if "turbulence_intensity" in wind_resource:
        # dict_to_site carries the resource TI into site.ds, which
        # run_simulation prefers over this entry; kept for interface parity
        # with the Weibull path.
        TI = wind_resource["turbulence_intensity"]["data"]
    else:
        TI = 0.06  # default when TI is absent from wind resource

    return {
        "site": dict_to_site(wind_resource),
        "ws": ws,
        "wd": wd,
        "TI": TI,
        "timeseries": False,
        "operating": np.ones((len(x_positions), 1)),
        "additional_heights": [],
        "cases_idx": np.ones(1).astype(bool),
    }


def _construct_weibull_site(resource_dat, hub_heights, x_positions, n_subsector=5):
    """Construct site from Weibull distribution data.

    Internal helper for construct_site().

    Parameters
    ----------
    resource_dat : dict
        Energy resource dictionary from windIO.
    hub_heights : dict
        Mapping of turbine type names to hub heights.
    x_positions : list
        Turbine x positions (for operating array sizing).
    n_subsector : int
        Number of sub-directions per wind direction sector.  Higher values
        smooth directional wake effects.  Default 5 (matching pywasp).
    """
    wind_resource = resource_dat["wind_resource"]
    A = wind_resource["weibull_a"]
    k = wind_resource["weibull_k"]
    wd_raw = wind_resource["wind_direction"]

    # --- Speedup computation ------------------------------------------------
    # Handle turbine-specific Weibull
    if "wind_turbine" in wind_resource["sector_probability"]["dims"]:
        mean_ws = np.array(A["data"]) * gamma(1 + 1.0 / np.array(k["data"]))
        wt_axis = list(A["dims"]).index("wind_turbine")
        max_mean = np.max(mean_ws, axis=wt_axis, keepdims=True)
        Speedup = mean_ws / max_mean
        wind_resource["Speedup"] = {
            "dims": list(A["dims"]),
            "data": Speedup.tolist(),
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

    # --- Flow case computation -----------------------------------------------
    # When wind_speed is absent from the windIO dict, WIFA computes optimal
    # flow cases: a Speedup-adjusted ws range and sub-sector wd values.
    # When wind_speed IS present, the user has chosen explicit flow cases
    # and both ws and wd are used as-is.
    ws = wind_resource.get("wind_speed", None)
    if ws is None:
        # -- Wind speed range from Weibull + Speedup --------------------------
        A_arr = np.asarray(A["data"], dtype=float)
        k_arr = np.asarray(k["data"], dtype=float)
        # Weibull inverse CDF at 99.9 %: ws = A * (-ln(0.001))^(1/k)
        ws_999 = A_arr * (-np.log(0.001)) ** (1.0 / k_arr)
        ws_max_local = float(np.max(ws_999))
        # Extend for speed-downs so the reference WS grid covers every
        # turbine's distribution after Speedup scaling
        if "Speedup" in wind_resource:
            min_speedup = float(np.min(wind_resource["Speedup"]["data"]))
            ws_max_ref = ws_max_local / max(min_speedup, 0.1)
        else:
            ws_max_ref = ws_max_local
        # Start at the first nonzero bin: a ws=0 reference case carries zero
        # energy for every model, but it is a degenerate flow case for the
        # WeightedSum superposition (Zong), whose convection-velocity iteration
        # divides by the convection speed and is undefined at zero wind speed.
        # Including ws=0 silently corrupts the WeightedSum AEP (collapsing the
        # apparent wake loss); dropping it is harmless for all other models.
        ws = np.arange(0.5, np.ceil(ws_max_ref) + 0.5, 0.5)

        # -- Wind direction sub-sectors ---------------------------------------
        wd = _subsector_wind_directions(wd_raw, n_subsector)
    else:
        # Explicit wind_speed provided: use original wd as-is
        wd = wd_raw

    # --- Site and TI --------------------------------------------------------
    site = dict_to_site(wind_resource)

    if "turbulence_intensity" in wind_resource:
        # Reuse the dataset already materialized inside ``site`` instead of
        # rebuilding it with a second dict_to_netcdf.  dict_to_site renames
        # turbulence_intensity -> TI and height -> h.
        ti_da = site.ds["TI"]
        if "x" in ti_da.dims:
            interpolated_ti = ti_da.interp(x=x_positions, y=x_positions)
            if "h" in interpolated_ti.dims:
                interpolated_ti = interpolated_ti.interp(h=hub_heights["0"])
            TI = np.array(
                [interpolated_ti.isel(x=i, y=i).values for i in range(len(x_positions))]
            )
        else:
            TI = wind_resource["turbulence_intensity"]["data"]
    else:
        TI = 0.06  # default when TI is absent from wind resource

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


def _interpolate_wind_dir(heights, wd, target_height):
    """Interpolate wind direction (degrees) to ``target_height``.

    Uses vector (sin/cos) interpolation so the 0/360 wrap-around is handled
    correctly — plain linear interpolation of the raw degrees would average,
    e.g., 350° and 10° to 180° instead of 0°.
    """
    if heights is None:
        return wd
    rad = np.deg2rad(np.asarray(wd, dtype=float))
    try:
        sin_i = interp1d(heights, np.sin(rad), axis=1, fill_value="extrapolate")(
            target_height
        )
        cos_i = interp1d(heights, np.cos(rad), axis=1, fill_value="extrapolate")(
            target_height
        )
    except ValueError:
        sin_i = interp1d(heights, np.sin(rad).T, axis=1, fill_value="extrapolate")(
            target_height
        )
        cos_i = interp1d(heights, np.cos(rad).T, axis=1, fill_value="extrapolate")(
            target_height
        )
    return np.mod(np.rad2deg(np.arctan2(sin_i, cos_i)), 360.0)


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


def configure_wake_model(
    system_dat,
    rotor_diameter,
    hub_height,
    resource_dat=None,
    turbine_geometries=None,
):
    """Configure the wake model components based on system configuration.

    Args:
        system_dat: System data dictionary
        rotor_diameter: Rotor diameter for FUGA LUT generation
        hub_height: Hub height for FUGA LUT generation
        resource_dat: Optional energy_resource dict; lets FUGA derive its LUT
            roughness/inversion height from the site (z0 from TI, zi from
            ABL_height) instead of using defaults.
        turbine_geometries: Optional list of (rotor_diameter, hub_height) for
            every turbine type; FUGA builds one LUT set per geometry so mixed
            farms interpolate over d_h. Defaults to the single
            (rotor_diameter, hub_height).

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

    wake_model_class, deficit_args, deficit_post_attrs = _configure_deficit_model(
        wind_deficit_data,
        analysis,
        rotor_diameter,
        hub_height,
        resource_dat,
        turbine_geometries,
    )
    deflection_model = _configure_deflection_model(deflection_data)
    turbulence_model = _configure_turbulence_model(turbulence_data)
    superposition_model = _configure_superposition_model(superposition_data)
    rotor_averaging = _configure_rotor_averaging(rotor_avg_data)
    blockage_model = _configure_blockage_model(blockage_data, deficit_args)

    # WeightedSum/CumulativeWakeSum impose two constraints that PyWake otherwise
    # enforces via bare AssertionErrors deep in a run. Fail fast with actionable
    # messages; windIO cannot express these cross-field constraints.
    from py_wake.deficit_models.deficit_model import ConvectionDeficitModel
    from py_wake.rotor_avg_models.rotor_avg_model import NodeRotorAvgModel
    from py_wake.superposition_models import CumulativeWakeSum, WeightedSum

    if isinstance(superposition_model, (WeightedSum, CumulativeWakeSum)):
        # 1. Requires a node-based rotor-averaging model, or None (rotor centre,
        #    which PyWake accepts; an explicit RotorCenter is rejected).
        if rotor_averaging is not None and not isinstance(
            rotor_averaging, NodeRotorAvgModel
        ):
            raise ValueError(
                "WeightedSum/CumulativeWakeSum superposition requires a node "
                "rotor-averaging model (grid/eq_grid/gq_grid/cgi) or 'none'; "
                "center, gaussian_overlap and area_overlap are not node models."
            )
        # 2. Requires a convection-based deficit (carries a convective velocity).
        if not issubclass(wake_model_class, ConvectionDeficitModel):
            raise ValueError(
                f"WeightedSum/CumulativeWakeSum superposition requires a "
                f"ConvectionDeficitModel-based deficit (e.g. Zong2020); "
                f"'{wind_deficit_data['name']}' "
                f"({wake_model_class.__name__}) is not one. SuperGaussian, "
                f"SuperGaussian2023 and GCL do not support it — use Linear."
            )

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
        "deficit_post_attrs": deficit_post_attrs,
        "wake_deficit_key": None,  # Deprecated: kept for API compatibility
        "deflection_model": deflection_model,
        "turbulence_model": turbulence_model,
        "superposition_model": superposition_model,
        "rotor_averaging": rotor_averaging,
        "blockage_model": blockage_model,
        "solver_class": solver_class,
        "solver_args": solver_args,
    }


# --- Fuga LUT generation -----------------------------------------------------
# Fuga is a linearised-RANS wake model that reads a precomputed look-up table
# (LUT). Historically those came from a Windows GUI; pyfuga (conda-forge, pure
# Python) now generates them, so WIFA can build a LUT on the fly for any farm.
#
# Fuga has NO turbulence-intensity input: ambient turbulence enters implicitly
# through the roughness z0 (which sets the neutral shear/mixing) and the
# stability zeta0. We therefore derive z0 from the site's representative TI via
# the same inversion PyWake uses at runtime (z0 = zref * exp(-1/TI), neutral),
# unless the windIO config or the resource supplies z0 directly.


def _fuga_default_lut_dir():
    """Persistent, shared LUT cache dir (override with $WIFA_FUGA_LUT_DIR).

    LUTs are content-addressed by filename, so a single shared dir is safe and
    lets the expensive preLUT stage be reused across runs and farms.
    """
    return Path(
        os.environ.get(
            "WIFA_FUGA_LUT_DIR", Path.home() / ".cache" / "wifa" / "fuga_luts"
        )
    )


def _resource_field_array(resource_dat, key):
    """Finite values of a windIO wind_resource field (dict-with-'data' or array).

    Returns a 1-D numpy array, or None if the field is absent/empty.
    """
    if resource_dat is None:
        return None
    field = resource_dat.get("wind_resource", {}).get(key)
    if field is None:
        return None
    data = np.asarray(field["data"] if isinstance(field, dict) else field, dtype=float)
    data = data[np.isfinite(data)].ravel()
    return data if data.size else None


def _mean_resource_field(resource_dat, key):
    """Finite mean of a windIO wind_resource field, or None."""
    data = _resource_field_array(resource_dat, key)
    return float(np.mean(data)) if data is not None else None


def _fuga_atmosphere(resource_dat, fuga_cfg, hub_height):
    """Resolve (z0, zi, zeta0, ti) for a Fuga LUT from config + site resource.

    Precedence for z0/zi: explicit fuga config > site resource field > default.
    z0 is derived from the site TI (Fuga has no TI knob) when not given.
    """
    from py_wake.utils import fuga_utils

    zeta0 = float(fuga_cfg.get("zeta0", 0.0))

    zi = fuga_cfg.get("zi")
    if zi is None:
        zi = _mean_resource_field(resource_dat, "ABL_height")
    if zi is None:
        zi = 500.0

    ti = _mean_resource_field(resource_dat, "turbulence_intensity")
    z0 = fuga_cfg.get("z0")
    if isinstance(z0, (list, tuple)):
        # An explicit z0 list is a sweep, handled by _fuga_z0_sweep; the single
        # scalar here is only a fallback, so don't treat the list as scalar.
        z0 = None
    if z0 is None:
        z0 = _mean_resource_field(resource_dat, "z0")
    if z0 is None and ti is not None and ti > 0:
        z0 = float(np.ravel(fuga_utils.z0(ti, hub_height, zeta0))[0])
    if z0 is None:
        z0 = 0.03  # open-farmland fallback if neither z0 nor TI is available
    return float(z0), float(zi), zeta0, ti


def _fuga_z0_sweep(resource_dat, fuga_cfg, hub_height, zeta0, z0_single):
    """z0 values for a TI-faithful multi-LUT, spanning the site TI distribution.

    Fuga reads TI off the LUT roughness, so a single mean-TI LUT evaluates the
    wake at loss(mean TI) and misses the low-TI tail that drives the deepest
    wakes. A sweep of LUTs across z0 lets FugaDeficit interpolate z0 = z0(TI)
    per flow case at run time, so the farm loss is integrated over the TI
    distribution instead of taken at its mean (cf. the GCL free-stream-TI gap).

    Returns a sorted list of distinct z0. Falls back to ``[z0_single]`` when TI
    data is unavailable, z0 is pinned in config, or n_z0 <= 1. Out-of-range TI
    is handled by FugaDeficit's bounds='limit' (clamped to the nearest LUT).
    """
    from py_wake.utils import fuga_utils

    if fuga_cfg.get("z0") is not None:
        z0s = fuga_cfg["z0"]
        z0s = z0s if isinstance(z0s, (list, tuple)) else [z0s]
        return sorted({float(z) for z in z0s})

    n = int(fuga_cfg.get("n_z0", 5))
    ti = _resource_field_array(resource_dat, "turbulence_intensity")
    if n <= 1 or ti is None:
        return [z0_single]
    ti = ti[ti > 0]
    if ti.size == 0:
        return [z0_single]
    lo, hi = np.quantile(
        ti, [fuga_cfg.get("ti_qlo", 0.05), fuga_cfg.get("ti_qhi", 0.95)]
    )
    # Clamp to a TI band that keeps the neutral-inversion z0 physical: the
    # mapping z0 = zhub*exp(-1/TI) sends high TI to absurd roughness (TI 0.30 ->
    # z0 ~2.8 m, well outside Fuga's linearisation). The low-TI tail drives the
    # deepest, most TI-sensitive wakes, so cover it; high-TI cases saturate to
    # shallow wakes and clamp to the roughest LUT via bounds='limit'. Band
    # [0.03, 0.18] keeps z0 in ~[1e-5, 0.3] m.
    ti_lo = float(fuga_cfg.get("ti_min", 0.03))
    ti_hi = float(fuga_cfg.get("ti_max", 0.18))
    lo, hi = float(np.clip(lo, ti_lo, ti_hi)), float(np.clip(hi, ti_lo, ti_hi))

    def _z0(ti_val):
        return round(float(np.ravel(fuga_utils.z0(ti_val, hub_height, zeta0))[0]), 8)

    if hi <= lo:
        # Whole TI distribution sits at/over a clamp bound -> a single LUT at
        # the clamped TI (still physical), not the unclamped mean-TI z0.
        return [_z0(hi)]
    return sorted({_z0(t) for t in np.linspace(lo, hi, n)})


def _ensure_fuga_luts(
    *,
    folder,
    zeta0,
    nkz0,
    nbeta,
    geometries,
    z0_list,
    zi,
    lut_vars,
    nx,
    ny,
    zlow=None,
    zhigh=None,
    dx=None,
    dy=None,
    n_cpu=None,
):
    """Generate/reuse a LUT for every (geometry, z0) pair; return the path list.

    All LUTs share the costly preLUT (which depends only on zeta0/nkz0/nbeta),
    so extra z0 values and turbine geometries add only the cheap per-LUT stage.
    FugaDeficit/FugaMultiLUTDeficit interpolate the resulting set over d_h
    (turbine geometry) and z0 (per-flow-case TI).

    zlow/zhigh/dx/dy default to each geometry's own hub height and D/4, D/16.
    Mixed-geometry layouts must pass a shared zlow/zhigh (spanning every hub
    height) and a shared dx/dy so FugaMultiLUTDeficit can merge the LUTs onto
    one z/x/y grid: merging single-height LUTs at different hub heights turns
    the whole table NaN (xarray cannot interpolate a size-1 z axis), which
    surfaced as zero power at every cross-type waked turbine.
    """
    paths = []
    for diameter, zhub in geometries:
        for z0 in z0_list:
            paths.append(
                _ensure_fuga_lut(
                    folder=folder,
                    zeta0=zeta0,
                    nkz0=nkz0,
                    nbeta=nbeta,
                    diameter=diameter,
                    zhub=zhub,
                    z0=z0,
                    zi=zi,
                    lut_vars=lut_vars,
                    nx=nx,
                    ny=ny,
                    zlow=zlow,
                    zhigh=zhigh,
                    dx=dx,
                    dy=dy,
                    n_cpu=n_cpu,
                )
            )
    return paths


def _ensure_fuga_lut(
    *,
    folder,
    zeta0,
    nkz0,
    nbeta,
    diameter,
    zhub,
    z0,
    zi,
    lut_vars,
    nx,
    ny,
    zlow=None,
    zhigh=None,
    dx=None,
    dy=None,
    n_cpu=None,
):
    """Generate (or reuse a cached) Fuga LUT; return its path.

    The LUT filename encodes every physical/grid parameter, so an existing file
    with the right name is a valid cache hit. pyfuga reuses the costly preLUT
    stage (which depends only on zeta0/nkz0/nbeta) across geometries.
    """
    from pyfuga import get_luts
    from pyfuga.paths import get_luts_path

    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    # pyfuga's own defaults; pass explicitly so the cache-probe path built by
    # get_luts_path matches the filename get_luts actually writes.
    if dx is None:
        dx = diameter / 4
    if dy is None:
        dy = diameter / 16
    # zlow == zhigh == zhub -> single hub-height level (the cheap path).
    if zlow is None:
        zlow = zhub
    if zhigh is None:
        zhigh = zhub
    lut_vars = list(lut_vars)
    lut_path = get_luts_path(
        folder,
        zeta0,
        nkz0,
        nbeta,
        diameter,
        zhub,
        z0,
        zi,
        zlow,
        zhigh,
        lut_vars,
        nx,
        ny,
        dx,
        dy,
    )
    if not lut_path.exists():
        get_luts(
            folder=folder,
            zeta0=zeta0,
            nkz0=nkz0,
            nbeta=nbeta,
            diameter=diameter,
            zhub=zhub,
            z0=z0,
            zi=zi,
            zlow=zlow,
            zhigh=zhigh,
            lut_vars=lut_vars,
            nx=nx,
            ny=ny,
            dx=dx,
            dy=dy,
            n_cpu=n_cpu,
        )
    return str(lut_path)


def _configure_deficit_model(
    wind_deficit_data,
    analysis,
    rotor_diameter,
    hub_height,
    resource_dat=None,
    turbine_geometries=None,
):
    """Configure the wind deficit model.

    Returns:
        tuple: (wake_model_class, deficit_args, deficit_post_attrs) where
            deficit_post_attrs is a dict of attributes to set on the built
            deficit after construction (e.g. TurbOPark's WS_key).
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
    from py_wake.deficit_models.noj import NOJDeficit, NOJLocalDeficit, TurboNOJDeficit

    model_name = wind_deficit_data["name"]
    normalized = _normalize_name(model_name)

    wind_deficit_cfg = analysis.get("wind_deficit_model", {})
    # Honor the windIO use_effective_ws flag (local vs free-stream inflow at the
    # waking turbine); deficits that don't accept it pop it below (NOJDeficit).
    deficit_args = {"use_effective_ws": wind_deficit_cfg.get("use_effective_ws", True)}
    wake_expansion = wind_deficit_cfg.get("wake_expansion_coefficient", {})

    GAUSSIAN_MODELS = {
        "bastankhah2014": BastankhahGaussianDeficit,
        "niayifar2016": NiayifarGaussianDeficit,
        "zong2020": ZongGaussianDeficit,
        "carbajofuertes2018": CarbajofuertesGaussianDeficit,
    }
    # Models that accept a=[k_a, k_b] instead of k (scalar)
    A_PARAM_MODELS = {"niayifar2016", "zong2020", "carbajofuertes2018"}
    # Deficits that expose a use_effective_ti param (TI-dependent expansion/width).
    # NOJLocalDeficit (Jensen) accepts it too: with a=[k_a, k_b] it references
    # effective TI, so honoring free_stream_ti lets a no-turbulence config use
    # ambient TI. GCLDeficit also accepts it (GCLLocal sets use_effective_ti=True).
    # Bastankhah2014, free-stream NOJDeficit and FUGA do not.
    TI_CAPABLE = {
        "jensen",
        "nojlocaldeficit",
        "niayifar2016",
        "carbajofuertes2018",
        "zong2020",
        "turbopark",
        "gcl",
        "supergaussian",
        "supergaussian2023",
        "eddyviscosity",
    }

    # windIO convention: k = k_a + k_b * TI (k_a constant, k_b multiplies TI).
    # PyWake's a-parametrized deficits compute k = a[0] * TI + a[1], so the
    # windIO pair maps to a = [k_b, k_a]. Scalar-k deficits take k_a and cannot
    # represent a nonzero k_b.
    if normalized in ("jensen", "nojlocaldeficit"):
        wake_model_class = NOJLocalDeficit
        if "k_a" in wake_expansion or "k_b" in wake_expansion:
            deficit_args["a"] = [
                wake_expansion.get("k_b", 0) or 0,
                wake_expansion.get("k_a", 0) or 0,
            ]

    elif normalized in ("jensen1983", "nojdeficit"):
        wake_model_class = NOJDeficit
        deficit_args.pop("use_effective_ws", None)
        # NOJDeficit takes a scalar k. windIO's wake_expansion_coefficient has
        # no scalar `k` field, so accept k_a (the constant) as well as `k`.
        if "k" in wake_expansion:
            deficit_args["k"] = wake_expansion["k"]
        elif "k_a" in wake_expansion:
            deficit_args["k"] = wake_expansion["k_a"]
        if wake_expansion.get("k_b"):
            warnings.warn(
                f"{model_name} takes a constant wake expansion k (= k_a); "
                f"the TI coefficient k_b={wake_expansion['k_b']} is ignored."
            )

    elif normalized in GAUSSIAN_MODELS:
        wake_model_class = GAUSSIAN_MODELS[normalized]
        if normalized in A_PARAM_MODELS:
            # Niayifar, Zong, Carbajofuertes: k = k_a + k_b*TI -> a=[k_b, k_a]
            if "k" in wake_expansion:
                warnings.warn(
                    f"{model_name} uses k_a/k_b (k = k_a + k_b*TI) for wake "
                    f"expansion, not scalar k. Scalar 'k' is ignored."
                )
            if "k_a" in wake_expansion or "k_b" in wake_expansion:
                if "k_b" not in wake_expansion:
                    warnings.warn(
                        f"k_b not specified for {model_name}, defaulting to 0 "
                        f"(TI-independent wake expansion)"
                    )
                deficit_args["a"] = [
                    wake_expansion.get("k_b", 0) or 0,
                    wake_expansion.get("k_a", 0) or 0,
                ]
        else:
            # Bastankhah2014 uses k (scalar) = the windIO constant k_a
            if "k_a" in wake_expansion:
                deficit_args["k"] = wake_expansion["k_a"]
            elif "k" in wake_expansion:
                deficit_args["k"] = wake_expansion["k"]
            if wake_expansion.get("k_b"):
                warnings.warn(
                    f"{model_name} takes a constant wake expansion k (= k_a); "
                    f"the TI coefficient k_b={wake_expansion['k_b']} is ignored."
                )
        # ceps maps to the deficit's near-wake epsilon coefficient. Bastankhah,
        # Niayifar and Carbajofuertes name it `ceps`; Zong names it `eps_coeff`.
        if "ceps" in wind_deficit_cfg:
            if normalized == "zong2020":
                deficit_args["eps_coeff"] = wind_deficit_cfg["ceps"]
            else:
                deficit_args["ceps"] = wind_deficit_cfg["ceps"]

    elif normalized == "supergaussian":
        wake_model_class = BlondelSuperGaussianDeficit2020

    elif normalized == "supergaussian2023":
        wake_model_class = BlondelSuperGaussianDeficit2023

    elif normalized == "turbopark":
        wake_model_class = TurboGaussianDeficit
        # Canonical Nygaard (2022) recipe (py_wake.literature.turbopark): a Mirror
        # ground model and ctlim=0.96 as constructor args; the WS_key='WS_jlk'
        # attribute (scale the deficit by the downstream turbine's ambient WS) is
        # applied post-construction via deficit_post below.
        from py_wake.ground_models.ground_models import Mirror
        from py_wake.superposition_models import SquaredSum

        deficit_args["groundModel"] = Mirror(superpositionModel=SquaredSum())
        deficit_args["ctlim"] = 0.96

    elif normalized == "turbonoj":
        wake_model_class = TurboNOJDeficit
        if "A" in wind_deficit_cfg:
            deficit_args["A"] = wind_deficit_cfg["A"]

    elif normalized == "gcl":
        wake_model_class = GCLDeficit

    elif normalized == "eddyviscosity":
        # Only on pyWake's unmerged EV branch (cj_add_eddy_viscosity_model);
        # import lazily so the released pyWake keeps working without it.
        from py_wake.deficit_models.eddy_viscosity import EddyViscosityDeficitModel

        wake_model_class = EddyViscosityDeficitModel
        # EV (Ainslie 1988) references deficits to the free-stream wind speed
        # and combines them with MaxSum (the WindFarmer convention), so unlike
        # the analytic deficits use_effective_ws defaults to False.
        deficit_args["use_effective_ws"] = wind_deficit_cfg.get(
            "use_effective_ws", False
        )

    elif normalized == "bastankhah2016":
        raise NotImplementedError(
            "Bastankhah2016 is not available in PyWake. Use flow_model 'foxes', "
            "or choose Bastankhah2014/Zong2020 for PyWake."
        )

    elif normalized == "fuga":
        wake_model_class = FugaDeficit
        # FugaDeficit reads a LUT instead of an analytic expansion; it takes no
        # use_effective_ws (it always uses the free-stream-referenced deficit).
        deficit_args.pop("use_effective_ws", None)
        fuga_cfg = wind_deficit_cfg.get("fuga", {}) or {}
        z0_single, zi, zeta0, _ti = _fuga_atmosphere(resource_dat, fuga_cfg, hub_height)
        # A z0 sweep across the site TI distribution + a LUT per turbine
        # geometry; FugaDeficit interpolates z0 (per-flow-case TI) and d_h at
        # run time. Degenerates to a single LUT for one geometry + n_z0<=1.
        z0_list = _fuga_z0_sweep(resource_dat, fuga_cfg, hub_height, zeta0, z0_single)
        geometries = turbine_geometries or [(rotor_diameter, hub_height)]
        # Dedupe: two turbine types with the same geometry share one LUT, and
        # duplicate d_h coordinates would break FugaMultiLUTDeficit's merge.
        geometries = list(dict.fromkeys((float(d), float(h)) for d, h in geometries))
        hub_heights = sorted({h for _, h in geometries})
        diameters = sorted({d for d, _ in geometries})
        # Mixed hub heights: every LUT must span all hub heights (a source
        # turbine's wake is evaluated at each target's hub height), and mixed
        # diameters need one shared x/y grid; otherwise FugaMultiLUTDeficit's
        # merge yields NaN deficits -> zero power at cross-type waked turbines.
        mixed_grid = {}
        if len(hub_heights) > 1:
            mixed_grid["zlow"] = hub_heights[0]
            mixed_grid["zhigh"] = hub_heights[-1]
            # Interpolate the merged LUTs at exactly the hub heights; pyfuga's
            # log-spaced z levels differ per z0, so without this the z-union
            # across a z0 sweep would reintroduce NaN edge cells.
            deficit_args["z_lst"] = hub_heights
        if len(diameters) > 1:
            # Finest natural resolution; identical x/y coords across LUTs.
            mixed_grid["dx"] = diameters[0] / 4
            mixed_grid["dy"] = diameters[0] / 16
        lut_paths = _ensure_fuga_luts(
            folder=fuga_cfg.get("cache_dir", _fuga_default_lut_dir()),
            zeta0=zeta0,
            nkz0=int(fuga_cfg.get("nkz0", 8)),
            nbeta=int(fuga_cfg.get("nbeta", 32)),
            geometries=geometries,
            z0_list=z0_list,
            zi=zi,
            lut_vars=fuga_cfg.get("lut_vars", ["UL"]),
            nx=int(fuga_cfg.get("nx", 2048)),
            ny=int(fuga_cfg.get("ny", 512)),
            n_cpu=fuga_cfg.get("n_cpu"),
            **mixed_grid,
        )
        # Single LUT -> plain path; multiple -> list (FugaDeficit globs/lists).
        deficit_args["LUT_path"] = lut_paths[0] if len(lut_paths) == 1 else lut_paths

    else:
        raise NotImplementedError(f"Wake model '{model_name}' is not supported")

    # TI reference: windIO carries this as the nested wake_expansion_coefficient
    # .free_stream_ti flag (foxes-compatible). PyWake's deficits expose the
    # inverse use_effective_ti param (use_effective_ti = not free_stream_ti),
    # but only the TI-dependent deficits accept it.
    if normalized in TI_CAPABLE and "free_stream_ti" in wake_expansion:
        deficit_args["use_effective_ti"] = not wake_expansion["free_stream_ti"]

    # Axial induction: windIO's axial_induction_model maps to PyWake's ct2a
    # (1D -> ct2a_mom1d, Madsen -> ct2a_madsen). Honor it on every deficit that
    # accepts a ct2a parameter; without this the deficit silently keeps its
    # ct2a_madsen default, so a "1D" request was previously ignored on the
    # pyWake path.
    import inspect

    from py_wake.deficit_models.utils import ct2a_madsen, ct2a_mom1d

    axial = analysis.get("axial_induction_model")
    if axial is not None:
        ct2a_fn = {"1d": ct2a_mom1d, "madsen": ct2a_madsen}.get(_normalize_name(axial))
        if (
            ct2a_fn is not None
            and "ct2a" in inspect.signature(wake_model_class.__init__).parameters
        ):
            deficit_args["ct2a"] = ct2a_fn

    # Attributes set on the deficit *after* construction (not constructor
    # kwargs), applied by run_simulation.
    deficit_post = {}
    if normalized == "turbopark":
        deficit_post["WS_key"] = "WS_jlk"

    return wake_model_class, deficit_args, deficit_post


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
        c = turbulence_data.get("c")
        if c is not None:
            # A paper's calibration (e.g. Niayifar 2016, Zong 2020): the
            # literature CrespoHernandez uses 1D induction and SqrMaxSum
            # added-turbulence summation alongside the calibrated coefficients.
            from py_wake.deficit_models.utils import ct2a_mom1d
            from py_wake.superposition_models import SqrMaxSum

            return CrespoHernandez(
                c=list(c),
                ct2a=ct2a_mom1d,
                addedTurbulenceSuperpositionModel=SqrMaxSum(),
            )
        return CrespoHernandez()
    if normalized == "gcl":
        return GCLTurbulence()
    if normalized in ("quartonandainslie", "modifiedquartonandainslie"):
        # Only on pyWake's unmerged EV branch; the Hassan (1992) modified
        # variant is the one the EV bundle uses. It carries its own added-TI
        # combination (SqrMaxSum default); windIO's ti_superposition is not
        # consulted on the pyWake path.
        from py_wake.turbulence_models.quarton_and_ainslie import (
            ModifiedQuartonAndAinslieTurbulenceModel,
        )

        return ModifiedQuartonAndAinslieTurbulenceModel()
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
    if normalized == "vector":
        raise NotImplementedError(
            "Vector superposition is foxes-only; not available in PyWake."
        )
    raise NotImplementedError(f"Superposition model '{name}' is not supported")


def _configure_rotor_averaging(rotor_avg_data):
    """Configure the rotor averaging model."""
    from py_wake.rotor_avg_models import (
        AreaOverlapAvgModel,
        CGIRotorAvg,
        EqGridRotorAvg,
        GaussianOverlapAvgModel,
        GQGridRotorAvg,
        GridRotorAvg,
        PolarGridRotorAvg,
        RotorCenter,
    )

    name = rotor_avg_data["name"]
    normalized = _normalize_name(name)

    if normalized == "none":
        # No rotor-averaging model. PyWake's Weighted superposition accepts this
        # (rotor centre) but rejects an explicit RotorCenter; the Zong (2020)
        # literature model uses None.
        return None
    if normalized == "center":
        return RotorCenter()
    # "grid" is the canonical windIO name; "avgdeficit" is a deprecated alias
    if normalized in ("grid", "avgdeficit"):
        return GridRotorAvg()
    if normalized == "gaussianoverlap":
        return GaussianOverlapAvgModel()
    if normalized == "areaoverlap":
        return AreaOverlapAvgModel()
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
    if normalized == "simplifiedgaussian":
        # Only on pyWake's unmerged EV branch: the LUT-based line average
        # across the rotor that the EV bundle pairs with its deficit.
        from py_wake.rotor_avg_models.simplified_gaussian_rotor_average_model import (
            SimplifiedGaussianRotorAverageModel,
        )

        return SimplifiedGaussianRotorAverageModel()
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
    from py_wake.superposition_models import LinearSum

    name = blockage_data["name"]
    if name is None:
        return None

    normalized = _normalize_name(name)
    if normalized == "none":
        return None

    # The analytic blockage models are calibrated without ground effects and
    # allow flow through the ground; ground_mirror enforces the slip boundary
    # condition with an image rotor (unlike wake models, which are calibrated
    # including ground effects and must not be mirrored).
    ground_model = None
    if blockage_data.get("ground_mirror", False):
        if normalized == "fuga":
            warnings.warn(
                "blockage_model.ground_mirror is ignored for FUGA: the Fuga "
                "LUTs come from a linearized RANS solver that already includes "
                "the ground."
            )
        else:
            from py_wake.ground_models.ground_models import Mirror

            ground_model = Mirror()

    # Models that take no constructor arguments
    SIMPLE_BLOCKAGE_MODELS = {
        "selfsimilaritydeficit": SelfSimilarityDeficit,
        "rankinehalfbody": RankineHalfBody,
        "rathmann": Rathmann,
        "vortexcylinder": VortexCylinder,
        "vortexdipole": VortexDipole,
        "hybridinduction": HybridInduction,
    }

    # Sum blockage deficits (across turbines, and real+image when mirrored)
    # linearly: potential-flow induction superposes linearly, and pyWake
    # otherwise falls back to the *wake* superposition model — SquaredSum
    # (e.g. the TurbOPark recipe) asserts on the speed-up (negative deficit)
    # regions every blockage model produces.
    blockage_superposition = LinearSum()

    if normalized == "selfsimilaritydeficit2020":
        return SelfSimilarityDeficit2020(
            ss_alpha=blockage_data.get("ss_alpha", 0.8888888888888888),
            groundModel=ground_model,
            superpositionModel=blockage_superposition,
        )
    if normalized in SIMPLE_BLOCKAGE_MODELS:
        return SIMPLE_BLOCKAGE_MODELS[normalized](
            groundModel=ground_model,
            superpositionModel=blockage_superposition,
        )
    if normalized == "fuga":
        return FugaDeficit(deficit_args["LUT_path"], z_lst=deficit_args.get("z_lst"))
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
    # Build deficit model. groundModel comes from deficit_args when a model needs
    # a specific one (e.g. TurbOPark's Mirror); otherwise the deficit's own
    # default (None) applies.
    deficit_args = dict(wake_config["deficit_args"])
    deficit_args.setdefault("groundModel", None)
    deficit_model = wake_config["wake_model_class"](
        rotorAvgModel=wake_config["rotor_averaging"],
        **deficit_args,
    )
    for attr, value in wake_config.get("deficit_post_attrs", {}).items():
        setattr(deficit_model, attr, value)

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
    # Collect every turbine geometry so FUGA can build a LUT set per type
    # (mixed farms interpolate over d_h); the first one drives single-type paths.
    if "turbines" in farm_dat:
        geometries = [
            (
                farm_dat["turbines"]["rotor_diameter"],
                farm_dat["turbines"]["hub_height"],
            )
        ]
    else:
        geometries = [
            (t["rotor_diameter"], t["hub_height"])
            for t in farm_dat["turbine_types"].values()
        ]
    rd, first_hh = geometries[0]

    wake_config = configure_wake_model(
        system_dat, rd, first_hh, resource_dat, geometries
    )

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
