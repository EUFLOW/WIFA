import argparse
import os
import warnings
from pathlib import Path

import numpy as np
import xarray as xr
import yaml
from scipy.interpolate import interp1d
from scipy.special import gamma
from windIO import dict_to_netcdf, load_yaml
from windIO import validate as validate_yaml

from wifa._optional import require

# Define default values for wind_deficit_model parameters
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


def _farm_layout(farm_dat):
    """Return the first layout of a farm (layouts may be a dict or a list)."""
    layouts = farm_dat["layouts"]
    return layouts[0] if isinstance(layouts, list) else layouts


def _farm_turbine_specs(farm_dat):
    """Return (turbine_specs, per_position_type_idx) for one farm.

    turbine_specs is the list of raw windIO turbine dicts used by the farm;
    per_position_type_idx maps each turbine position to an index into that list.
    """
    layout = _farm_layout(farm_dat)
    n_positions = len(layout["coordinates"]["x"])
    if "turbines" in farm_dat:
        specs = [farm_dat["turbines"]]
        per_pos = [0] * n_positions
    else:
        type_map = farm_dat["turbine_types"]
        keys = list(type_map.keys())
        specs = [type_map[k] for k in keys]
        if "turbine_types" in layout:
            # Layout entries are keys into the turbine_types mapping (windIO
            # allows arbitrary keys, e.g. 1-based). Resolution rules:
            # 1. If every entry resolves as a mapping key (exact match first;
            #    a string/number coercion like 1 vs "1" is consulted only when
            #    it does not shadow a different real key), use key semantics.
            # 2. Otherwise, if every entry is a valid 0-based integer index,
            #    interpret positionally (the windIO schema also describes the
            #    entries as integer indices) and warn with the per-entry
            #    assignment so nothing is reinterpreted silently.
            # 3. Otherwise raise. YAML booleans are rejected outright: they
            #    would alias integer keys/indices 0 and 1.
            entries = list(layout["turbine_types"])
            if any(isinstance(k, (bool, np.bool_)) for k in entries):
                raise ValueError(
                    f"Layout turbine_types of farm "
                    f"'{farm_dat.get('name', '?')}' contains boolean entries "
                    f"({entries!r}); booleans are neither mapping keys nor "
                    "positional indices"
                )
            exact = {k: i for i, k in enumerate(keys)}
            coerced = {}
            for i, k in enumerate(keys):
                s = str(k)
                if s in exact and exact[s] != i:
                    continue  # a different key claims this literal
                if s not in coerced:
                    coerced[s] = i

            def _resolve(entry):
                idx = exact.get(entry)
                if idx is None:
                    idx = coerced.get(str(entry))
                return idx

            resolved = [_resolve(k) for k in entries]
            positional_ok = all(
                isinstance(k, (int, np.integer)) and 0 <= k < len(specs)
                for k in entries
            )
            if all(idx is not None for idx in resolved):
                per_pos = resolved
            elif positional_ok:
                per_pos = [int(k) for k in entries]
                assignment = {
                    str(k): specs[i].get("name", i) for k, i in zip(entries, per_pos)
                }
                warnings.warn(
                    f"Farm '{farm_dat.get('name', '?')}': layout turbine_types "
                    "entries do not all match the turbine_types mapping keys; "
                    f"interpreting them as 0-based positional indices "
                    f"({assignment})"
                )
            else:
                bad = sorted(
                    {k for k, idx in zip(entries, resolved) if idx is None},
                    key=str,
                )
                raise ValueError(
                    f"Layout turbine types {bad!r} of farm "
                    f"'{farm_dat.get('name', '?')}' neither match its "
                    f"turbine_types mapping keys ({keys}) nor form valid "
                    "0-based positional indices"
                )
        elif len(specs) == 1:
            per_pos = [0] * n_positions
        else:
            raise ValueError(
                f"Farm '{farm_dat.get('name', '?')}' has multiple turbine_types "
                "but its layout does not specify a turbine type per position"
            )
    if len(per_pos) != n_positions:
        raise ValueError(
            f"Farm '{farm_dat.get('name', '?')}' has {n_positions} turbine "
            f"positions but {len(per_pos)} layout turbine_types entries"
        )
    return specs, per_pos


def create_turbine(turbine_dat):
    """Create a PyWake WindTurbine from a windIO turbine dict."""
    from py_wake.wind_turbines import WindTurbine
    from py_wake.wind_turbines.power_ct_functions import (
        PowerCtFunctionList,
        PowerCtTabular,
    )

    hh = turbine_dat["hub_height"]
    rd = turbine_dat["rotor_diameter"]

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
        raise ValueError("Missing Cp_curve or power_curve in turbine performance data")

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
    return this_turbine


def _specs_equal(a, b):
    """Compare two turbine spec fragments, tolerating numpy arrays and NaNs.

    Plain dict equality raises "ambiguous truth value" when dict-input callers
    supply curves as numpy arrays, and NaN placeholders must compare equal to
    themselves so identical specs are recognized as identical.
    """
    if a is b:
        return True
    if isinstance(a, np.ndarray) and a.ndim == 0:
        a = a.item()
    if isinstance(b, np.ndarray) and b.ndim == 0:
        b = b.item()
    if isinstance(a, dict) and isinstance(b, dict):
        return a.keys() == b.keys() and all(_specs_equal(a[k], b[k]) for k in a)
    a_is_seq = isinstance(a, (list, tuple, np.ndarray))
    b_is_seq = isinstance(b, (list, tuple, np.ndarray))
    if a_is_seq or b_is_seq:
        if not (a_is_seq and b_is_seq):
            return False
        a_list, b_list = list(a), list(b)
        return len(a_list) == len(b_list) and all(
            _specs_equal(x, y) for x, y in zip(a_list, b_list)
        )
    try:
        if np.isnan(a) and np.isnan(b):
            return True
    except (TypeError, ValueError):
        pass
    return bool(a == b)


def _build_multifarm_turbines(farms):
    """Build turbine objects spanning one or more farms.

    Turbine specs that appear in several farms are merged when they are
    identical; reusing a turbine name for a different spec is an error.

    Returns:
        turbine: WindTurbine (single type) or WindTurbines (multi-type)
        turbine_types: 0 or array length sum(N_i) — global type index per turbine
        hub_heights: dict mapping global type index (str) to hub height
        farm_slices: list of slice() objects, one per farm, indexing the global turbine axis
        rotor_diameter: rotor diameter of the first turbine type
    """
    from py_wake.wind_turbines import WindTurbines

    merged_specs = []
    type_indices = []
    farm_slices = []
    cursor = 0

    for farm_dat in farms:
        specs, per_pos = _farm_turbine_specs(farm_dat)

        local_to_global = []
        for spec in specs:
            for global_idx, seen in enumerate(merged_specs):
                if seen["name"] == spec["name"]:
                    if not _specs_equal(seen, spec):
                        raise ValueError(
                            f"Turbine '{spec['name']}' is defined differently "
                            "in different farms"
                        )
                    local_to_global.append(global_idx)
                    break
            else:
                merged_specs.append(spec)
                local_to_global.append(len(merged_specs) - 1)

        type_indices.extend(local_to_global[i] for i in per_pos)
        farm_slices.append(slice(cursor, cursor + len(per_pos)))
        cursor += len(per_pos)

    turbines = [create_turbine(spec) for spec in merged_specs]
    hub_heights = {str(i): spec["hub_height"] for i, spec in enumerate(merged_specs)}
    rotor_diameter = merged_specs[0]["rotor_diameter"]

    if len(turbines) == 1:
        return turbines[0], 0, hub_heights, farm_slices, rotor_diameter
    return (
        WindTurbines.from_WindTurbine_lst(turbines),
        np.asarray(type_indices),
        hub_heights,
        farm_slices,
        rotor_diameter,
    )


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
    print("making site with ", resource_ds)
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
    from py_wake.examples.data.hornsrev1 import Hornsrev1Site
    from py_wake.site import XRSite
    from windIO import dict_to_netcdf

    # Get flow field bounds from config or site boundaries
    boundaries = system_dat["site"]["boundaries"]["polygons"][0]
    WFXLB = np.min(boundaries["x"])
    WFXUB = np.max(boundaries["x"])
    WFYLB = np.min(boundaries["y"])
    WFYUB = np.max(boundaries["y"])

    # Override with explicit flow field bounds if specified
    WFXLB = get_flow_field_param(system_dat, "xlb", WFXLB)
    WFXUB = get_flow_field_param(system_dat, "xub", WFXUB)
    WFYLB = get_flow_field_param(system_dat, "ylb", WFYLB)
    WFYUB = get_flow_field_param(system_dat, "yub", WFYUB)
    WFDX = get_flow_field_param(system_dat, "dx", (WFXUB - WFXLB) / 100)
    WFDY = get_flow_field_param(system_dat, "dy", (WFYUB - WFYLB) / 100)

    flow_bounds = {
        "xlb": WFXLB,
        "xub": WFXUB,
        "ylb": WFYLB,
        "yub": WFYUB,
        "dx": WFDX,
        "dy": WFDY,
    }

    # Determine site type and construct accordingly
    if "time" in resource_dat["wind_resource"]:
        # Timeseries site
        result = _construct_timeseries_site(
            system_dat, resource_dat, hub_heights, x_positions
        )
        result["flow_bounds"] = flow_bounds
        return result

    elif "weibull_k" in resource_dat["wind_resource"]:
        # Weibull distribution site
        result = _construct_weibull_site(resource_dat, hub_heights, x_positions)
        result["flow_bounds"] = flow_bounds
        return result

    else:
        # Simple probability-based site
        ws = resource_dat["wind_resource"]["wind_speed"]
        wd = resource_dat["wind_resource"]["wind_direction"]
        site = dict_to_site(resource_dat["wind_resource"])
        TI = resource_dat["wind_resource"]["turbulence_intensity"]["data"]

        return {
            "site": site,
            "ws": ws,
            "wd": wd,
            "TI": TI,
            "timeseries": False,
            "operating": np.ones((len(x_positions), 1)),
            "additional_heights": [],
            "cases_idx": np.ones(1).astype(bool),
            "flow_bounds": flow_bounds,
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
    output_spec = system_dat["attributes"].get("model_outputs_specification", {})
    if "run_configuration" in output_spec:
        run_config = output_spec["run_configuration"]
        if "times_run" in run_config and not run_config["times_run"].get(
            "all_occurences", True
        ):
            if "subset" in run_config["times_run"]:
                cases_idx = run_config["times_run"]["subset"]

    heights = wind_resource.get("height")
    n_cases = len(np.arange(len(times))[cases_idx])

    # Helper to get time-subset data and dimensions safely, honoring the
    # declared dims order rather than assuming time is the leading axis
    def get_resource_data(var_name, default_dims=("time",)):
        data_obj = wind_resource[var_name]
        vals = np.array(data_obj["data"])
        dims = list(data_obj.get("dims", default_dims))
        if heights is not None and vals.ndim == len(dims) + 1 and "height" not in dims:
            # dims underspecified (e.g. 2-D data with no dims declared): the
            # extra axis is taken as height only when its length matches the
            # height coordinate and cannot be the turbine axis instead;
            # anything ambiguous must declare dims explicitly
            n_h = len(heights)
            n_wt = len(x_positions)
            if vals.shape[-1] == n_h and n_h != n_wt:
                dims = dims + ["height"]
            elif vals.shape[0] == n_h and vals.shape[-1] != n_h and n_h != n_wt:
                dims = ["height"] + dims
            else:
                raise ValueError(
                    f"Cannot infer the dimensions of wind_resource variable "
                    f"'{var_name}' with shape {vals.shape}: declare its "
                    "'dims' explicitly"
                )
        if "time" in dims:
            time_sel = np.arange(len(times))[cases_idx]
            vals = np.take(vals, time_sel, axis=dims.index("time"))
        return vals, dims

    def get_density_series(target_height):
        density_vals, density_dims = _mean_over_turbines(*get_resource_data("density"))
        return _interp_along_height(density_vals, density_dims, heights, target_height)

    # Extract raw data (time-subset)
    ws_vals, ws_dims = get_resource_data("wind_speed")
    wd_vals, wd_dims = get_resource_data("wind_direction")

    # Prepare reference arrays - average across turbines if turbine-specific
    ws, ws_dims_eff = _mean_over_turbines(ws_vals, ws_dims)

    if "wind_turbine" in wd_dims:
        # Vector mean for direction to handle 360/0 boundary
        rads = np.deg2rad(wd_vals)
        wt_axis = wd_dims.index("wind_turbine")
        mean_sin = np.mean(np.sin(rads), axis=wt_axis)
        mean_cos = np.mean(np.cos(rads), axis=wt_axis)
        wd = np.mod(np.rad2deg(np.arctan2(mean_sin, mean_cos)), 360)
        wd_dims_eff = [d for d in wd_dims if d != "wind_turbine"]
    else:
        wd = wd_vals
        wd_dims_eff = wd_dims

    # Handle operating status
    if "operating" in wind_resource:
        op_vals, op_dims = get_resource_data("operating", ("time", "wind_turbine"))
        if "wind_turbine" in op_dims and op_dims.index("wind_turbine") != 0:
            op_vals = op_vals.T
        operating = op_vals
        assert operating.shape[0] == len(x_positions)
    else:
        operating = np.ones((len(x_positions), n_cases))

    # Handle multi-height interpolation
    additional_heights = []
    hh = first_hh = list(hub_heights.values())[0]
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
            speeds.append(
                _interp_along_height(ws, ws_dims_eff, heights, hh, min_val=0.0)
            )
            dirs.append(_interp_direction_along_height(wd, wd_dims_eff, heights, hh))

        ws, wd = speeds[-1], dirs[-1]

        # Handle TI interpolation: one entry per unique height in `seen`
        if "turbulence_intensity" in wind_resource:
            TI_data, ti_dims = get_resource_data("turbulence_intensity")
            if "wind_turbine" in ti_dims:
                warnings.warn(
                    "Per-turbine turbulence_intensity is averaged across "
                    "turbines for the multi-hub-height time-series site; the "
                    "farm-mean TI is used at every height"
                )
            TI_data, ti_dims = _mean_over_turbines(TI_data, ti_dims)
            for hh in seen:
                TIs.append(
                    _interp_along_height(TI_data, ti_dims, heights, hh, min_val=0.02)
                )
            TI = TIs[-1]
        else:
            TI = 0.02
            TIs = [np.full(np.shape(speeds[0]), TI) for _ in seen]

        data_vars = {
            "WS": (["h", "time"], np.array(speeds)),
            "WD": (["h", "time"], np.array(dirs)),
            "TI": (["h", "time"], np.array(TIs)),
            "P": 1,
        }
        if "density" in wind_resource:
            data_vars["Air_density"] = (["time"], get_density_series(first_hh))
        site = XRSite(
            xr.Dataset(
                data_vars=data_vars,
                coords={"h": seen, "time": np.arange(np.shape(speeds[0])[0])},
            )
        )
    else:
        # Single turbine type
        ws = _interp_along_height(ws, ws_dims_eff, heights, hh, min_val=0.0)
        wd = _interp_direction_along_height(wd, wd_dims_eff, heights, hh)

        assert len(np.array(times)[cases_idx]) == len(ws)
        assert len(wd) == len(ws)

        if "wind_turbine" in ws_dims or "wind_turbine" in wd_dims:
            site = dict_to_site(wind_resource)
        else:
            site = Hornsrev1Site()
            if "density" in wind_resource:
                site.ds["Air_density"] = (("time",), get_density_series(hh))

        # Handle TI (kept per-turbine here; PyWake accepts turbine-resolved TI)
        if "turbulence_intensity" not in wind_resource:
            TI = 0.02
        else:
            TI, ti_dims = get_resource_data("turbulence_intensity")
            TI = _interp_along_height(TI, ti_dims, heights, hh, min_val=0.02)

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


def _mean_over_turbines(vals, dims):
    """Average a resource variable over its declared wind_turbine dim.

    Returns the (possibly reduced) values and the dims list without the
    wind_turbine entry, so downstream axis lookups stay consistent.
    """
    if "wind_turbine" in dims:
        vals = np.mean(vals, axis=dims.index("wind_turbine"))
        dims = [d for d in dims if d != "wind_turbine"]
    return vals, dims


def _interp_along_height(vals, dims, heights, target_height, min_val=None):
    """Interpolate a resource variable to target_height along its declared
    height dim; variables without a height dim pass through unchanged.

    min_val is a floor applied only when target_height lies outside the
    declared height range (i.e. only to extrapolated values); in-range data
    is never clamped.
    """
    if "height" not in dims:
        return vals
    if heights is None:
        raise ValueError(
            "A wind_resource variable declares a 'height' dim but the "
            "wind_resource has no 'height' coordinate"
        )
    out = interp1d(heights, vals, axis=dims.index("height"), fill_value="extrapolate")(
        target_height
    )
    if min_val is not None and not (
        np.min(heights) <= target_height <= np.max(heights)
    ):
        out = np.maximum(out, min_val)
    return out


def _interp_direction_along_height(vals, dims, heights, target_height):
    """Interpolate wind direction to target_height via its sine/cosine
    components, so the 0/360 wraparound cannot produce spurious directions."""
    if "height" not in dims:
        return vals
    rads = np.deg2rad(np.asarray(vals, dtype=float))
    sin_int = _interp_along_height(np.sin(rads), dims, heights, target_height)
    cos_int = _interp_along_height(np.cos(rads), dims, heights, target_height)
    return np.mod(np.rad2deg(np.arctan2(sin_int, cos_int)), 360)


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
    from py_wake.deficit_models import SelfSimilarityDeficit2020
    from py_wake.deficit_models.fuga import FugaDeficit
    from py_wake.deficit_models.gaussian import (
        BastankhahGaussianDeficit,
        BlondelSuperGaussianDeficit2020,
        TurboGaussianDeficit,
    )
    from py_wake.deficit_models.noj import NOJLocalDeficit
    from py_wake.deflection_models import JimenezWakeDeflection
    from py_wake.rotor_avg_models import GridRotorAvg, RotorCenter
    from py_wake.superposition_models import LinearSum, SquaredSum
    from py_wake.turbulence_models import (
        CrespoHernandez,
        STF2005TurbulenceModel,
        STF2017TurbulenceModel,
    )
    from py_wake.wind_farm_models import All2AllIterative, PropagateDownwind

    analysis = system_dat["attributes"]["analysis"]

    # Get model configurations with defaults
    wind_deficit_data = get_with_default(analysis, "wind_deficit_model", DEFAULTS)
    deflection_data = get_with_default(analysis, "deflection_model", DEFAULTS)
    turbulence_data = get_with_default(analysis, "turbulence_model", DEFAULTS)
    superposition_data = get_with_default(analysis, "superposition_model", DEFAULTS)
    rotor_avg_data = get_with_default(analysis, "rotor_averaging", DEFAULTS)
    blockage_data = get_with_default(analysis, "blockage_model", DEFAULTS)

    # Configure wind deficit model
    deficit_args = {"use_effective_ws": True}
    wake_deficit_key = None

    print("Running deficit ", wind_deficit_data)

    wake_model_class, deficit_args, wake_deficit_key = _configure_deficit_model(
        wind_deficit_data, analysis, rotor_diameter, hub_height, deficit_args
    )

    print("deficit args ", deficit_args)

    # Configure deflection model
    deflection_model = _configure_deflection_model(deflection_data)

    # Configure turbulence model
    turbulence_model = _configure_turbulence_model(turbulence_data)

    # Configure superposition model
    superposition_model = _configure_superposition_model(superposition_data)
    print("using superposition ", superposition_data)

    # Configure rotor averaging
    rotor_averaging = _configure_rotor_averaging(rotor_avg_data)

    # Configure blockage model
    blockage_model = _configure_blockage_model(blockage_data, deficit_args)

    # Determine solver based on blockage
    solver_args = {}
    if blockage_model is not None:
        solver_class = All2AllIterative
        solver_args["blockage_deficitModel"] = blockage_model
    else:
        solver_class = PropagateDownwind

    return {
        "wake_model_class": wake_model_class,
        "deficit_args": deficit_args,
        "wake_deficit_key": wake_deficit_key,
        "deflection_model": deflection_model,
        "turbulence_model": turbulence_model,
        "superposition_model": superposition_model,
        "rotor_averaging": rotor_averaging,
        "blockage_model": blockage_model,
        "solver_class": solver_class,
        "solver_args": solver_args,
    }


def _configure_deficit_model(
    wind_deficit_data, analysis, rotor_diameter, hub_height, deficit_args
):
    """Configure the wind deficit model.

    Returns:
        tuple: (wake_model_class, deficit_args, wake_deficit_key)
    """
    from py_wake.deficit_models.fuga import FugaDeficit
    from py_wake.deficit_models.gaussian import (
        BastankhahGaussianDeficit,
        BlondelSuperGaussianDeficit2020,
        TurboGaussianDeficit,
    )
    from py_wake.deficit_models.noj import NOJLocalDeficit

    wake_deficit_key = None
    model_name = wind_deficit_data["name"]

    if model_name == "Jensen":
        wake_model_class = NOJLocalDeficit
        wake_expansion = analysis.get("wind_deficit_model", {}).get(
            "wake_expansion_coefficient", {}
        )
        if "k_b" in wake_expansion:
            k_a = wake_expansion.get("k_a", 0)
            k_b = wake_expansion["k_b"]
            deficit_args["a"] = [k_a, k_b]

    elif model_name.lower() == "bastankhah2014":
        wake_model_class = BastankhahGaussianDeficit
        wake_expansion = analysis.get("wind_deficit_model", {}).get(
            "wake_expansion_coefficient", {}
        )
        if "k_b" in wake_expansion:
            deficit_args["k"] = wake_expansion["k_b"]
        elif "k" in wake_expansion:
            deficit_args["k"] = wake_expansion["k"]
        if "ceps" in analysis.get("wind_deficit_model", {}):
            deficit_args["ceps"] = analysis["wind_deficit_model"]["ceps"]

    elif model_name == "SuperGaussian":
        wake_model_class = BlondelSuperGaussianDeficit2020

    elif model_name == "TurbOPark":
        wake_model_class = TurboGaussianDeficit

    elif model_name.upper() == "FUGA":
        wake_model_class = FugaDeficit
        from pyfuga import get_luts

        lut = get_luts(
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

    # Handle k/k2 format conversion
    if "k2" in deficit_args:
        k = deficit_args.pop("k")
        k2 = deficit_args.pop("k2")
        deficit_args["a"] = [k2, k]

    return wake_model_class, deficit_args, wake_deficit_key


def _configure_deflection_model(deflection_data):
    """Configure the wake deflection model."""
    from py_wake.deflection_models import JimenezWakeDeflection

    name = deflection_data["name"].lower()
    if name == "none":
        return None
    elif name == "jimenez":
        return JimenezWakeDeflection(beta=deflection_data["beta"])
    else:
        raise NotImplementedError(
            f"Deflection model '{deflection_data['name']}' is not supported"
        )


def _configure_turbulence_model(turbulence_data):
    """Configure the turbulence model."""
    from py_wake.turbulence_models import (
        CrespoHernandez,
        STF2005TurbulenceModel,
        STF2017TurbulenceModel,
    )

    name = turbulence_data["name"].upper()
    if turbulence_data["name"].lower() == "none":
        return None
    elif name == "STF2005":
        return STF2005TurbulenceModel(c=[turbulence_data["c1"], turbulence_data["c2"]])
    elif name == "STF2017":
        return STF2017TurbulenceModel(c=[turbulence_data["c1"], turbulence_data["c2"]])
    elif name == "CRESPOHERNANDEZ":
        return CrespoHernandez()
    else:
        raise NotImplementedError(
            f"Turbulence model '{turbulence_data['name']}' is not supported"
        )


def _configure_superposition_model(superposition_data):
    """Configure the superposition model."""
    from py_wake.superposition_models import LinearSum, SquaredSum

    name = superposition_data["ws_superposition"].lower()
    if name == "linear":
        return LinearSum()
    elif name == "squared":
        return SquaredSum()
    else:
        raise NotImplementedError(
            f"Superposition model '{superposition_data['ws_superposition']}' is not supported"
        )


def _configure_rotor_averaging(rotor_avg_data):
    """Configure the rotor averaging model."""
    from py_wake.rotor_avg_models import GridRotorAvg, RotorCenter

    name = rotor_avg_data["name"].lower()
    if name == "center":
        print("Using Center Average")
        return RotorCenter()
    elif name == "avg_deficit":
        return GridRotorAvg()
    else:
        raise NotImplementedError(
            f"Rotor averaging model '{rotor_avg_data['name']}' is not supported"
        )


def _configure_blockage_model(blockage_data, deficit_args):
    """Configure the blockage model."""
    from py_wake.deficit_models import SelfSimilarityDeficit2020
    from py_wake.deficit_models.fuga import FugaDeficit

    name = blockage_data["name"]
    if name == "None" or name is None:
        return None
    elif name == "SelfSimilarityDeficit2020":
        return SelfSimilarityDeficit2020(ss_alpha=blockage_data["ss_alpha"])
    elif name.upper() == "FUGA":
        return FugaDeficit(deficit_args["LUT_path"])
    else:
        raise ValueError(f"Unknown blockage model: {name}")


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
    print("Running ", wake_config["wake_model_class"], wake_config["deficit_args"])
    deficit_model = wake_config["wake_model_class"](
        rotorAvgModel=wake_config["rotor_averaging"],
        groundModel=None,
        **wake_config["deficit_args"],
    )

    if wake_config["wake_deficit_key"]:
        deficit_model.WS_key = wake_config["wake_deficit_key"]

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
    aep = sim_res.aep(normalize_probabilities=not site_data["timeseries"]).sum()
    print("aep is ", aep, "GWh")

    # Calculate per-turbine AEP on the same normalization basis as the total
    aep_per_turbine = (
        sim_res.aep(normalize_probabilities=not site_data["timeseries"])
        .sum(["time"] if site_data["timeseries"] else ["ws", "wd"])
        .to_numpy()
    )

    print(sim_res)

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

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Write turbine outputs if requested
    output_spec = system_dat["attributes"].get("model_outputs_specification", {})
    if "turbine_outputs" in output_spec:
        sim_res_formatted = sim_res[["Power", "WS_eff"]].rename(
            {"Power": "power", "WS_eff": "effective_wind_speed", "wt": "turbine"}
        )
        turbine_nc_filename = str(
            output_spec.get("turbine_outputs", {}).get(
                "turbine_nc_filename", "PowerTable.nc"
            )
        )
        turbine_nc_filepath = Path(output_dir) / turbine_nc_filename
        sim_res_formatted.to_netcdf(turbine_nc_filepath)

    # Flow field handling
    flow_map = _generate_flow_field(
        sim_res, system_dat, site_data, hub_heights, flow_bounds
    )

    if flow_map:
        flow_map = flow_map[["WS_eff", "TI_eff"]].rename(
            {
                "h": "z",
                "WS_eff": "wind_speed",
                "TI_eff": "turbulence_intensity",
            }
        )
        flow_map.to_netcdf(Path(output_dir) / "FarmFlow.nc")

    # Write YAML output
    _write_yaml_output(output_dir)

    return sim_results["aep"]


def _generate_flow_field(sim_res, system_dat, site_data, hub_heights, flow_bounds):
    """Generate flow field data if requested.

    Returns:
        Flow map xarray or None
    """
    output_spec = system_dat["attributes"].get("model_outputs_specification", {})
    timeseries = site_data["timeseries"]

    WFXLB, WFXUB = flow_bounds["xlb"], flow_bounds["xub"]
    WFYLB, WFYUB = flow_bounds["ylb"], flow_bounds["yub"]
    WFDX, WFDY = flow_bounds["dx"], flow_bounds["dy"]

    flow_map = None

    if "flow_field" in output_spec and not timeseries:
        flow_map = sim_res.flow_box(
            x=np.arange(WFXLB, WFXUB + WFDX, WFDX),
            y=np.arange(WFYLB, WFYUB + WFDY, WFDY),
            h=list(hub_heights.values()),
        )

        # Warn if user requests unsupported outputs
        requested_vars = output_spec["flow_field"].get("output_variables", [])
        if any(
            var not in ["velocity_u", "turbulence_intensity"] for var in requested_vars
        ):
            warnings.warn("PyWake can only output velocity_u and turbulence_intensity")

    elif "flow_field" in output_spec and timeseries:
        flow_field_spec = output_spec["flow_field"]
        if flow_field_spec.get("report") is not False:
            z_list = flow_field_spec.get("z_list", sorted(list(hub_heights.values())))
            flow_map = sim_res.flow_box(
                x=np.arange(WFXLB, WFXUB + WFDX, WFDX),
                y=np.arange(WFYLB, WFYUB + WFDY, WFDY),
                h=z_list,
                time=sim_res.time.values,
            )

    return flow_map


def _write_yaml_output(output_dir):
    """Write the output YAML file with include directives."""
    data = {
        "wind_energy_system": "INCLUDE_YAML_PLACEHOLDER",
        "power_table": "INCLUDE_POWER_TABLE_PLACEHOLDER",
        "flow_field": "INCLUDE_FLOW_FIELD_PLACEHOLDER",
    }

    output_yaml_name = Path(output_dir) / "output.yaml"
    with open(output_yaml_name, "w") as file:
        yaml.dump(data, file, default_flow_style=False, allow_unicode=True)

    # Replace placeholders with include directives
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

    with open(output_yaml_name, "w") as file:
        file.write(yaml_content)


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
        float: Total AEP in GWh (single farm), or
        list[float]: AEP in GWh per farm, in input order (multi-farm input)
    """
    # Step 1: Load and validate configuration
    require("py_wake")

    system_dat, output_dir = load_and_validate_config(yaml_input, output_dir)

    # Step 2: Create turbine objects (multi-farm aware)
    farm_entry = system_dat["wind_farm"]
    multi_farm = isinstance(farm_entry, list)
    farms = farm_entry if multi_farm else [farm_entry]

    (
        turbine,
        turbine_types,
        hub_heights,
        farm_slices,
        rotor_diameter,
    ) = _build_multifarm_turbines(farms)

    # Get turbine positions across all farms
    x, y = [], []
    for farm_dat in farms:
        coords = _farm_layout(farm_dat)["coordinates"]
        x.extend(coords["x"])
        y.extend(coords["y"])

    # Step 3: Construct site
    resource_dat = system_dat["site"]["energy_resource"]
    site_data = construct_site(system_dat, resource_dat, hub_heights, x)
    site = site_data["site"]

    # Step 4: Configure wake model
    # First turbine type's dimensions are used for FUGA LUT generation if needed
    first_hh = list(hub_heights.values())[0]
    wake_config = configure_wake_model(system_dat, rotor_diameter, first_hh)

    # Step 5: Run simulation
    sim_results = run_simulation(
        site, turbine, wake_config, site_data, x, y, turbine_types
    )

    # Step 6: Generate outputs
    aep = generate_outputs(sim_results, system_dat, site_data, hub_heights, output_dir)

    if multi_farm:
        per_turbine = sim_results["aep_per_turbine"]
        return [float(np.sum(per_turbine[s])) for s in farm_slices]
    return aep


def run():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_yaml", help="The input yaml file")
    args = parser.parse_args()

    run_pywake(args.input_yaml)


if __name__ == "__main__":
    run()
