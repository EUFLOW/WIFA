from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List

import numpy as np
import windIO
from windIO import load_yaml

from wifa._optional import require

if TYPE_CHECKING:
    from floris import FlorisModel


def run_floris(yaml_input):
    """
    Run FLORIS model based on windIO YAML input and extract specified outputs.

    Parameters
    ----------
    yaml_input: str, Path, or dict
        Path to the windIO YAML input file or a dictionary representing the input data.

    Notes
    -----
    WindIO is not yet suppoeted by FLOIS main branch.
    Currently refer to https://github.com/lejeunemax/floris/tree/windIO.

    Some features may not be fully supported, these include:
    - Subsetting of wind conditions (times_run, wind_speeds_run, directions_run)
    - Only "power", "thrust_coefficient", "axial_induction", and "turbulence_intensity"
      are currently supported as turbine output variables.
    - Some advanced windIO properties may not be fully supported. These include:
      blockage_model, rotor_averaging, and axial_induction_model.
    """
    require("floris", "floris")

    from floris import FlorisModel
    from floris.read_windio import TrackedDict

    if isinstance(yaml_input, (str, Path)):
        windio_dict = windIO.load_yaml(yaml_input)
    elif isinstance(yaml_input, dict):
        windio_dict = yaml_input
    else:
        raise TypeError(
            f"yaml_input must be a file path or dict, got {type(yaml_input)}"
        )

    fmodel = FlorisModel.from_windio(windio_dict)
    fmodel.run()

    windio_dict["_context"] = ""

    with TrackedDict.from_parent(windio_dict, "attributes") as attrs:
        attrs.untrack("analysis")
        attrs.untrack("flow_model")

        out_specs = attrs["model_outputs_specification"]
        run_config = out_specs["run_configuration"]

        slice_selection = None
        times_run: Dict = run_config.get("times_run", None)
        wind_speeds_run: Dict = run_config.get("wind_speeds_run", None)
        directions_run: Dict = run_config.get("directions_run", None)

        if times_run:
            all_occurences = times_run.get("all_occurences", False)

            if all_occurences:
                slice_selection = slice(None)
            else:
                slice_selection = times_run.get("subset", None)

        elif wind_speeds_run and directions_run:
            if not wind_speeds_run.get("all_values", False):
                raise NotImplementedError(
                    f"Wind speed and direction subsets are not yet supported, got {wind_speeds_run.name} {wind_speeds_run}"
                )
            if not directions_run.get("all_values", False):
                raise NotImplementedError(
                    f"Wind speed and direction subsets are not yet supported, got {directions_run.name} {directions_run}"
                )

        output_dir = out_specs.get("output_folder", ".")
        Path(output_dir).mkdir(parents=True, exist_ok=True)

        if "turbine_outputs" in out_specs:
            _read_turbines(fmodel, output_dir, out_specs, slice_selection)

        if "flow_field" in out_specs:
            _read_fields(fmodel, output_dir, out_specs, slice_selection)


def _read_turbines(fmodel: "FlorisModel", output_dir, out_specs, slice_selection):
    """
    Extract turbine output data from FLORIS model and write to NetCDF file.

    Args:
        fmodel: FlorisModel instance that has been run
        output_dir: Output directory path
        out_specs: TrackedDict containing turbine_outputs configuration from out_specs["turbine_outputs"]
        slice_selection: Time slice selection for filtering output

    Returns:
        xarray.Dataset containing turbine output data
    """
    import xarray as xr

    turbine_outputs = out_specs["turbine_outputs"]

    output_vars: List[str] = turbine_outputs["output_variables"]
    turbine_nc_filename = turbine_outputs.get(
        "turbine_nc_filename", "turbine_outputs.nc"
    )

    # Make a copy to avoid modifying the original list
    requested_vars = output_vars.copy()
    data_vars: Dict[str, tuple[List[str], Any]] = {}
    coords = {}

    # Determine the wind data type to understand output structure
    wind_data_type = None
    if fmodel.wind_data is not None:
        from floris.wind_data import WindRose, WindRoseWRG, WindTIRose

        if isinstance(fmodel.wind_data, WindTIRose):
            wind_data_type = "wind_ti_rose"
        elif isinstance(fmodel.wind_data, (WindRose, WindRoseWRG)):
            wind_data_type = "wind_rose"

    # Define handlers for each variable type
    # Format: variable_name -> (extractor_func, is_coord)
    variable_handlers = {
        "time": (lambda: fmodel.get_metadata("wind_resource", {}).get("time"), True),
        "x": (lambda: fmodel.layout_x, False),
        "y": (lambda: fmodel.layout_y, False),
        "z": (lambda: fmodel.get_turbine_layout(z=True)[2], False),
        "power": (lambda: fmodel.get_turbine_powers(), False),
        "thrust_coefficient": (lambda: fmodel.get_turbine_thrust_coefficients(), False),
        "axial_induction": (lambda: fmodel.get_turbine_ais(), False),
        "turbulence_intensity": (lambda: fmodel.get_turbine_TIs(), False),
    }

    # Process each requested variable
    unhandled_vars = []
    for var_name in requested_vars:
        if var_name in variable_handlers:
            try:
                extractor, is_coord = variable_handlers[var_name]
                data = extractor()

                if data is not None:
                    data = np.asarray(data)

                    # Determine dimensions based on data shape
                    if is_coord:
                        # Coordinate variables are 1D
                        if data.ndim == 0:
                            data = np.array([data])
                        coords[var_name] = data
                    else:
                        # Data variables - infer dimensions based on shape and wind data type
                        if data.ndim == 1:
                            # 1D: turbine dimension only
                            data_vars[var_name] = (["turbine"], data)
                        elif data.ndim == 2:
                            # 2D: (time, turbine) - TimeSeries or flattened
                            data_vars[var_name] = (["time", "turbine"], data)
                        elif data.ndim == 3:
                            # 3D: (wind_direction, wind_speed, turbine) - WindRose
                            data_vars[var_name] = (
                                ["wind_direction", "wind_speed", "turbine"],
                                data,
                            )
                        elif data.ndim == 4:
                            # 4D: (wind_direction, wind_speed, turbulence_intensity, turbine) - WindTIRose
                            data_vars[var_name] = (
                                [
                                    "wind_direction",
                                    "wind_speed",
                                    "turbulence_intensity",
                                    "turbine",
                                ],
                                data,
                            )
                        else:
                            print(
                                f"Warning: Unexpected data shape for '{var_name}': {data.shape}; skipping."
                            )
                            continue
                else:
                    print(f"Warning: '{var_name}' data not available; skipping.")
            except Exception as e:
                print(f"Warning: Failed to extract '{var_name}': {e}; skipping.")
        else:
            unhandled_vars.append(var_name)

    # Warn about unhandled variables
    if unhandled_vars:
        print(f"Warning: Unhandled output variables: {unhandled_vars}")

    # Build coordinates based on data structure
    final_coords = {}

    # Add turbine coordinate
    n_turbines = len(fmodel.layout_x)
    final_coords["turbine"] = np.arange(n_turbines)

    # Determine primary dimension and add appropriate coordinates
    if "time" in coords and coords["time"] is not None:
        # TimeSeries case - use time coordinate
        final_coords["time"] = coords["time"]
        # Rename findex to time in data_vars
        renamed_data_vars = {}
        for var_name, (dims, data) in data_vars.items():
            new_dims = ["time" if d == "findex" else d for d in dims]
            renamed_data_vars[var_name] = (new_dims, data)
        data_vars = renamed_data_vars
    elif wind_data_type in ["wind_rose", "wind_ti_rose"]:
        # WindRose case - add wind direction and speed coordinates
        if fmodel.wind_data is not None:
            final_coords["wind_direction"] = fmodel.wind_data.wind_directions
            final_coords["wind_speed"] = fmodel.wind_data.wind_speeds
            if wind_data_type == "wind_ti_rose":
                final_coords[
                    "turbulence_intensity"
                ] = fmodel.wind_data.turbulence_intensities
    else:
        # Fallback: infer findex dimension
        n_findex = None
        for var_name, (dims, data) in data_vars.items():
            if "findex" in dims:
                n_findex = data.shape[dims.index("findex")]
                break
        if n_findex is not None:
            final_coords["findex"] = np.arange(n_findex)

    # Create dataset
    ds = xr.Dataset(data_vars=data_vars, coords=final_coords)

    # Add metadata attributes
    ds.attrs["description"] = "FLORIS turbine output data"
    ds.attrs["source"] = "FLORIS wind farm wake model"

    # Add units and descriptions for known variables
    var_attrs = {
        "x": {"units": "m", "long_name": "turbine x-coordinate"},
        "y": {"units": "m", "long_name": "turbine y-coordinate"},
        "z": {"units": "m", "long_name": "turbine hub height"},
        "power": {"units": "W", "long_name": "turbine power output"},
        "thrust_coefficient": {"units": "-", "long_name": "turbine thrust coefficient"},
        "axial_induction": {
            "units": "-",
            "long_name": "turbine axial induction factor",
        },
        "turbulence_intensity": {
            "units": "-",
            "long_name": "turbulence intensity at turbine",
        },
    }

    for var_name, attrs in var_attrs.items():
        if var_name in ds.data_vars:
            ds[var_name].attrs.update(attrs)

    # Apply time selection if specified
    if slice_selection is not None:
        # Determine which dimension to slice
        if "time" in ds.dims:
            ds = ds.isel({"time": slice_selection})
        elif "findex" in ds.dims:
            ds = ds.isel({"findex": slice_selection})
        elif "wind_direction" in ds.dims:
            # For WindRose, could slice wind_direction or wind_speed
            # This would need more sophisticated handling based on user intent
            print("Warning: slice_selection not applied for WindRose data")

    # Write to NetCDF
    ds.to_netcdf(f"{output_dir}/{turbine_nc_filename}")
    print(f"Turbine output data written to {output_dir}/{turbine_nc_filename}")

    return ds


def _read_fields(fmodel: "FlorisModel", output_dir, out_specs, slice_selection):
    """
    Extract flow field data from FLORIS model and write to NetCDF file.

    Args:
        fmodel: FlorisModel instance that has been run
        output_dir: Output directory path
        out_specs: TrackedDict containing flow_field configuration from out_specs["flow_field"]
        slice_selection: Time slice selection for filtering output

    Returns:
        xarray.Dataset containing flow field data
    """
    import xarray as xr

    flow_field = out_specs["flow_field"]

    # Check if flow field reporting is enabled
    report = flow_field.get("report", True)

    if not report:
        return None

    # Get output configuration
    output_vars: List[str] = flow_field["output_variables"]
    flow_nc_filename = flow_field.get("flow_nc_filename", "flow_field.nc")

    # Get z_planes configuration
    z_planes = flow_field["z_planes"]
    z_sampling = z_planes["z_sampling"]
    xy_sampling = z_planes.get("xy_sampling", "grid")

    # Determine z values to sample
    if z_sampling == "hub_heights":
        z_list = fmodel.core.farm.hub_heights
        # Get unique hub heights
        z_list = np.unique(z_list)
    elif z_sampling == "plane_list":
        z_list = z_planes["z_list"]
    else:
        raise ValueError(f"Unknown z_sampling type: {z_sampling}")

    # Get xy sampling parameters
    if xy_sampling == "grid":
        x_bounds = z_planes["x_bounds"]
        y_bounds = z_planes["y_bounds"]

        # Support either Nx/Ny or dx/dy specification
        if "Nx" in z_planes and "Ny" in z_planes:
            Nx = z_planes["Nx"]
            Ny = z_planes["Ny"]
        elif "dx" in z_planes and "dy" in z_planes:
            dx = z_planes["dx"]
            dy = z_planes["dy"]
            Nx = int((x_bounds[1] - x_bounds[0]) / dx) + 1
            Ny = int((y_bounds[1] - y_bounds[0]) / dy) + 1
        else:
            # Default resolution
            Nx = 200
            Ny = 200

    elif xy_sampling == "original_grid":
        # Use the grid from the model
        x_bounds = None
        y_bounds = None
        Nx = None
        Ny = None
    else:
        raise ValueError(f"Unknown xy_sampling type: {xy_sampling}")

    # Calculate horizontal planes for each z level
    data_vars = {}

    # Prepare coordinate arrays
    n_findex = fmodel.core.flow_field.n_findex
    if slice_selection is not None:
        if isinstance(slice_selection, slice):
            findex_range = range(n_findex)[slice_selection]
        else:
            findex_range = slice_selection
    else:
        findex_range = range(n_findex)

    # For each z level, calculate the horizontal plane
    for z_idx, z_height in enumerate(z_list):
        # Calculate horizontal plane at this height
        # Use keep_inertial_frame=True to avoid rotation for WIFA output
        horizontal_plane = fmodel.calculate_horizontal_plane(
            height=z_height,
            x_resolution=Nx,
            y_resolution=Ny,
            x_bounds=x_bounds,
            y_bounds=y_bounds,
            findex_for_viz=0,  # Will loop over findices if needed
            keep_inertial_frame=True,
        )

        # Extract coordinates from first plane
        if z_idx == 0:
            # Get unique x and y coordinates
            x_coords = np.unique(horizontal_plane.df.x1.values)
            y_coords = np.unique(horizontal_plane.df.x2.values)

            # Initialize data arrays
            if "u" in output_vars:
                u_data = np.zeros(
                    (len(findex_range), len(z_list), len(y_coords), len(x_coords))
                )
            if "v" in output_vars:
                v_data = np.zeros(
                    (len(findex_range), len(z_list), len(y_coords), len(x_coords))
                )
            if "w" in output_vars:
                w_data = np.zeros(
                    (len(findex_range), len(z_list), len(y_coords), len(x_coords))
                )

        # Reshape plane data to grid
        # The dataframe has x1, x2, x3, u, v, w columns
        for t_idx, findex in enumerate(findex_range):
            # Recalculate for this findex if multiple conditions
            if findex > 0:
                horizontal_plane = fmodel.calculate_horizontal_plane(
                    height=z_height,
                    x_resolution=Nx,
                    y_resolution=Ny,
                    x_bounds=x_bounds,
                    y_bounds=y_bounds,
                    findex_for_viz=findex,
                    keep_inertial_frame=True,
                )

            # Reshape to 2D grid (y, x)
            u_grid = horizontal_plane.df.u.values.reshape(len(y_coords), len(x_coords))
            v_grid = horizontal_plane.df.v.values.reshape(len(y_coords), len(x_coords))
            w_grid = horizontal_plane.df.w.values.reshape(len(y_coords), len(x_coords))

            if "u" in output_vars:
                u_data[t_idx, z_idx, :, :] = u_grid
            if "v" in output_vars:
                v_data[t_idx, z_idx, :, :] = v_grid
            if "w" in output_vars:
                w_data[t_idx, z_idx, :, :] = w_grid

    # Create xarray dataset
    coords = {
        "time": list(findex_range),
        "z": z_list,
        "y": y_coords,
        "x": x_coords,
    }

    if "u" in output_vars:
        data_vars["u"] = (("time", "z", "y", "x"), u_data)
    if "v" in output_vars:
        data_vars["v"] = (("time", "z", "y", "x"), v_data)
    if "w" in output_vars:
        data_vars["w"] = (("time", "z", "y", "x"), w_data)

    ds = xr.Dataset(data_vars=data_vars, coords=coords)

    # Add metadata
    ds.attrs["description"] = "FLORIS flow field output"
    ds.x.attrs["units"] = "m"
    ds.y.attrs["units"] = "m"
    ds.z.attrs["units"] = "m"

    if "u" in output_vars:
        ds.u.attrs["units"] = "m/s"
        ds.u.attrs["long_name"] = "x-component of velocity"
    if "v" in output_vars:
        ds.v.attrs["units"] = "m/s"
        ds.v.attrs["long_name"] = "y-component of velocity"
    if "w" in output_vars:
        ds.w.attrs["units"] = "m/s"
        ds.w.attrs["long_name"] = "z-component of velocity"

    # Write to NetCDF
    ds.to_netcdf(f"{output_dir}/{flow_nc_filename}")
    print(f"Flow field data written to {output_dir}/{flow_nc_filename}")

    return ds
