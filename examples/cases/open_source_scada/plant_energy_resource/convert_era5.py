#!/usr/bin/env python
# ERA5 Wind Data Conversion Script
# Converts ERA5 wind dataset to simplified format with time, wind_speed, and wind_direction

import numpy as np
import xarray as xr


def convert_era5_wind_data(input_file, output_file, height_level=3):
    """
    Convert ERA5 wind dataset to simplified format with time, wind_speed, and wind_direction.

    Parameters:
    -----------
    input_file : str
        Path to the input ERA5 netCDF file
    output_file : str
        Path to save the output netCDF file
    height_level : int, optional
        Index of the height level to extract (default: 3, which is typically 75m)
    """
    # Load the ERA5 dataset
    print(f"Loading ERA5 dataset from {input_file}...")
    era5_data = xr.load_dataset(input_file)

    # Print available heights for reference
    print(f"Available heights: {era5_data.height.values} meters")
    print(
        f"Using height level {height_level}: {era5_data.height.values[height_level]} meters"
    )

    # Create a new dataset with the desired variables
    # Extract wind speed and direction at the specified height level
    wind_speed = era5_data.WS.sel(height=era5_data.height.values[height_level])
    wind_direction = era5_data.WD.sel(height=era5_data.height.values[height_level])

    # Create a new time coordinate that's just sequential integers
    new_time = np.arange(len(era5_data.time))

    # Create the new dataset
    new_data = xr.Dataset(
        data_vars={
            "wind_speed": ("time", wind_speed.values),
            "wind_direction": ("time", wind_direction.values),
        },
        coords={
            "time": new_time,
        },
    )

    # Save the new dataset
    print(f"Saving converted dataset to {output_file}...")
    new_data.to_netcdf(output_file)
    print("Conversion complete!")

    # Display a summary of the converted dataset
    print("\nConverted dataset summary:")
    print(new_data)

    return new_data


if __name__ == "__main__":
    # Example usage
    input_file = "./ERA5_wind_timeseries_flow_testpark_v2.nc"
    output_file = "./simplified_wind_data.nc"

    # Default is height index 3 (75m), but you can change this as needed
    # Heights in the original dataset are: 10.0, 25.0, 50.0, 75.0, 100.0, 150.0, 200.0, 250.0 meters
    convert_era5_wind_data(input_file, output_file, height_level=3)
