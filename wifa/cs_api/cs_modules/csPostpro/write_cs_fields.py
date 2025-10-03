import sys
from os import path, sep

import netCDF4 as nc
import numpy as np
from netCDF4 import Dataset
from scipy.interpolate import griddata

import wifa.cs_api.cs_modules.csPostpro.cs_postprocess_utils as cs_pp


def get_output_at_plane_and_time(ens, output_varname, origin, normal, time_step):
    """
    evaluate the model for given variable at specific heigh and time.
    """
    #
    data = cs_pp.read_ensight_data_at_time(ens, inst=time_step)
    plane_cut = cs_pp.extract_plane_from_ensight(data, origin=origin, normal=normal)
    output_variable = cs_pp.get_field_from_ensight(plane_cut, output_varname)

    return output_variable


def get_plane_info(ens, origin, normal):
    """
    evaluate the model for given variable at specific heigh and time.
    """
    #
    data = cs_pp.read_ensight_data_at_time(ens, inst=0)
    coupe = cs_pp.extract_plane_from_ensight(data, origin=origin, normal=normal)
    triang = cs_pp.extract_saturne_triangulation(coupe, normal)
    points = data.GetPoints()

    return triang


def find_closest(z, value):
    diff = np.abs(z - value)
    min_diff = np.min(diff)
    closest_indices = np.where(diff == min_diff)
    return closest_indices[0]


postpro_dir = sys.argv[1]
print(postpro_dir)
file_name = sys.argv[2]
print(file_name)
case_name = sys.argv[3]
print(case_name)
case_dir = sys.argv[4]
print(case_dir)

postprocess_heights_file = open(postpro_dir + "/postprocessing_heights.csv", "r")
postprocess_cases_file = open(postpro_dir + "/postprocessing_cases.csv", "r")
postprocess_fields_file = open(postpro_dir + "/postprocessing_fields.csv", "r")

cases = np.atleast_1d(np.genfromtxt(postprocess_cases_file, delimiter=","))
zplot = np.array(np.genfromtxt(postprocess_heights_file, delimiter=","), ndmin=1)
fields = np.genfromtxt(postprocess_fields_file, delimiter=",", dtype=str)

#############CS FIELD DATA###################
# File creation

if sys.argv[5] == "None":
    rootgrp = Dataset(postpro_dir + sep + file_name, "w", format="NETCDF4")

    for i, casei in enumerate(cases):
        print(str(i) + "/" + str(len(cases)), end="\r")
        # TODO: string formatting for id "%05d"
        case_name_id = str(int(1000000 + casei + 1))[1:]
        result_dir = case_name + "_" + case_name_id
        ens = cs_pp.load_ensight_data(
            case_dir + sep + "RESU" + sep + result_dir + sep + "postprocessing",
            "RESULTS_FLUID_DOMAIN.case",
        )
        data = cs_pp.read_ensight_data_at_time(ens, inst=0)
        if i == 0:
            print(case_dir + sep + "RESU" + sep + result_dir + sep + "postprocessing")
            cell_centers = cs_pp.get_cell_centers(data)
            print(cell_centers)
            num_points = cell_centers.GetNumberOfPoints()
            x = np.empty(num_points)
            y = np.empty(num_points)
            z = np.empty(num_points)
            for n in range(num_points):
                x[n], y[n], z[n] = cell_centers.GetPoint(n)

            indices_double = find_closest(z, zplot[0])
            indices = indices_double[::2]  # Take 1 indice sur 2
            x0 = x[indices]
            y0 = y[indices]
            z0 = z[indices]
            npoints = len(indices)
            #
            points = rootgrp.createDimension("points", npoints)
            nc_cases = rootgrp.createDimension("time", len(cases))
            altitudes = rootgrp.createDimension("z", len(zplot))
            rootgrp.createVariable("x", "f8", ("points",))
            rootgrp.createVariable("y", "f8", ("points",))
            rootgrp.createVariable("z", "f8", ("z",))
            rootgrp.createVariable("time", "f8", ("time",))
            rootgrp.variables["x"][:] = x0
            rootgrp.variables["y"][:] = y0
            rootgrp.variables["z"][:] = zplot
            rootgrp.variables["time"][:] = cases
            for field in fields:
                rootgrp.createVariable(
                    field,
                    "f8",
                    (
                        "points",
                        "z",
                        "time",
                    ),
                )
        for j, zj in enumerate(zplot):
            indices_double = find_closest(z, zj)
            indices = indices_double[::2]
            if "wind_speed" or ("wind_direction" in fields):
                velocity = np.empty([len(x0), 3])
                field = data.GetCellData().GetArray("Velocity")
                for m, index in enumerate(indices):
                    velocity[m, :] = field.GetTuple3(index)
                if "wind_speed" in fields:
                    speed = np.sqrt(pow(velocity[:, 0], 2.0) + pow(velocity[:, 1], 2.0))
                    rootgrp.variables["wind_speed"][:, j, i] = speed
                if "wind_direction" in fields:
                    direction = (
                        np.arctan(velocity[:, 1] / velocity[:, 0]) * 360 / (2 * np.pi)
                        + 270
                    )
                    rootgrp.variables["wind_direction"][:, j, i] = direction
            if "pressure" in fields:
                pressure = np.empty(len(x0))
                field = data.GetCellData().GetArray("total_pressure")
                for m, index in enumerate(indices):
                    pressure[m] = field.GetTuple(index)
                rootgrp.variables["pressure"][:, j, i] = pressure
            if "tke" in fields:
                tke = np.empty(len(x0))
                field = data.GetCellData().GetArray("k")
                for m, index in enumerate(indices):
                    tke[m] = field.GetTuple(index)
                rootgrp.variables["tke"][:, j, i] = tke

    rootgrp.close()


#########################################################################
###############     Interpolation Results        #######################
#########################################################################

if sys.argv[5] != "None":
    x_bounds = [float(sys.argv[5]), float(sys.argv[6])]
    y_bounds = [float(sys.argv[7]), float(sys.argv[8])]
    dx = float(sys.argv[9])
    dy = float(sys.argv[10])

    rootgrp = nc.Dataset(postpro_dir + sep + file_name, "w", format="NETCDF4")

    # Define the structured grid
    x_list = np.arange(x_bounds[0], x_bounds[1] + dx, dx)
    y_list = np.arange(y_bounds[0], y_bounds[1] + dy, dy)
    x_grid, y_grid = np.meshgrid(x_list, y_list)

    # Create dimensions for the structured grid
    nx = x_grid.shape[1]
    ny = y_grid.shape[0]
    rootgrp.createDimension("x", nx)
    rootgrp.createDimension("y", ny)
    rootgrp.createDimension("z", len(zplot))
    rootgrp.createDimension("time", len(cases))

    # Create variables for the structured grid
    rootgrp.createVariable("x", "f8", ("x",))
    rootgrp.createVariable("y", "f8", ("y",))
    rootgrp.createVariable("z", "f8", ("z",))
    rootgrp.createVariable("time", "f8", ("time",))
    rootgrp.variables["x"][:] = x_list
    rootgrp.variables["y"][:] = y_list
    rootgrp.variables["z"][:] = zplot
    rootgrp.variables["time"][:] = cases

    # Create variables for the fields
    for field in fields:
        rootgrp.createVariable(
            field,
            "f8",
            (
                "time",
                "z",
                "x",
                "y",
            ),
        )

    # Process each case
    for i, casei in enumerate(cases):
        print(f"{i}/{len(cases)}", end="\r")
        case_name_id = str(int(1000000 + casei + 1))[1:]
        result_dir = case_name + "_" + case_name_id
        ens = cs_pp.load_ensight_data(
            case_dir + sep + "RESU" + sep + result_dir + sep + "postprocessing",
            "RESULTS_FLUID_DOMAIN.case",
        )
        data = cs_pp.read_ensight_data_at_time(ens, inst=0)

        if i == 0:
            print(case_dir + sep + "RESU" + sep + result_dir + sep + "postprocessing")
            cell_centers = cs_pp.get_cell_centers(data)
            num_points = cell_centers.GetNumberOfPoints()
            x = np.empty(num_points)
            y = np.empty(num_points)
            z = np.empty(num_points)
            for n in range(num_points):
                x[n], y[n], z[n] = cell_centers.GetPoint(n)

            indices_double = find_closest(z, zplot[0])
            indices = indices_double[::2]  # Take 1 indice sur 2
            x0 = x[indices]
            y0 = y[indices]
            z0 = z[indices]
            npoints = len(indices)

        for j, zj in enumerate(zplot):
            indices_double = find_closest(z, zj)
            indices = indices_double[::2]
            if "wind_speed" in fields or "wind_direction" in fields:
                velocity = np.empty([len(x0), 3])
                field = data.GetCellData().GetArray("Velocity")
                for m, index in enumerate(indices):
                    velocity[m, :] = field.GetTuple3(index)

                if "wind_speed" in fields:
                    speed = np.sqrt(velocity[:, 0] ** 2 + velocity[:, 1] ** 2)
                    speed_interp = griddata(
                        (x0, y0), speed, (x_grid, y_grid), method="linear"
                    )
                    rootgrp.variables["wind_speed"][i, j, :, :] = speed_interp.T

                if "wind_direction" in fields:
                    direction = (
                        np.arctan2(velocity[:, 1], velocity[:, 0]) * 180 / np.pi + 270
                    )
                    direction_interp = griddata(
                        (x0, y0), direction, (x_grid, y_grid), method="linear"
                    )
                    rootgrp.variables["wind_direction"][i, j, :, :] = direction_interp.T

            if "pressure" in fields:
                pressure = np.empty(len(x0))
                field = data.GetCellData().GetArray("total_pressure")
                for m, index in enumerate(indices):
                    pressure[m] = field.GetTuple(index)
                pressure_interp = griddata(
                    (x0, y0), pressure, (x_grid, y_grid), method="linear"
                )
                rootgrp.variables["pressure"][i, j, :, :] = pressure_interp.T

            if "tke" in fields:
                tke = np.empty(len(x0))
                field = data.GetCellData().GetArray("k")
                for m, index in enumerate(indices):
                    tke[m] = field.GetTuple(index)
                tke_interp = griddata((x0, y0), tke, (x_grid, y_grid), method="linear")
                rootgrp.variables["tke"][i, j, :, :] = tke_interp.T

    # Close the NetCDF file
    rootgrp.close()
