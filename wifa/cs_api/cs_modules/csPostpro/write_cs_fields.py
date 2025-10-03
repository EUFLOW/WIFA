import sys
from os import path, sep

import numpy as np
from netCDF4 import Dataset

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


postpro_dir = sys.argv[1]
file_name = sys.argv[2]
case_name = sys.argv[5]
case_dir = sys.argv[6]

postprocess_heights_file = open(postpro_dir + "/postprocessing_heights.csv", "r")
postprocess_cases_file = open(postpro_dir + "/postprocessing_cases.csv", "r")
postprocess_fields_file = open(postpro_dir + "/postprocessing_fields.csv", "r")

cases = np.atleast_1d(np.genfromtxt(postprocess_cases_file, delimiter=","))
zplot = np.array(np.genfromtxt(postprocess_heights_file, delimiter=","), ndmin=1)
fields = np.genfromtxt(postprocess_fields_file, delimiter=",", dtype=str)

###############CS FIELD DATA###################
# File creation
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
    if i == 0:
        x, y, z, triang = get_plane_info(ens, (0, 0, zplot[0]), (0, 0, 1))
        npoints = triang.triangles.shape[0]
        xtri = np.zeros((npoints))
        ytri = np.zeros((npoints))
        for k in range(npoints):
            i1, i2, i3 = triang.triangles[k]
            xtri[k] = (x[i1] + x[i2] + x[i3]) / 3.0
            ytri[k] = (y[i1] + y[i2] + y[i3]) / 3.0

        #
        points = rootgrp.createDimension("points", npoints)
        nc_cases = rootgrp.createDimension("time", len(cases))
        altitudes = rootgrp.createDimension("altitudes", len(zplot))
        rootgrp.createVariable("x", "f8", ("points",))
        rootgrp.createVariable("y", "f8", ("points",))
        rootgrp.createVariable("z", "f8", ("altitudes",))
        rootgrp.createVariable("time", "f8", ("time",))
        rootgrp.variables["x"][:] = xtri
        rootgrp.variables["y"][:] = ytri
        rootgrp.variables["z"][:] = zplot
        rootgrp.variables["time"][:] = cases
        for field in fields:
            rootgrp.createVariable(
                field,
                "f8",
                (
                    "points",
                    "altitudes",
                    "time",
                ),
            )
    for j, zj in enumerate(zplot):
        zplane_center = (0, 0, zj)

        if "wind_speed" or ("wind_direction" in fields):
            velocity = get_output_at_plane_and_time(
                ens, "Velocity", zplane_center, (0, 0, 1), 0
            )
            if "wind_speed" in fields:
                speed = np.sqrt(pow(velocity[:, 0], 2.0) + pow(velocity[:, 1], 2.0))
                rootgrp.variables["wind_speed"][:, j, i] = speed
            if "wind_direction" in fields:
                direction = (
                    np.arctan(velocity[:, 1] / velocity[:, 0]) * 360 / (2 * np.pi) + 270
                )
                rootgrp.variables["wind_direction"][:, j, i] = direction
        if "pressure" in fields:
            pressure = get_output_at_plane_and_time(
                ens, "total_pressure", zplane_center, (0, 0, 1), 0
            )
            rootgrp.variables["pressure"][:, j, i] = pressure
        if "tke" in fields:
            tke = get_output_at_plane_and_time(ens, "k", zplane_center, (0, 0, 1), 0)
            rootgrp.variables["tke"][:, j, i] = tke

rootgrp.close()
