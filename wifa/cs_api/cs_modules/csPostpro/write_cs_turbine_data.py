from netCDF4 import Dataset
import sys
from os import path, sep
import numpy as np
import wifa.cs_api.cs_modules.csPostpro.cs_postprocess_utils as cs_pp


postpro_dir = sys.argv[1]
turbine_file_name = sys.argv[2]
turbine_number = int(sys.argv[3])
case_name = sys.argv[4]
case_dir = sys.argv[5]
ntmax_str = sys.argv[6]

if not path.exists(postpro_dir):
    sys.mkdir(postpro_dir)

postprocess_cases_file = open(postpro_dir + "/postprocessing_cases.csv", "r")

cases = np.atleast_1d(np.genfromtxt(postprocess_cases_file, delimiter=","))

################TURBINE DATA###################
# File creation
rootgrp = Dataset(postpro_dir + sep + turbine_file_name, "w", format="NETCDF4")
# Dimensions
turbines = rootgrp.createDimension("turbine", turbine_number)
time = rootgrp.createDimension("time", len(cases))
u_file = rootgrp.createVariable(
    "rotor_effective_velocity",
    "f8",
    (
        "turbine",
        "time",
    ),
)
dir_file = rootgrp.createVariable(
    "wind_direction",
    "f8",
    (
        "turbine",
        "time",
    ),
)
powercp_file = rootgrp.createVariable(
    "power",
    "f8",
    (
        "turbine",
        "time",
    ),
)

thrust_file = rootgrp.createVariable(
    "thrust",
    "f8",
    (
        "turbine",
        "time",
    ),
)

turbines_file = rootgrp.createVariable(
    "turbine",
    "f8",
    ("turbine",
    ),
)
time_file = rootgrp.createVariable(
    "time",
    "f8",
    ( "time",
    ),
)
turbines_file[:]=np.arange(turbine_number)
time_file[:]=cases
for j, casei in enumerate(cases):
    print(str(j) + "/" + str(len(cases)), end="\r")
    # TODO: string formatting for id "%05d"
    case_name_id = str(int(1000000 + casei + 1))[1:]
    result_dir = case_name + "_" + case_name_id
    power_file = (
        case_dir
        + sep
        + "RESU"
        + sep
        + result_dir
        + sep
        + "power_iter"
        + ntmax_str
        + ".csv"
    )
    total_power = []
    x_coords = []
    y_coords = []
    z_hub = []
    diameters = []
    ux = []
    uy = []
    uz = []
    u = []
    u_hub = []
    dir_table = []
    ctstar = []
    cpstar = []
    thrust = []
    power_table_cpstar = []

    # Get turbine and power output info
    with open(power_file, "r") as file:
        # Read the first line
        first_line = file.readline()
        total_power.append(float(first_line.split(" ")[-1]))
        second_line = file.readlines()[0]
        var_name = second_line.replace(" ", "").replace("\n", "").split(",")
    power_file_table = np.atleast_2d(np.genfromtxt(power_file, delimiter=",", skip_header=2))
    x_coords.append(
        power_file_table[:, var_name.index("xhub")]
        - np.mean(power_file_table[:, var_name.index("xhub")])
    )
    y_coords.append(
        power_file_table[:, var_name.index("yhub")]
        - np.mean(power_file_table[:, var_name.index("yhub")])
    )
    z_hub.append(power_file_table[:, var_name.index("zhub")])
    diameters.append(power_file_table[:, var_name.index("turbine_diameter")])
    ux.append(power_file_table[:, var_name.index("ux")])
    uy.append(power_file_table[:, var_name.index("uy")])
    uz.append(power_file_table[:, var_name.index("uz")])
    u.append(power_file_table[:, var_name.index("u")])
    u_hub.append(power_file_table[:, var_name.index("u_hub")])
    dir_table.append(power_file_table[:, var_name.index("dir")])
    ctstar.append(power_file_table[:, var_name.index("ct*")])
    cpstar.append(power_file_table[:, var_name.index("cp*")])
    thrust.append(power_file_table[:, var_name.index("thrust")])
    power_table_cpstar.append(power_file_table[:, var_name.index("power_cpstar")])

    u_file[:, j] = u[0]
    dir_file[:, j] = dir_table[0]
    powercp_file[:, j] = power_table_cpstar[0]
    thrust_file[:, j] = thrust[0]
rootgrp.close()
