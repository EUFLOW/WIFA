import os as os
import shutil
from datetime import datetime
from os import mkdir, sep, walk

import numpy as np

from wifa.cs_api import (
    cs_api_path,
    cs_env_command,
    cs_exe_path,
    python_scripts_env_command,
    python_scripts_exe,
    salome_env_command,
    salome_exe_path,
)
from wifa.cs_api.cs_modules.csLaunch.cs_launch_modules import *


def initialize_cs_case_from_windio(windio_input, output_dir):
    """Runner to code_saturne for the FLOW api

    Parameters:
    windio_input (str): main windio file path

    """

    # TODO: paths from windio or a config file in the api during install
    if output_dir is None:
        cs_run_folder = "cs_run_" + datetime.now().strftime("%Y%m%d_%H-%M-%S")
    else:
        cs_run_folder = output_dir
    #
    windfarm_study = CS_study(
        cs_run_folder=cs_run_folder,
        case_dir="Farm",
        cs_path=cs_exe_path,
        cs_api_path=cs_api_path,
        salome_path=salome_exe_path,
        python_env_command=python_scripts_env_command,
        salome_env_command=salome_env_command,
        cs_env_command=cs_env_command,
        python_exe=python_scripts_exe,
    )

    # Example 1 : get data from windio files
    windfarm_study.set_windio(windio_input)
    windfarm_study.get_windio_data()
    return windfarm_study


def run_cs_windfarm_study(windfarm_study, test_mode=False, postprocess_only=False):
    """Run the wind farm study

    Parameters:
    windfarm_study (class): main CS_study object from params

    """

    if postprocess_only:
        if not os.path.exists(windfarm_study.cs_run_folder):
            raise ValueError("Result directory does not exist.")
    else:
        # Creates the temporary folders for the run
        windfarm_study.set_run_env()

    # TODO: get case_name from windio
    windfarm_study.case_name = "wf"
    ntmax = 5000
    WTntprint = 100
    precntmax = 1000000

    if test_mode:
        windfarm_study.case_name = "test"
        windfarm_study.run_wall_time_hours = 1.0 / 6.0  # 10minutes for testing
        ntmax = 5  # max number of timesteps for cs run
        WTntprint = 1  # period of turbine output printout for cs run
        precntmax = 5000  # max precursor number of timesteps. Multiples of 5000 only.
        windfarm_study.run_wall_time_hours = 1.0 / 6.0  # 10minutes
        windfarm_study.prec_wall_time_hours = 1.0 / 6.0  # 10minutes

    # =============USER SIMULATION INFO ===============
    # Default that will be changed below depending on windio case
    farm_notebook_parameters = {
        "meteo_profile": 1,
        #
        "teta": 270,
        #
        "st_method": 0,  # 0 : homogeneous AD
        "isol": 0,  # 0 if full farm, i>0 if turbine i in isolation
        #
        "energy": 1,  # 1 to solve energy equation, 0 for constant density
        #
        "Coriolis": 1,  # 1 if coriolis, 0 if not
        #
        "damping": 0,  # activate a dumping layer for gravity waves
        #
        "control": 0,  # activate control of wind turbines with local wind dir
        #
        "ntmax": ntmax,  # max number of timesteps. Default is 5000.
        "WTntprint": WTntprint,  # period of turbine output printout. Default is 100.
    }
    prec_notebook_parameters = {
        "z0": 0.0001,
        "zref": windfarm_study.farm.hub_heights[0],  # TODO : for multiple turbines
        "ureff": 10.0,
        #
        "t0": 293.15,
        "Lmoinv": 0.0,
        #
        "lat": 55,  # latitude for coriolis
        #
        "precntmax": precntmax,  # max precursor number of timesteps. Default is 1000000.
    }

    # cs files can be written before loop
    # as long as layout and turbines do not change
    windfarm_study.write_cs_input_files()

    if windfarm_study.inflow.capping_inversion == True:
        damping_layer = True
        windfarm_study.mesh.damping_length = 10000.0
        windfarm_study.mesh.domain_height = 25000.0
        turbine_control = True
    else:
        windfarm_study.mesh.domain_height = 800.0
        damping_layer = False
        turbine_control = False  # TODO : depends on Coriolis
    if windfarm_study.inflow.run_precursor:
        precursor = True
    else:
        precursor = False

    if not (windfarm_study.inflow.coriolis):
        farm_notebook_parameters["Coriolis"] = 0
    #
    remesh = windfarm_study.mesh.remesh
    mesh_file_name = windfarm_study.mesh.mesh_file_name

    # Run
    if postprocess_only:
        windfarm_study.set_notebook_param_from_dictionary(
            farm_notebook_parameters, prec_notebook_parameters
        )
        launch_file_name = "postprocess.sh"
        windfarm_study.postprocess(
            standalone=True, launch_file_name="postprocess.sh", log_folder="logs"
        )
        os.system(
            "cd " + windfarm_study.cs_run_folder + " ; sbatch " + launch_file_name
        )
    else:
        launch_file_name = "launch_farm.sh"
        # ===================Time loop====================
        for j in windfarm_study.inflow.run_times:
            #
            case_name_id = str(int(1000000 + j + 1))[1:]
            windfarm_study.inflow.time_iter = j
            # =============Set code_saturne parameters================
            # modify parameter values
            prec_notebook_parameters["z0"] = np.round(
                windfarm_study.inflow.roughness_height[j], 6
            )
            prec_notebook_parameters["lat"] = windfarm_study.inflow.latitude[j]
            #
            if windfarm_study.inflow.data_type == "timeseries_hub":
                prec_notebook_parameters["ureff"] = np.round(
                    windfarm_study.inflow.wind_velocity[j], 2
                )
                prec_notebook_parameters[
                    "t0"
                ] = 293.15  # 20 deg Celcius, arbitrary #TODO: user choice?
                farm_notebook_parameters["teta"] = np.round(
                    windfarm_study.inflow.wind_dir[j], 2
                )
                windfarm_study.wind_origin = np.round(
                    windfarm_study.inflow.wind_dir[j], 2
                )
                if windfarm_study.inflow.capping_inversion:
                    farm_notebook_parameters["meteo_profile"] = 1
                    #
                    meteo_file_name = "meteo_files" + sep + "meteo_file_" + case_name_id
                    precursor_meteo_file_name = (
                        "precursor_meteo_files" + sep + "meteo_file_" + case_name_id
                    )
                    #
                    precursor = True
                    #
                    if np.abs(windfarm_study.inflow.LMO_values[j]) > 1000:
                        windfarm_study.inflow.pottemp[
                            :, j
                        ] = windfarm_study.generate_temp_CNBL(j)
                        windfarm_study.inflow.u[:, j] = (
                            np.ones(len(windfarm_study.inflow.heights))
                            * windfarm_study.inflow.wind_velocity[j]
                        )
                        windfarm_study.inflow.v[:, j] = np.zeros(
                            len(windfarm_study.inflow.heights)
                        )
                        windfarm_study.inflow.tke[:, j] = (
                            np.ones(len(windfarm_study.inflow.heights)) * 0.1
                        )
                        windfarm_study.inflow.epsilon[:, j] = (
                            np.ones(len(windfarm_study.inflow.heights)) * 0.003
                        )
                        windfarm_study.write_cs_meteo_file(
                            precursor_meteo_file_name, j, precursor
                        )
                        #
                        prec_notebook_parameters["Lmoinv"] = 0.0  # neutral
                        #
                        #
                        prec_notebook_parameters["ustar"] = 0
                        prec_notebook_parameters["tstar"] = 0
                        prec_notebook_parameters["zi"] = 0
                    elif windfarm_study.inflow.LMO_values[j] > 0:
                        (
                            pottemp,
                            u,
                            v,
                            tke,
                            epsilon,
                            ustar,
                            tstar,
                            zi,
                        ) = windfarm_study.generate_prof_stable(j)
                        windfarm_study.inflow.pottemp[:, j] = pottemp
                        windfarm_study.inflow.u[:, j] = u
                        windfarm_study.inflow.v[:, j] = v
                        windfarm_study.inflow.tke[:, j] = tke
                        windfarm_study.inflow.epsilon[:, j] = epsilon
                        windfarm_study.write_cs_meteo_file(
                            precursor_meteo_file_name, j, precursor
                        )
                        #
                        prec_notebook_parameters["Lmoinv"] = (
                            1.0 / windfarm_study.inflow.LMO_values[j]
                        )
                        if not test_mode:
                            prec_notebook_parameters["precntmax"] = 20000
                        #
                        prec_notebook_parameters["ustar"] = ustar
                        prec_notebook_parameters["tstar"] = tstar
                        prec_notebook_parameters["zi"] = zi
                    else:
                        raise ValueError(
                            "Ongoing devs - can not handle unstable stratification with capping inversion"
                        )
                else:
                    farm_notebook_parameters["meteo_profile"] = 2
                    farm_notebook_parameters["energy"] = 0
                    meteo_file_name = ""
                    precursor_meteo_file_name = ""
                    if np.abs(windfarm_study.inflow.LMO_values[j]) > 1000:
                        prec_notebook_parameters["Lmoinv"] = 0.0  # neutral
                    else:
                        prec_notebook_parameters["Lmoinv"] = np.round(
                            1.0 / windfarm_study.inflow.LMO_values[j], 3
                        )
            #
            elif windfarm_study.inflow.data_type == "timeseries_profile":
                farm_notebook_parameters["meteo_profile"] = 1
                farm_notebook_parameters["energy"] = 0
                #
                meteo_file_name = "meteo_files" + sep + "meteo_file_" + case_name_id
                if precursor:
                    precursor_meteo_file_name = (
                        "precursor_meteo_files" + sep + "meteo_file_" + case_name_id
                    )
                    windfarm_study.write_cs_meteo_file(precursor_meteo_file_name, j)

                    windfarm_study.wind_origin = (
                        windfarm_study.get_wind_dir_from_meteo_file(
                            precursor_meteo_file_name,
                            zmeteo=windfarm_study.farm.hub_heights[0],
                        )
                    )
                else:
                    precursor_meteo_file_name = None
                    windfarm_study.write_cs_meteo_file(meteo_file_name, j)

                    windfarm_study.wind_origin = (
                        windfarm_study.get_wind_dir_from_meteo_file(
                            meteo_file_name, zmeteo=windfarm_study.farm.hub_heights[0]
                        )
                    )
                farm_notebook_parameters["teta"] = windfarm_study.wind_origin

            else:
                raise ValueError(
                    "Can not handle this inflow type : "
                    + windfarm_study.inflow.data_type
                )

            if damping_layer:
                windfarm_study.inflow.damping_lapse = (
                    windfarm_study.get_lapse_rate_from_meteo_file(meteo_file_name)
                )
                farm_notebook_parameters["damping"] = 1
                farm_notebook_parameters["gamma"] = windfarm_study.inflow.damping_lapse
                farm_notebook_parameters["nura"] = 3.0
                farm_notebook_parameters["Lra"] = windfarm_study.mesh.damping_length
                farm_notebook_parameters["Sra"] = 2.0
                farm_notebook_parameters[
                    "start_rad"
                ] = 20000  # default, modified with domain size in run_case
            if turbine_control:
                farm_notebook_parameters["control"] = 1

            # =================================================
            windfarm_study.set_notebook_param_from_dictionary(
                farm_notebook_parameters, prec_notebook_parameters
            )
            #
            # ============RUN SIMULATION AND LOG===============
            windfarm_study.set_result_dir(windfarm_study.case_name + "_" + case_name_id)
            # mesh_file_name = "mesh_"+case_name_id+".med"
            job_name = windfarm_study.case_name + "_" + str(j + 1)
            #
            if j == windfarm_study.inflow.run_times[0]:
                first_case = True
            else:
                first_case = False
            windfarm_study.run_case(
                first_case=first_case,
                remesh=remesh,
                turbine_control=turbine_control,
                damping_layer=damping_layer,
                mesh_file_name=mesh_file_name,
                launch_file_name=launch_file_name,
                job_name=job_name,
                log_folder="logs",
                meteo_file_name=meteo_file_name,
                precursor=precursor,
                precursor_meteo_file_name=precursor_meteo_file_name,
            )
            #
            # TODO : add verbose option
            # print('Wrote and launched "'+launch_file_name+'". Waiting for job to finish')
            # if(remesh):
            #    print("After meshing is finished, the mesh will be stored as "+"MESH"+sep+mesh_file_name)
            # else:
            #    print("Used existing mesh is "+"MESH"+sep+mesh_file_name)
            # print("After cs simulation is finished, power output will be stored in "+windfarm_study.case_dir+sep+"RESU"+sep+windfarm_study.result_dir+sep+"power.txt file")

        windfarm_study.postprocess(launch_file_name=launch_file_name, log_folder="logs")
        # sys.exit()
        os.system(
            "cd " + windfarm_study.cs_run_folder + " ; sbatch " + launch_file_name
        )
        # =================================================


def run_code_saturne(
    windio_input, test_mode=False, output_dir=None, postprocess_only=False
):
    """Runner to code_saturne for the FLOW api

    Parameters:
    windio_input (str): main windio file path

    """
    windfarm_study = initialize_cs_case_from_windio(windio_input, output_dir)
    run_cs_windfarm_study(
        windfarm_study, test_mode=test_mode, postprocess_only=postprocess_only
    )


def validate_yaml_code_saturne(windio_input):
    """Validation of the windio content for code_saturne
    - Temperature profile
    - Coherence between asked stability and LMO value
    - Warning for no ABL and Default choice
    - Etc.

    Parameters:
    windio_input (str): main windio file path
    """
    # TODO : check temperature profile

    return 0


def run():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_yaml", help="The input yaml file")
    args = parser.parse_args()

    run_code_saturne(
        args.input_yaml, test_mode=False, output_dir=None, postprocess_only=False
    )


if __name__ == "__main__":
    run()
