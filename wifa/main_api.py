import argparse
import os
import sys

import windIO
import yaml
from windIO import load_yaml
from windIO import validate as validate_yaml

from .cs_api.cs_modules.csLaunch.cs_run_function import run_code_saturne
from .floris_api import run_floris
from .foxes_api import run_foxes
from .pywake_api import run_pywake
from .wayve_api import run_wayve

sys.path.append(windIO.__path__[0])


def run_api(yaml_input):
    # validate input
    validate_yaml(yaml_input, windIO.__path__[0] + "/plant/wind_energy_system.yaml")

    # get number of turbines
    if isinstance(yaml_input, dict):
        yaml_dat = yaml_input
    else:
        yaml_dat = load_yaml(yaml_input)

    model_name = yaml_dat["attributes"]["flow_model"]["name"]

    if model_name.lower() == "pywake":
        pywake_aep = run_pywake(yaml_input)

    elif model_name.lower() == "foxes":
        foxes_aep = run_foxes(yaml_input)

    elif model_name.lower() == "floris":
        floris_aep = run_floris(yaml_input)

    elif model_name.lower() == "wayve":
        # Output directory
        # yaml_input_no_ext = os.path.splitext(yaml_input)[0]  # Remove the file extension
        # output_dir_name = 'output_wayve' + yaml_input_no_ext.replace(os.sep, '_')  # Replace directory separators
        output_dir_name = yaml_dat["attributes"]["model_outputs_specification"][
            "output_folder"
        ]
        if not os.path.exists(output_dir_name):
            os.makedirs(output_dir_name)

        run_wayve(yaml_input, output_dir_name)

    elif model_name.lower() == "codesaturne":
        run_code_saturne(yaml_input, test_mode=True)

    else:
        print("Invalid Model")


def run():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_yaml", help="The input yaml file")
    args = parser.parse_args()

    run_api(args.input_yaml)


if __name__ == "__main__":
    run()
