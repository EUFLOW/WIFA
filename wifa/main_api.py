import argparse
import os
import sys

import windIO
from windIO import validate as validate_yaml

from .cs_api.cs_modules.csLaunch.cs_run_function import run_code_saturne
from .floris_api import run_floris
from .foxes_api import run_foxes
from .pywake_api import run_pywake
from .wayve_api import run_wayve

sys.path.append(windIO.__path__[0])


def run_api(yaml_input):
    if isinstance(yaml_input, dict):
        yaml_dat = yaml_input
    else:
        yaml_dat = validate_yaml(yaml_input, "plant/wind_energy_system")

    model_name = yaml_dat["attributes"]["flow_model"]["name"]

    if model_name.lower() == "pywake":
        run_pywake(yaml_dat)

    elif model_name.lower() == "foxes":
        run_foxes(yaml_dat)

    elif model_name.lower() == "floris":
        run_floris(yaml_dat)

    elif model_name.lower() == "wayve":
        output_dir_name = yaml_dat["attributes"]["model_outputs_specification"][
            "output_folder"
        ]
        if not os.path.exists(output_dir_name):
            os.makedirs(output_dir_name)

        run_wayve(yaml_dat, output_dir_name)

    elif model_name.lower() == "codesaturne":
        run_code_saturne(yaml_dat, test_mode=True)

    else:
        print("Invalid Model")


def run():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_yaml", help="The input yaml file")
    args = parser.parse_args()

    run_api(args.input_yaml)


if __name__ == "__main__":
    run()
