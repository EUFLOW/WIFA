import os
from pathlib import Path

from windIO import __path__ as wiop
from windIO import validate as validate_yaml

from wifa.cs_api.cs_modules.csLaunch.cs_run_function import run_code_saturne

test_path = Path(os.path.dirname(__file__))
windIO_path = Path(wiop[0])


def _run_cs(wes_dir, output_dir):
    i = 1
    for yaml_input in wes_dir.glob("system*"):
        print("\nRUNNING CODE_SATURNE ON", yaml_input, "\n")
        validate_yaml(yaml_input, Path("plant/wind_energy_system"))
        run_code_saturne(
            yaml_input,
            test_mode=False,
            output_dir="output_test_" + output_dir + "_" + str(i),
        )
        i += 1


def test_cs_KUL():
    wes_dir = test_path / "../examples/cases/KUL_LES/wind_energy_system/"
    _run_cs(wes_dir, "KUL")


def test_cs_4wts():
    wes_dir = test_path / "../examples/cases/windio_4turbines/wind_energy_system/"
    _run_cs(wes_dir, "4wts")


def test_cs_abl():
    wes_dir = test_path / "../examples/cases/windio_4turbines_ABL/wind_energy_system/"
    _run_cs(wes_dir, "abl")


def test_cs_abl_stable():
    wes_dir = (
        test_path / "../examples/cases/windio_4turbines_ABL_stable/wind_energy_system/"
    )
    _run_cs(wes_dir, "abl_stable")


def test_cs_profiles():
    wes_dir = (
        test_path
        / "../examples/cases/windio_4turbines_profiles_stable/wind_energy_system/"
    )
    _run_cs(wes_dir, "profiles")


if __name__ == "__main__":
    test_cs_KUL()
    test_cs_4wts()
    test_cs_abl()
    test_cs_abl_stable()
    test_cs_profiles()
