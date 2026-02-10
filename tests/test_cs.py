import os
from pathlib import Path

from windIO import __path__ as wiop
from windIO import validate as validate_yaml

from wifa.cs_api.cs_modules.csLaunch.cs_run_function import run_code_saturne

test_path = Path(os.path.dirname(__file__))
windIO_path = Path(wiop[0])


def _run_cs(wes_dir, output_dir):
    """Run Code Saturne on all system* files in the given directory."""
    i = 1
    for yaml_input in wes_dir.glob("system*"):
        print("\nRUNNING CODE_SATURNE ON", yaml_input, "\n")
        validate_yaml(yaml_input, Path("plant/wind_energy_system"))
        # Pass subdirectory path - run_code_saturne will create it
        sub_output_dir = output_dir / f"run_{i}"
        run_code_saturne(
            yaml_input,
            test_mode=False,
            output_dir=str(sub_output_dir),
        )
        i += 1


def test_cs_KUL(output_dir):
    wes_dir = test_path / "../examples/cases/KUL_LES/wind_energy_system/"
    _run_cs(wes_dir, output_dir)


def test_cs_4wts(output_dir):
    wes_dir = test_path / "../examples/cases/windio_4turbines/wind_energy_system/"
    _run_cs(wes_dir, output_dir)


def test_cs_abl(output_dir):
    wes_dir = test_path / "../examples/cases/windio_4turbines_ABL/wind_energy_system/"
    _run_cs(wes_dir, output_dir)


def test_cs_abl_stable(output_dir):
    wes_dir = (
        test_path / "../examples/cases/windio_4turbines_ABL_stable/wind_energy_system/"
    )
    _run_cs(wes_dir, output_dir)


def test_cs_profiles(output_dir):
    wes_dir = (
        test_path
        / "../examples/cases/windio_4turbines_profiles_stable/wind_energy_system/"
    )
    _run_cs(wes_dir, output_dir)
