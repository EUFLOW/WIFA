import os
from pathlib import Path
from shutil import rmtree

from windIO import __path__ as wiop
from windIO import validate as validate_yaml

from wifa.foxes_api import run_foxes

test_path = Path(os.path.dirname(__file__))
windIO_path = Path(wiop[0])


def _run_foxes(wes_dir, output_dir):
    """Run FOXES on all system.yaml files in the given directory."""
    assert wes_dir.is_dir(), f"{wes_dir} is not a directory"

    for yaml_input in wes_dir.glob("system.yaml"):
        print("\nRUNNING FOXES ON", yaml_input, "\n")
        validate_yaml(yaml_input, Path("plant/wind_energy_system"))
        run_foxes(yaml_input, output_dir=output_dir)


def test_foxes_KUL(output_dir):
    wes_dir = test_path / "../examples/cases/KUL_LES/wind_energy_system/"
    _run_foxes(wes_dir, output_dir)


def test_foxes_4wts(output_dir):
    wes_dir = test_path / "../examples/cases/windio_4turbines/wind_energy_system/"
    _run_foxes(wes_dir, output_dir)


def test_foxes_abl(output_dir):
    wes_dir = test_path / "../examples/cases/windio_4turbines_ABL/wind_energy_system/"
    _run_foxes(wes_dir, output_dir)


def test_foxes_abl_stable(output_dir):
    wes_dir = (
        test_path / "../examples/cases/windio_4turbines_ABL_stable/wind_energy_system/"
    )
    _run_foxes(wes_dir, output_dir)


def test_foxes_profiles(output_dir):
    wes_dir = (
        test_path
        / "../examples/cases/windio_4turbines_profiles_stable/wind_energy_system/"
    )
    _run_foxes(wes_dir, output_dir)


def test_foxes_heterogeneous_wind_rose_at_turbines(output_dir):
    wes_dir = (
        test_path
        / "../examples/cases/heterogeneous_wind_rose_at_turbines/wind_energy_system/"
    )
    _run_foxes(wes_dir, output_dir)


def test_foxes_heterogeneous_wind_rose_map(output_dir):
    wes_dir = (
        test_path / "../examples/cases/heterogeneous_wind_rose_map/wind_energy_system/"
    )
    _run_foxes(wes_dir, output_dir)


def test_foxes_simple_wind_rose(output_dir):
    wes_dir = test_path / "../examples/cases/simple_wind_rose/wind_energy_system/"
    _run_foxes(wes_dir, output_dir)


def test_foxes_timeseries_with_operating_flag(output_dir):
    wes_dir = (
        test_path
        / "../examples/cases/timeseries_with_operating_flag/wind_energy_system/"
    )
    _run_foxes(wes_dir, output_dir)
