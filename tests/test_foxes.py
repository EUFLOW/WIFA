import os
from pathlib import Path
from shutil import rmtree

import numpy as np
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


def test_timeseries_per_turbine_with_density(output_dir):
    import foxes.variables as FV
    from conftest import make_timeseries_per_turbine_system_dict

    # Run with density
    system_dict = make_timeseries_per_turbine_system_dict("foxes")
    out_with = output_dir / "output_foxes_ts"
    farm_results = run_foxes(system_dict, verbosity=0, output_dir=str(out_with))[0]
    farmP_with = farm_results[FV.P].sum()
    assert np.isfinite(farmP_with) and farmP_with > 0

    # Run without density — same config but density removed
    system_dict_no = make_timeseries_per_turbine_system_dict("foxes")
    del system_dict_no["site"]["energy_resource"]["wind_resource"]["density"]
    out_without = output_dir / "output_foxes_ts_no_density"
    farm_results_no = run_foxes(
        system_dict_no, verbosity=0, output_dir=str(out_without)
    )[0]
    farmP_without = farm_results_no[FV.P].sum()

    # Density correction should change AEP (test data varies around 1.225)
    assert farmP_with != farmP_without
