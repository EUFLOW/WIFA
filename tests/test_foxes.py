import os
import numpy as np
from pathlib import Path
from shutil import rmtree

from windIO import __path__ as wiop
from windIO import validate as validate_yaml

from wifa.foxes_api import run_foxes

test_path = Path(os.path.dirname(__file__))
windIO_path = Path(wiop[0])


def _run_foxes(wes_dir):
    assert wes_dir.is_dir(), f"{wes_dir} is not a directory"

    for yaml_input in wes_dir.glob("system.yaml"):
        print("\nRUNNING FOXES ON", yaml_input, "\n")
        validate_yaml(yaml_input, Path("plant/wind_energy_system"))
        output_dir_name = Path("output_test_foxes")
        output_dir_name.mkdir(parents=True, exist_ok=True)
        run_foxes(yaml_input, output_dir=output_dir_name)
        rmtree(output_dir_name)


def test_foxes_KUL():
    wes_dir = test_path / "../examples/cases/KUL_LES/wind_energy_system/"
    _run_foxes(wes_dir)


def test_foxes_4wts():
    wes_dir = test_path / "../examples/cases/windio_4turbines/wind_energy_system/"
    _run_foxes(wes_dir)


def test_foxes_abl():
    wes_dir = test_path / "../examples/cases/windio_4turbines_ABL/wind_energy_system/"
    _run_foxes(wes_dir)


def test_foxes_abl_stable():
    wes_dir = (
        test_path / "../examples/cases/windio_4turbines_ABL_stable/wind_energy_system/"
    )
    _run_foxes(wes_dir)


def test_foxes_profiles():
    wes_dir = (
        test_path
        / "../examples/cases/windio_4turbines_profiles_stable/wind_energy_system/"
    )
    _run_foxes(wes_dir)


def test_foxes_heterogeneous_wind_rose_at_turbines():
    wes_dir = (
        test_path
        / "../examples/cases/heterogeneous_wind_rose_at_turbines/wind_energy_system/"
    )
    _run_foxes(wes_dir)


def test_foxes_heterogeneous_wind_rose_map():
    wes_dir = (
        test_path / "../examples/cases/heterogeneous_wind_rose_map/wind_energy_system/"
    )
    _run_foxes(wes_dir)


def test_foxes_simple_wind_rose():
    wes_dir = test_path / "../examples/cases/simple_wind_rose/wind_energy_system/"
    _run_foxes(wes_dir)


def test_foxes_timeseries_with_operating_flag():
    wes_dir = (
        test_path
        / "../examples/cases/timeseries_with_operating_flag/wind_energy_system/"
    )
    _run_foxes(wes_dir)

def test_timeseries_per_turbine_with_density(tmp_path=Path(".")):
    from conftest import make_timeseries_per_turbine_system_dict
    import foxes.variables as FV

    # Run with density
    system_dict = make_timeseries_per_turbine_system_dict("foxes")
    output_dir = tmp_path / "output_foxes_ts"
    farm_results = run_foxes(system_dict, verbosity=0, output_dir=str(output_dir))[0]
    farmP_with = farm_results[FV.P].sum()
    print("Farm power with density:", farmP_with)
    assert np.isfinite(farmP_with) and farmP_with > 0

    # Run without density — same config but density removed
    system_dict_no = make_timeseries_per_turbine_system_dict("foxes")
    del system_dict_no["site"]["energy_resource"]["wind_resource"]["density"]
    output_dir_no = tmp_path / "output_foxes_ts_no_density"
    farm_results_no = run_foxes(system_dict_no, verbosity=0, output_dir=str(output_dir_no))[0]
    farmP_without = farm_results_no[FV.P].sum()
    print("Farm power without density:", farmP_without)

    rmtree(output_dir)
    rmtree(output_dir_no)

    # Density correction should change AEP (test data varies around 1.225)
    assert farmP_with != farmP_without

if __name__ == "__main__":
    #test_foxes_KUL()
    #test_foxes_4wts()
    #test_foxes_abl()
    #test_foxes_abl_stable()
    #test_foxes_profiles()
    #test_foxes_heterogeneous_wind_rose_at_turbines()
    #test_foxes_heterogeneous_wind_rose_map()
    #test_foxes_simple_wind_rose()
    test_timeseries_per_turbine_with_density()
