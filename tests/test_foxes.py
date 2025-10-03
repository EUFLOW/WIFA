import os
from pathlib import Path

from foxes import Engine, reset_engine
from windIO import __path__ as wiop
from windIO.utils.yml_utils import validate_yaml

from wifa.foxes_api import run_foxes

test_path = Path(os.path.dirname(__file__))
windIO_path = Path(wiop[0])

engine = None


def _run_foxes(wes_dir):
    assert wes_dir.is_dir(), f"{wes_dir} is not a directory"

    global engine
    if engine is None:
        engine = Engine.new("process", verbosity=0)
        engine.initialize()
        print("SETTING ENGINE:", engine)

    for yaml_input in wes_dir.glob("system*.yaml"):
        if "_noXYgrid" not in str(yaml_input):
            print("\nRUNNING FOXES ON", yaml_input, "\n")
            validate_yaml(yaml_input, windIO_path / "plant/wind_energy_system.yaml")
            output_dir_name = Path("output_test_foxes")
            output_dir_name.mkdir(parents=True, exist_ok=True)
            try:
                run_foxes(yaml_input, output_dir=output_dir_name, engine=None)
            except Exception as e:
                reset_engine()
                engine = None
                raise e


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


if __name__ == "__main__":
    test_foxes_KUL()
    test_foxes_4wts()
    test_foxes_abl()
    test_foxes_abl_stable()
    test_foxes_profiles()
