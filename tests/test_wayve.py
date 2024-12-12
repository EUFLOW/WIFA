from flow_api.wayve_api import run_wayve
from windIO.utils.yml_utils import validate_yaml
import os
from pathlib import Path
from windIO import __path__ as wiop


test_path = Path(os.path.dirname(__file__))
windIO_path = Path(wiop[0])


def test_wayve_4wts():
    yaml_input = test_path / "../examples/cases/windio_4turbines/wind_energy_system/system.yaml"
    validate_yaml(yaml_input, windIO_path / "plant/wind_energy_system.yaml")
    output_dir_name = Path("output_test_wayve")
    output_dir_name.mkdir(parents=True, exist_ok=True)
    run_wayve(yaml_input, output_dir=output_dir_name, debug_mode=True)


if __name__ == "__main__":
    test_wayve_4wts()
