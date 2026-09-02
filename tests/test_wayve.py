import os
from pathlib import Path

import pytest

pytest.importorskip(
    "wayve", reason="wayve not installed, install with: pip install wifa[wayve]"
)

from windIO import __path__ as wiop
from windIO import validate as validate_yaml

from wifa.wayve_api import run_wayve

test_path = Path(os.path.dirname(__file__))
windIO_path = Path(wiop[0])


# @pytest.mark.skip()
def test_wayve_4wts():
    yaml_input = (
        test_path / "../examples/cases/windio_4turbines/wind_energy_system/system.yaml"
    )
    validate_yaml(yaml_input, Path("plant/wind_energy_system"))
    output_dir_name = Path("output_test_wayve")
    output_dir_name.mkdir(parents=True, exist_ok=True)
    run_wayve(yaml_input, output_dir=output_dir_name, debug_mode=True)


def test_wayve_multifarm_not_implemented():
    """A wind_farm list must raise a clear NotImplementedError, not a TypeError."""
    yaml_input = (
        test_path
        / "../examples/cases/multiple_wind_farms/wind_energy_system/system.yaml"
    )
    with pytest.raises(NotImplementedError, match="multiple wind farms"):
        run_wayve(yaml_input, output_dir="output_test_wayve_multifarm")


if __name__ == "__main__":
    test_wayve_4wts()
