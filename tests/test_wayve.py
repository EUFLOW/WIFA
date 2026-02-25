import os
from pathlib import Path

import pytest
from windIO import __path__ as wiop
from windIO import validate as validate_yaml

from wifa.wayve_api import run_wayve

test_path = Path(os.path.dirname(__file__))
windIO_path = Path(wiop[0])


@pytest.mark.slow
def test_wayve_4wts(output_dir):
    yaml_input = (
        test_path / "../examples/cases/windio_4turbines/wind_energy_system/system.yaml"
    )
    validate_yaml(yaml_input, Path("plant/wind_energy_system"))
    run_wayve(yaml_input, output_dir=output_dir, debug_mode=True)
