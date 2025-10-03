from windIO.utils import plant_schemas_path
from windIO.utils.yml_utils import load_yaml, validate_yaml

validate_yaml(
    "cases/windio_4turbines/wind_energy_system/system.yaml",
    plant_schemas_path + "wind_energy_system.yaml",
)
