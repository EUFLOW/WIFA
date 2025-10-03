import os
from pathlib import Path

import pytest
from windIO import __path__ as wiop
from windIO import load_yaml
from windIO import validate as validate_yaml


# Setup fixtures
@pytest.fixture
def base_path():
    return Path("examples/cases")


@pytest.fixture
def schema_path():
    return Path("plant/wind_energy_system")


@pytest.fixture
def windio_path():
    return Path(wiop[0])


def validate_system_yaml(
    base_path, case_path, windio_path, schema_path, system_fn="system.yaml"
):
    """Helper function to validate a system.yaml file"""
    yaml_fn = "wind_energy_system/" + system_fn
    system_yaml = base_path / case_path / yaml_fn
    validate_yaml(system_yaml, schema_path)


# Test cases
def test_awaken_schema(base_path, windio_path, schema_path):
    """Test AWAKEN case schema"""
    validate_system_yaml(base_path, "AWAKEN", windio_path, schema_path)


def test_kul_les_schema(base_path, windio_path, schema_path):
    """Test KUL_LES case schema"""
    validate_system_yaml(base_path, "KUL_LES", windio_path, schema_path)


def test_kul_les_pywake_schema(base_path, windio_path, schema_path):
    """Test KUL_LES pywake case schema"""
    validate_system_yaml(
        base_path, "KUL_LES", windio_path, schema_path, system_fn="system_pywake.yaml"
    )


def test_open_source_scada_schema(base_path, windio_path, schema_path):
    """Test open source SCADA case schema"""
    validate_system_yaml(base_path, "open_source_scada", windio_path, schema_path)


def test_windio_4turbines_basic_schema(base_path, windio_path, schema_path):
    """Test basic 4 turbines case schema"""
    # Test main system.yaml
    validate_system_yaml(base_path, "windio_4turbines", windio_path, schema_path)

    # Test no_field variant
    system_files = [
        "system.yaml",
        "system_noFlowField.yaml",
        "system_noXYgrid.yaml",
        "system_runAll.yaml",
        "system_zplaneList.yaml",
    ]
    for system in system_files:
        no_field_yaml = base_path / Path(
            f"windio_4turbines/wind_energy_system/{system}"
        )
        validate_yaml(no_field_yaml, schema_path)


def test_windio_4turbines_abl_schema(base_path, windio_path, schema_path):
    """Test 4 turbines ABL case schema"""
    validate_system_yaml(base_path, "windio_4turbines_ABL", windio_path, schema_path)


def test_windio_4turbines_abl_stable_schema(base_path, windio_path, schema_path):
    """Test 4 turbines ABL stable case schema"""
    validate_system_yaml(
        base_path, "windio_4turbines_ABL_stable", windio_path, schema_path
    )


def test_windio_4turbines_multiple_schema(base_path, windio_path, schema_path):
    """Test 4 turbines multiple turbines case schema"""
    validate_system_yaml(
        base_path, "windio_4turbines_multipleTurbines", windio_path, schema_path
    )


def test_windio_4turbines_profiles_stable_schema(base_path, windio_path, schema_path):
    """Test 4 turbines profiles stable case schema"""
    # Test base system.yaml
    base_yaml = (
        base_path / "windio_4turbines_profiles_stable/wind_energy_system/system.yaml"
    )
    validate_yaml(base_yaml, schema_path)

    # Test grid variant
    grid_yaml = (
        base_path
        / "windio_4turbines_profiles_stable/wind_energy_system/system_grid.yaml"
    )
    validate_yaml(grid_yaml, schema_path)


def test_operating_flag_timeseries(base_path, windio_path, schema_path):
    base_yaml = (
        base_path / "timeseries_with_operating_flag/wind_energy_system/system.yaml"
    )
    validate_yaml(base_yaml, schema_path)


def test_simple_wind_rose(base_path, windio_path, schema_path):
    base_yaml = base_path / "simple_wind_rose/wind_energy_system/system.yaml"
    validate_yaml(base_yaml, schema_path)


def test_hetero_wind_rose(base_path, windio_path, schema_path):
    base_yaml = (
        base_path / "heterogeneous_wind_rose_at_turbines/wind_energy_system/system.yaml"
    )
    validate_yaml(base_yaml, schema_path)

    base_yaml = base_path / "heterogeneous_wind_rose_map/wind_energy_system/system.yaml"
    validate_yaml(base_yaml, schema_path)


if __name__ == "__main__":
    # Setup default values that would normally come from fixtures
    base_path = Path("examples/cases")
    schema_path = Path("plant/wind_energy_system")
    windio_path = Path(wiop[0])

    print("Running all schema validation tests...")
    try:
        print("Testing AWAKEN...")
        test_awaken_schema(base_path, windio_path, schema_path)
        print("✓ AWAKEN schema validated")

        print("Testing KUL LES...")
        test_kul_les_schema(base_path, windio_path, schema_path)
        print("✓ KUL_LES schema validated")

        print("Testing KUL LES...")
        test_kul_les_pywake_schema(base_path, windio_path, schema_path)
        print("✓ KUL_LES pywake schema validated")

        print("Testing Open source SCADA...")
        test_open_source_scada_schema(base_path, windio_path, schema_path)
        print("✓ Open source SCADA schema validated")

        print("Testing 4 turbines basic schema...")
        test_windio_4turbines_basic_schema(base_path, windio_path, schema_path)
        print("✓ 4 turbines basic schema validated")

        print("Testing 4 turbines ABL schema...")
        test_windio_4turbines_abl_schema(base_path, windio_path, schema_path)
        print("✓ 4 turbines ABL schema validated")

        print("Testing 4 turbines ABL Stable schema...")
        test_windio_4turbines_abl_stable_schema(base_path, windio_path, schema_path)
        print("✓ 4 turbines ABL stable schema validated")

        print("Testing 4 turbines multiple turbines schema...")
        test_windio_4turbines_multiple_schema(base_path, windio_path, schema_path)
        print("✓ 4 turbines multiple schema validated")

        print("Testing 4 turbines profiles stable schema...")
        test_windio_4turbines_profiles_stable_schema(
            base_path, windio_path, schema_path
        )
        print("✓ 4 turbines profiles stable schema validated")

        print("\nAll schemas validated successfully!")

    except Exception as e:
        print(f"\n❌ Validation failed: {str(e)}")
        raise
