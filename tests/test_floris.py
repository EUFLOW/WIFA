"""
FLORIS Testing Suite for WIFA Integration

Validates WIFA-generated results against direct FLORIS calculations.
Tests single/multiple turbine farms, wind roses, and operating flags.

Run with: pytest tests/test_floris.py -v
"""

import os
import shutil
from pathlib import Path

import numpy as np
import pytest
import xarray as xr

try:
    import floris
    from floris.turbine_library import build_cosine_loss_turbine_dict
except (ImportError, TypeError):
    pytest.skip(
        "floris not available or incompatible with this Python version",
        allow_module_level=True,
    )
from windIO import __path__ as wiop
from windIO import load_yaml
from windIO import validate as validate_yaml

from wifa.floris_api import run_floris

# ============================================================================ #
# CONFIGURATION
# ============================================================================ #

CWD = Path(__file__).parent.resolve()
save_dir = CWD / ".output_test_floris"
test_path = Path(os.path.dirname(__file__))
windIO_path = Path(wiop[0])

# ============================================================================ #
# PYTEST FIXTURES
# ============================================================================ #


@pytest.fixture(scope="session")
def floris_config():
    """Default FLORIS configuration for all tests."""
    return get_default_floris_dict()


@pytest.fixture
def cleanup_temp_files():
    """Clean up temporary files after each test."""
    yield
    # Cleanup after test
    if save_dir.exists():
        shutil.rmtree(save_dir)


# ============================================================================ #
# UTILITY FUNCTIONS
# ============================================================================ #


def recursive_dict_update(original, updates):
    """Recursively update dictionary with nested structure preservation."""
    for key, value in updates.items():
        if isinstance(value, dict) and key in original:
            recursive_dict_update(original[key], value)
        else:
            original[key] = value
    return original


def get_default_floris_dict():
    """Load default FLORIS configuration."""
    floris_dict_default_yaml = Path(floris.__path__[0]) / "default_inputs.yaml"
    return load_yaml(floris_dict_default_yaml)


def create_output_dir():
    """Create test output directory."""
    output_dir_name = Path("output_test_floris")
    output_dir_name.mkdir(parents=True, exist_ok=True)
    return output_dir_name


def _run_floris(yaml_path, drop_fields=True):
    """Run FLORIS via WIFA, load outputs, and clean up."""
    print(f"\nRUNNING FLORIS ON {yaml_path}\n")

    # Validate and load configuration
    validate_yaml(yaml_path, Path("plant/wind_energy_system"))
    yaml_input = load_yaml(yaml_path)
    yaml_input["attributes"]["model_outputs_specification"]["output_folder"] = save_dir
    yaml_input["attributes"]["flow_model"]["name"] = "floris"

    if drop_fields:
        yaml_input["attributes"]["model_outputs_specification"].pop("flow_field", None)

    # Run simulation
    save_dir.mkdir(parents=True, exist_ok=True)
    run_floris(yaml_input)

    # Load NetCDF outputs
    output_data = {}
    if save_dir.exists():
        for file_path in save_dir.iterdir():
            if file_path.is_file() and file_path.suffix == ".nc":
                output_data[file_path.name] = xr.load_dataset(file_path)

    # Clean up
    if save_dir.exists():
        shutil.rmtree(save_dir)

    return output_data


# ============================================================================ #
# TEST ENVIRONMENT SETUP
# ============================================================================ #


def create_turbine(turbine_name="dtu_10mw"):
    """Create FLORIS turbine dictionary from YAML. Supports DTU and IEA turbines."""
    # Determine YAML path
    if turbine_name == "dtu_10mw":
        turbine_yaml_path = (
            test_path
            / "../examples/cases/windio_4turbines/plant_energy_turbine/DTU_10MW_turbine.yaml"
        )
    elif turbine_name in ["IEA_10MW", "IEA_15MW", "IEA_22MW"]:
        turbine_yaml_path = (
            test_path
            / f"../examples/cases/windio_4turbines_multipleTurbines/plant_energy_turbine/{turbine_name}_turbine.yaml"
        )
    else:
        raise ValueError(
            f"Unknown turbine: {turbine_name}. Options: 'dtu_10mw', 'IEA_10MW', 'IEA_15MW', 'IEA_22MW'"
        )

    turbine_data = load_yaml(turbine_yaml_path)

    # Handle different data formats
    if "power_curve" in turbine_data["performance"]:
        # Direct power curve (DTU)
        wind_speeds = turbine_data["performance"]["power_curve"]["power_wind_speeds"]
        power_values = turbine_data["performance"]["power_curve"]["power_values"]
        ct_values = turbine_data["performance"]["Ct_curve"]["Ct_values"]

        # Auto-convert W to kW
        power_values_kw = (
            [p / 1000.0 for p in power_values]
            if max(power_values) > 100000
            else power_values
        )

    elif "Cp_curve" in turbine_data["performance"]:
        # Calculate from Cp curve (IEA)
        wind_speeds = turbine_data["performance"]["Cp_curve"]["Cp_wind_speeds"]
        cp_values = turbine_data["performance"]["Cp_curve"]["Cp_values"]
        ct_values = turbine_data["performance"]["Ct_curve"]["Ct_values"]

        # P = 0.5 * ρ * A * V³ * Cp
        air_density = 1.225
        rotor_area = np.pi * (turbine_data["rotor_diameter"] / 2) ** 2
        power_values_kw = [
            0.5 * air_density * rotor_area * (ws**3) * cp / 1000
            for ws, cp in zip(wind_speeds, cp_values)
        ]
    else:
        raise ValueError("YAML must contain 'power_curve' or 'Cp_curve'")

    power_thrust_table = {
        "power": power_values_kw,
        "thrust_coefficient": ct_values,
        "wind_speed": wind_speeds,
    }

    return build_cosine_loss_turbine_dict(
        power_thrust_table,
        turbine_name.lower(),
        generator_efficiency=1.0,
        hub_height=turbine_data["hub_height"],
        rotor_diameter=turbine_data["rotor_diameter"],
        TSR=7.0,
    )


def create_four_turbine_farm(layout_type: str):
    """Create 4-turbine farm. 'same' = all DTU, 'multiple' = mixed IEA types."""
    if layout_type == "same":
        turbine_dict = create_turbine("dtu_10mw")
        turbine_dicts = [turbine_dict] * 4
    elif layout_type == "multiple":
        turbine_types = ["IEA_10MW", "IEA_15MW", "IEA_15MW", "IEA_22MW"]
        turbine_dicts = [create_turbine(t_type) for t_type in turbine_types]
    else:
        raise ValueError(
            f"Unknown layout_type: {layout_type}. Use 'same' or 'multiple'."
        )

    return {
        "layout_x": [0, 1248.1, 2496.2, 3744.3],
        "layout_y": [0.0, 0.0, 0.0, 0.0],
        "turbine_type": turbine_dicts,
    }


def create_timeseries():
    """Create 3-step wind time series for testing."""
    return {
        "wind_speeds": [10.091022491455078, 8.797999382019043, 10.307791709899902],
        "wind_directions": [271.8246154785156, 268.6852111816406, 271.0501403808594],
        "turbulence_intensities": [
            2.6189446449279785,
            1.6507437229156494,
            3.7473068237304688,
        ],
    }


def create_timeseries_with_operating_flag():
    """Create time series with turbine availability flags."""
    return {
        "wind_directions": [271.8246, 266.2015, 268.6852, 273.6164, 263.4558, 271.0501],
        "wind_speeds": [10.09102, 10.23302, 8.797999, 9.662098, 9.78371, 10.30779],
        "turbulence_intensities": [
            2.618945,
            1.537051,
            1.650744,
            3.149817,
            0.5826571,
            3.747307,
        ],
        "disable_turbines": np.logical_not(
            np.array(
                [
                    [0, 1, 1, 1],
                    [0, 1, 1, 1],
                    [0, 1, 1, 1],  # T0 offline
                    [0, 1, 1, 1],
                    [1, 1, 1, 1],
                    [1, 1, 1, 1],  # All online
                ]
            )
        ),
    }


def create_simple_wind_rose_energy_resource():
    """Create wind rose with Weibull distributions for 12 sectors."""
    # Sector probabilities and Weibull parameters
    sector_probability = [
        0.0359715203597152,
        0.0394868203948682,
        0.0516739505167395,
        0.0700015407000154,
        0.0836454708364547,
        0.064348500643485,
        0.0864319408643194,
        0.117705101177051,
        0.151575701515757,
        0.147379201473792,
        0.100120501001205,
        0.0516597505165975,
        0.0359715203597152,
    ]

    weibull_a = [
        9.176929,
        9.782334,
        9.531809,
        9.909545,
        10.04269,
        9.593921,
        9.584007,
        10.51499,
        11.39895,
        11.68746,
        11.63732,
        10.08803,
        9.176929,
    ]

    weibull_k = [
        2.392578,
        2.447266,
        2.412109,
        2.591797,
        2.755859,
        2.595703,
        2.583984,
        2.548828,
        2.470703,
        2.607422,
        2.626953,
        2.326172,
        2.392578,
    ]

    wind_direction = [0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360]
    wind_speed = [
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
        17,
        18,
        19,
        20,
        21,
        22,
        23,
        24,
        25,
    ]

    # Normalize probabilities and calculate frequency table
    sector_probability = np.array(sector_probability) / np.sum(sector_probability)
    wind_speeds = np.array(wind_speed)
    freq_table = np.zeros((len(wind_direction), len(wind_speed)))

    for i_dir in range(len(wind_direction)):
        A, k = weibull_a[i_dir], weibull_k[i_dir]
        wind_speed_edges = np.arange(wind_speeds[0] - 0.5, wind_speeds[-1] + 1.5, 1.0)

        # Weibull CDF: F(x) = 1 - exp(-(x/A)^k)
        cdf_edges = 1.0 - np.exp(-((wind_speed_edges / A) ** k))
        cdf_edges[wind_speed_edges < 0] = 0.0

        freq = np.diff(cdf_edges)
        freq_table[i_dir, :] = (freq / freq.sum()) * sector_probability[i_dir]

    wr = floris.WindRose(
        wind_directions=np.array(wind_direction),
        wind_speeds=wind_speeds,
        ti_table=0.1,
        freq_table=freq_table,
    )
    return {"wind_data": wr}


def create_wake_model_config():
    """Standard wake model configuration for testing."""
    return {
        "enable_active_wake_mixing": False,
        "enable_secondary_steering": False,
        "enable_transverse_velocities": False,
        "enable_yaw_added_recovery": False,
        "model_strings": {
            "combination_model": "fls",
            "deflection_model": "none",
            "velocity_model": "gauss",
            "turbulence_model": "crespo_hernandez",
        },
        "wake_velocity_parameters": {"gauss": {"ka": 0.0, "kb": 0.04}},
    }


# ============================================================================ #
# TEST FUNCTIONS
# ============================================================================ #


def test_floris_4wts(floris_config):
    """Test homogeneous 4-turbine farm. Max diff < 1e-6."""
    wes_dir = (
        test_path / "../examples/cases/windio_4turbines/wind_energy_system/system.yaml"
    )
    output_data = _run_floris(wes_dir)
    powers_wifa = output_data["turbine_data.nc"]["power"].values

    floris_dict = floris_config.copy()
    recursive_dict_update(floris_dict["farm"], create_four_turbine_farm("same"))
    recursive_dict_update(floris_dict["flow_field"], create_timeseries())
    recursive_dict_update(floris_dict["wake"], create_wake_model_config())

    fmodel = floris.FlorisModel(floris_dict)
    fmodel.run()
    powers_floris = fmodel.get_turbine_powers()

    # Use pytest.approx for better numerical comparison
    assert powers_wifa == pytest.approx(
        powers_floris, rel=1e-6
    ), f"Power outputs don't match within 1e-6 tolerance"


def test_floris_simple_wind_rose(floris_config):
    """Test wind rose analysis. Max diff < 1e-2."""
    wes_dir = (
        test_path / "../examples/cases/simple_wind_rose/wind_energy_system/system.yaml"
    )
    output_data = _run_floris(wes_dir, drop_fields=True)
    powers_wifa = output_data["PowerTable.nc"]["power"].values

    floris_dict = floris_config.copy()
    recursive_dict_update(floris_dict["farm"], create_four_turbine_farm("same"))
    recursive_dict_update(floris_dict["wake"], create_wake_model_config())

    fmodel = floris.FlorisModel(floris_dict)
    fmodel.set(**create_simple_wind_rose_energy_resource())
    fmodel.run()
    powers_floris = fmodel.get_turbine_powers()

    assert powers_wifa == pytest.approx(
        powers_floris, rel=1e-2
    ), f"Wind rose power outputs don't match within 1e-2 tolerance"


def test_floris_timeseries_with_operating_flag(floris_config):
    """Test turbine availability flags. Max diff < 1e-2."""
    wes_dir = (
        test_path
        / "../examples/cases/timeseries_with_operating_flag/wind_energy_system/system.yaml"
    )
    output_data = _run_floris(wes_dir, drop_fields=True)
    powers_wifa = output_data["turbine_data.nc"]["power"].values

    floris_dict = floris_config.copy()
    recursive_dict_update(floris_dict["farm"], create_four_turbine_farm("same"))
    recursive_dict_update(floris_dict["wake"], create_wake_model_config())

    fmodel = floris.FlorisModel(floris_dict)
    fmodel.set_operation_model("mixed")
    fmodel.set(**create_timeseries_with_operating_flag())
    fmodel.run()
    powers_floris = fmodel.get_turbine_powers()

    assert powers_wifa == pytest.approx(
        powers_floris, rel=1e-2
    ), f"Operating flag power outputs don't match within 1e-2 tolerance"


def test_floris_multiple_turbines(floris_config):
    """Test mixed turbine types. Max diff < 0.1. Marked as slow test."""
    wes_dir = (
        test_path
        / "../examples/cases/windio_4turbines_multipleTurbines/wind_energy_system/system.yaml"
    )
    output_data = _run_floris(wes_dir, drop_fields=True)
    powers_wifa = output_data["turbine_data.nc"]["power"].values

    floris_dict = floris_config.copy()
    recursive_dict_update(floris_dict["farm"], create_four_turbine_farm("multiple"))
    recursive_dict_update(floris_dict["flow_field"], create_timeseries())
    recursive_dict_update(floris_dict["wake"], create_wake_model_config())

    # Set reference height to mean hub height
    mean_hub_height = np.mean(
        [t["hub_height"] for t in floris_dict["farm"]["turbine_type"]]
    )
    floris_dict["flow_field"]["reference_wind_height"] = mean_hub_height

    fmodel = floris.FlorisModel(floris_dict)
    fmodel.run()
    powers_floris = fmodel.get_turbine_powers()

    max_diff = np.max(np.abs(powers_wifa - powers_floris)) / np.max(np.abs(powers_wifa))
    print(f"Max relative difference: {max_diff:.6f}")
    assert powers_wifa == pytest.approx(
        powers_floris, rel=0.1
    ), f"Mixed turbine power outputs don't match within 0.1 tolerance"


# ============================================================================ #
# PYTEST FIXTURES
# ============================================================================ #


@pytest.fixture(scope="session")
def test_output_dir():
    """Shared test output directory fixture with automatic cleanup."""
    output_dir = Path("output_test_floris")
    output_dir.mkdir(parents=True, exist_ok=True)
    yield output_dir
    if output_dir.exists():
        shutil.rmtree(output_dir)


if __name__ == "__main__":
    pytest.main([str(CWD / "test_floris.py"), "-v"])
