"""
Pytest configuration and fixtures for WIFA tests.

Provides:
- Pre-test cleanup of leftover output directories
- Output directory fixtures with conditional cleanup (preserved on failure)
"""

import shutil
from pathlib import Path

import numpy as np
import pytest

# Store test outcomes for conditional cleanup
_test_outcomes = {}


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_makereport(item, call):
    """Track test outcomes to conditionally preserve output on failure."""
    outcome = yield
    report = outcome.get_result()
    if report.failed:
        _test_outcomes[item.nodeid] = True


def pytest_configure(config):
    """Register custom markers."""
    config.addinivalue_line(
        "markers", "slow: marks tests as slow (deselect with '-m \"not slow\"')"
    )


@pytest.fixture(scope="session", autouse=True)
def cleanup_old_outputs():
    """Remove test output directories from previous test runs at session start."""
    patterns = [
        "output_pywake_*",
        "output_test_*",
    ]
    for pattern in patterns:
        for path in Path(".").glob(pattern):
            if path.is_dir():
                shutil.rmtree(path)
            elif path.is_file():
                path.unlink()
    yield


@pytest.fixture
def output_dir(request, tmp_path):
    """
    Provide a unique temporary output directory for tests.

    Cleans up automatically on test success, preserves on failure for debugging.
    """
    yield tmp_path

    # Only clean up if test passed
    test_failed = _test_outcomes.get(request.node.nodeid, False)
    if not test_failed:
        if tmp_path.exists():
            shutil.rmtree(tmp_path)


@pytest.fixture
def named_output_dir(request):
    """
    Provide a named output directory based on test name.

    Use this when tests need output in the current working directory
    (e.g., for Code Saturne integration).

    Cleans up automatically on test success, preserves on failure.
    """
    test_name = request.node.name.replace("[", "_").replace("]", "_").rstrip("_")
    output_path = Path(f"output_test_{test_name}")
    output_path.mkdir(parents=True, exist_ok=True)

    yield output_path

    # Only clean up if test passed
    test_failed = _test_outcomes.get(request.node.nodeid, False)
    if not test_failed:
        if output_path.exists():
            shutil.rmtree(output_path)


@pytest.fixture
def cleanup_output_dir(request):
    """
    Cleanup fixture for tests that write to the default 'output/' directory.

    Use this when tests can't control their output location (e.g., when YAML
    specifies output_folder). Only cleans up directories created during the test.
    """
    output_path = Path("output")
    existed_before = output_path.exists()

    yield output_path

    # Only clean up if test passed and the directory was created by the test
    test_failed = _test_outcomes.get(request.node.nodeid, False)
    if not test_failed and not existed_before:
        if output_path.exists():
            shutil.rmtree(output_path)


# DTU 10MW turbine data
# (from examples/cases/windio_4turbines/plant_energy_turbine/DTU_10MW_turbine.yaml)
_TURBINE = {
    "name": "DTU 10MW Offshore Reference Turbine",
    "hub_height": 119.0,
    "rotor_diameter": 178.3,
    "performance": {
        "power_curve": {
            "power_wind_speeds": [
                4.0,
                5.0,
                6.0,
                7.0,
                8.0,
                9.0,
                10.0,
                11.0,
                12.0,
                13.0,
                14.0,
                15.0,
                16.0,
                17.0,
                18.0,
                19.0,
                20.0,
                21.0,
                22.0,
                23.0,
                24.0,
                25.0,
            ],
            "power_values": [
                263388.0,
                751154.0,
                1440738.0,
                2355734.0,
                3506858.0,
                4993092.0,
                6849310.0,
                9116402.0,
                10000754.0,
                10009590.0,
                10000942.0,
                10042678.0,
                10003480.0,
                10001600.0,
                10001506.0,
                10013632.0,
                10007428.0,
                10005360.0,
                10002728.0,
                10001130.0,
                10004984.0,
                9997558.0,
            ],
        },
        "Ct_curve": {
            "Ct_wind_speeds": [
                4.0,
                5.0,
                6.0,
                7.0,
                8.0,
                9.0,
                10.0,
                11.0,
                12.0,
                13.0,
                14.0,
                15.0,
                16.0,
                17.0,
                18.0,
                19.0,
                20.0,
                21.0,
                22.0,
                23.0,
                24.0,
                25.0,
            ],
            "Ct_values": [
                0.923,
                0.919,
                0.904,
                0.858,
                0.814,
                0.814,
                0.814,
                0.814,
                0.577,
                0.419,
                0.323,
                0.259,
                0.211,
                0.175,
                0.148,
                0.126,
                0.109,
                0.095,
                0.084,
                0.074,
                0.066,
                0.059,
            ],
        },
    },
}

# 3 turbines in a row, ~7D spacing (7 * 178.3 = 1248.1)
_LAYOUT_X = [0, 1248.1, 2496.2]
_LAYOUT_Y = [0, 0, 0]

# Jensen wake model, minimal analysis config
_ANALYSIS = {
    "wind_deficit_model": {
        "name": "Jensen",
        "wake_expansion_coefficient": {"k_a": 0.0, "k_b": 0.04},
    },
    "deflection_model": {"name": "None"},
    "turbulence_model": {"name": "STF2005", "c1": 1.0, "c2": 1.0},
    "superposition_model": {
        "ws_superposition": "Linear",
        "ti_superposition": "Squared",
    },
    "rotor_averaging": {"name": "Center"},
    "blockage_model": {"name": "None"},
}


def make_timeseries_per_turbine_system_dict(flow_model_name):
    """Build a complete system dict with per-turbine timeseries data including density.

    3 turbines, 6 timesteps, all variables have dims ["time", "wind_turbine"].
    """
    n_times = 6
    n_turbines = 3

    # fmt: off
    ws_data = [
        [9.0, 8.5, 9.2],
        [10.0, 9.8, 10.5],
        [7.5, 7.2, 7.8],
        [8.0, 7.6, 8.3],
        [11.0, 10.8, 11.2],
        [9.5, 9.0, 9.8],
    ]
    wd_data = [
        [270.0, 269.5, 270.5],
        [268.0, 267.5, 268.8],
        [272.0, 271.0, 272.5],
        [270.5, 270.0, 271.0],
        [269.0, 268.5, 269.5],
        [271.0, 270.5, 271.5],
    ]
    ti_data = [
        [0.06, 0.07, 0.05],
        [0.08, 0.09, 0.07],
        [0.05, 0.06, 0.05],
        [0.07, 0.08, 0.06],
        [0.10, 0.09, 0.08],
        [0.06, 0.07, 0.06],
    ]
    # Turbine 0 off at timesteps 2-3
    operating_data = [
        [1, 1, 1],
        [1, 1, 1],
        [0, 1, 1],
        [0, 1, 1],
        [1, 1, 1],
        [1, 1, 1],
    ]
    density_data = [
        [1.225, 1.223, 1.227],
        [1.220, 1.218, 1.222],
        [1.230, 1.228, 1.232],
        [1.228, 1.226, 1.230],
        [1.218, 1.220, 1.222],
        [1.235, 1.233, 1.230],
    ]
    # fmt: on

    return {
        "name": "Dict test: timeseries per-turbine with density",
        "site": {
            "name": "Test site",
            "boundaries": {
                "polygons": [{"x": [-90, 2600, 2600, -90], "y": [90, 90, -90, -90]}]
            },
            "energy_resource": {
                "name": "Test resource",
                "wind_resource": {
                    "time": list(range(n_times)),
                    "wind_turbine": list(range(n_turbines)),
                    "wind_speed": {
                        "data": ws_data,
                        "dims": ["time", "wind_turbine"],
                    },
                    "wind_direction": {
                        "data": wd_data,
                        "dims": ["time", "wind_turbine"],
                    },
                    "turbulence_intensity": {
                        "data": ti_data,
                        "dims": ["time", "wind_turbine"],
                    },
                    "operating": {
                        "data": operating_data,
                        "dims": ["time", "wind_turbine"],
                    },
                    "density": {
                        "data": density_data,
                        "dims": ["time", "wind_turbine"],
                    },
                },
            },
        },
        "wind_farm": {
            "name": "Test farm",
            "layouts": [{"coordinates": {"x": _LAYOUT_X, "y": _LAYOUT_Y}}],
            "turbines": _TURBINE,
        },
        "attributes": {
            "flow_model": {"name": flow_model_name},
            "analysis": _ANALYSIS,
            "model_outputs_specification": {
                "turbine_outputs": {
                    "turbine_nc_filename": "turbine_data.nc",
                    "output_variables": ["power"],
                },
            },
        },
    }
