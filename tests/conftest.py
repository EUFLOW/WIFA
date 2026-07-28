import shutil
from pathlib import Path

import numpy as np
import pytest


@pytest.fixture(autouse=True)
def cleanup_test_outputs():
    """Clean up output directories created during tests."""
    yield
    for pattern in ["output_pywake_*", "output_test_*"]:
        for path in Path(".").glob(pattern):
            if path.is_dir():
                shutil.rmtree(path)
    output = Path("output")
    if output.is_dir():
        shutil.rmtree(output)


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
        "wake_expansion_coefficient": {"k_a": 0.04, "k_b": 0.0},
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


def make_mixed_type_timeseries_system_dict(flow_model_name):
    """Build a system dict with TWO turbine types at TWO hub heights.

    3x3 layout (9 turbines), per-turbine timeseries inflow with dims
    ["time", "wind_turbine"] and NO vertical ``height`` profile (the wind
    speeds are already at each turbine's own hub height, as a microscale
    terrain-speedup model would produce). This exercises the mixed
    hub-height + per-turbine path in the pyWake API.

    Type 0 is the DTU 10MW at 119 m; type 1 is a smaller turbine (half the
    power, 120 m rotor, 90 m hub) so both the type assignment and the two
    hub heights are genuinely distinct. Wind speeds are kept in 8-11 m/s,
    away from cut-in/out, so the comparison is unaffected by curve edges.
    """
    n_times = 5
    n_turbines = 9

    # 3x3 grid, ~5D spacing of the larger rotor (5 * 178.3)
    spacing = 5 * _TURBINE["rotor_diameter"]
    axis = [0.0, spacing, 2 * spacing]
    grid_x = [axis[i] for _ in range(3) for i in range(3)]
    grid_y = [axis[j] for j in range(3) for _ in range(3)]

    # Alternating type assignment -> [0, 1, 0, 1, 0, 1, 0, 1, 0]
    turbine_type_idx = [k % 2 for k in range(n_turbines)]

    perf = _TURBINE["performance"]
    pc_ws = perf["power_curve"]["power_wind_speeds"]
    pc_p = perf["power_curve"]["power_values"]
    ct_ws = perf["Ct_curve"]["Ct_wind_speeds"]
    ct_v = perf["Ct_curve"]["Ct_values"]
    type1_power = [0.5 * p for p in pc_p]

    # Deterministic per-turbine timeseries (no RNG), wind from ~270 deg so the
    # rows along x interact.
    base_ws = [8.0, 9.0, 10.0, 8.5, 11.0]
    base_wd = [270.0, 268.0, 272.0, 270.5, 269.0]
    ws_data, wd_data, ti_data = [], [], []
    for t in range(n_times):
        ws_data.append([base_ws[t] + 0.1 * i for i in range(n_turbines)])
        wd_data.append([base_wd[t] + 0.2 * (i - 4) for i in range(n_turbines)])
        ti_data.append([0.06 + 0.002 * i for i in range(n_turbines)])

    return {
        "name": "Dict test: mixed turbine types and hub heights",
        "site": {
            "name": "Test site",
            "boundaries": {
                "polygons": [
                    {
                        "x": [-spacing, 3 * spacing, 3 * spacing, -spacing],
                        "y": [3 * spacing, 3 * spacing, -spacing, -spacing],
                    }
                ]
            },
            "energy_resource": {
                "name": "Test resource",
                "wind_resource": {
                    "time": list(range(n_times)),
                    "wind_turbine": list(range(n_turbines)),
                    "wind_speed": {"data": ws_data, "dims": ["time", "wind_turbine"]},
                    "wind_direction": {
                        "data": wd_data,
                        "dims": ["time", "wind_turbine"],
                    },
                    "turbulence_intensity": {
                        "data": ti_data,
                        "dims": ["time", "wind_turbine"],
                    },
                },
            },
        },
        "wind_farm": {
            "name": "Test farm",
            "layouts": [
                {
                    "coordinates": {"x": grid_x, "y": grid_y},
                    "turbine_types": turbine_type_idx,
                }
            ],
            "turbine_types": {
                0: {
                    "name": "type_0_DTU10MW",
                    "hub_height": 119.0,
                    "rotor_diameter": 178.3,
                    "performance": {
                        "power_curve": {
                            "power_wind_speeds": pc_ws,
                            "power_values": pc_p,
                        },
                        "Ct_curve": {"Ct_wind_speeds": ct_ws, "Ct_values": ct_v},
                    },
                },
                1: {
                    "name": "type_1_small",
                    "hub_height": 90.0,
                    "rotor_diameter": 120.0,
                    "performance": {
                        "power_curve": {
                            "power_wind_speeds": pc_ws,
                            "power_values": type1_power,
                        },
                        "Ct_curve": {"Ct_wind_speeds": ct_ws, "Ct_values": ct_v},
                    },
                },
            },
        },
        "attributes": {
            "flow_model": {"name": flow_model_name},
            "analysis": {
                "wind_deficit_model": {
                    "name": "Bastankhah2014",
                    "wake_expansion_coefficient": {
                        "k_a": 0.04,
                        "k_b": 0.0,
                        "free_stream_ti": False,
                    },
                    "ceps": 0.2,
                    "use_effective_ws": True,
                },
                "axial_induction_model": "Madsen",
                "deflection_model": {"name": "None"},
                "turbulence_model": {"name": "CrespoHernandez"},
                "superposition_model": {
                    "ws_superposition": "Linear",
                    "ti_superposition": "Squared",
                },
                "rotor_averaging": {"name": "Center"},
                "blockage_model": {"name": "None"},
            },
            "model_outputs_specification": {
                "turbine_outputs": {
                    "turbine_nc_filename": "turbine_data.nc",
                    "output_variables": ["power"],
                },
            },
        },
    }


def make_mixed_type_profile_system_dict(flow_model_name):
    """Mixed types/hub heights with a (time, height) vertical wind profile.

    Same farm as :func:`make_mixed_type_timeseries_system_dict` (two types at
    119 m and 90 m), but the inflow is given as a vertical profile over a
    ``height`` dimension that brackets both hub heights, so WIFA must
    interpolate the profile to each turbine's hub height. The wind direction
    veers around 360 deg to exercise circular interpolation.
    """
    system = make_mixed_type_timeseries_system_dict(flow_model_name)

    n_times = 5
    heights = [60.0, 90.0, 120.0, 150.0]
    href, alpha = 119.0, 0.14  # power-law shear reference + exponent
    ws_base = [8.0, 9.0, 10.0, 8.5, 11.0]
    wd_base = [358.0, 359.0, 357.0, 0.5, 359.5]  # straddles 360 deg

    ws_data, wd_data, ti_data = [], [], []
    for t in range(n_times):
        ws_data.append([ws_base[t] * (h / href) ** alpha for h in heights])
        wd_data.append([(wd_base[t] + 0.05 * (h - href)) % 360 for h in heights])
        ti_data.append([max(0.08 - 0.0001 * (h - 60.0), 0.02) for h in heights])

    system["site"]["energy_resource"]["wind_resource"] = {
        "time": list(range(n_times)),
        "height": heights,
        "wind_speed": {"data": ws_data, "dims": ["time", "height"]},
        "wind_direction": {"data": wd_data, "dims": ["time", "height"]},
        "turbulence_intensity": {"data": ti_data, "dims": ["time", "height"]},
    }
    return system


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
