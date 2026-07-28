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


def test_wayve_k_convention_lanzilao():
    """windIO defines k = k_a + k_b*TI; wayve's Lanzilao computes
    kwake = ka*TI + kb, so the reader must map swapped: ka <- k_b, kb <- k_a.
    Regression for the unswapped pass-through that silently turned the
    constant expansion into a TI multiplier (and vice versa)."""
    from wifa.wayve_api import wake_model_setup

    wm = wake_model_setup(
        {
            "wind_deficit_model": {
                "wake_expansion_coefficient": {"k_a": 0.004, "k_b": 0.38},
                "ceps": 0.2,
            }
        }
    )
    assert wm.ka == 0.38  # Lanzilao TI multiplier holds windIO k_b
    assert wm.kb == 0.004  # Lanzilao constant holds windIO k_a


def test_wayve_scalar_k_is_constant_expansion():
    """A scalar k is a constant expansion: Lanzilao kb, with no TI term."""
    from wifa.wayve_api import wake_model_setup

    wm = wake_model_setup(
        {
            "wind_deficit_model": {
                "wake_expansion_coefficient": {"k": 0.05},
                "ceps": 0.2,
            }
        }
    )
    assert wm.ka == 0.0
    assert wm.kb == 0.05


_LANZILAO_ANALYSIS = {
    "wind_deficit_model": {
        "wake_expansion_coefficient": {"k_a": 0.004, "k_b": 0.38},
        "ceps": 0.2,
    }
}


def test_wake_tool_read_from_wm_coupling():
    """The windIO schema nests `wake_tool` under `wm_coupling` (a top-level
    `analysis.wake_tool` fails validation). Reading only the analysis level made
    the foxes coupling unreachable from any schema-valid file."""
    from wifa.wayve_api import wake_model_setup

    analysis = {**_LANZILAO_ANALYSIS, "wm_coupling": {"wake_tool": "wayve"}}
    wm = wake_model_setup(analysis)
    assert type(wm).__name__ == "Lanzilao"

    analysis = {**_LANZILAO_ANALYSIS, "wm_coupling": {"wake_tool": "nonsense"}}
    with pytest.raises(NotImplementedError, match="nonsense"):
        wake_model_setup(analysis)


def test_wake_tool_analysis_level_fallback():
    """Legacy hand-written yamls put `wake_tool` at the analysis level."""
    from wifa.wayve_api import wake_model_setup

    with pytest.raises(NotImplementedError, match="nonsense"):
        wake_model_setup({**_LANZILAO_ANALYSIS, "wake_tool": "nonsense"})


def test_wake_tool_defaults_to_wayve():
    from wifa.wayve_api import wake_model_setup

    wm = wake_model_setup(dict(_LANZILAO_ANALYSIS))
    assert type(wm).__name__ == "Lanzilao"


def _scalar_resource(ti):
    """Minimal single-point (no `height`) wind resource for flow_io_abl."""
    return {
        "wind_speed": {"data": [9.0]},
        "wind_direction": {"data": [270.0]},
        "turbulence_intensity": {"data": [ti]},
        "z0": {"data": [0.03]},
    }


def test_scalar_resource_ti_is_a_fraction():
    """windIO carries TI as a fraction. The scalar branch used to divide it by
    100 (and read it only when `z0` happened to be present), collapsing the
    TI-proportional wake expansion by two orders of magnitude."""
    from wifa.wayve_api import flow_io_abl

    abl, _ = flow_io_abl(_scalar_resource(0.08), 0, zh=78.0, h1=156.0)
    assert abl.TI == pytest.approx(0.08)


def test_scalar_resource_ti_read_without_z0():
    """TI must be read whenever it is present, not only when `z0` is."""
    from wifa.wayve_api import flow_io_abl

    resource = _scalar_resource(0.08)
    del resource["z0"]
    abl, _ = flow_io_abl(resource, 0, zh=78.0, h1=156.0)
    assert abl.TI == pytest.approx(0.08)


def test_scalar_resource_ti_defaults_when_absent():
    from wifa.wayve_api import flow_io_abl

    resource = _scalar_resource(0.08)
    del resource["turbulence_intensity"]
    abl, _ = flow_io_abl(resource, 0, zh=78.0, h1=156.0)
    assert abl.TI == pytest.approx(0.04)


def test_scalar_resource_air_density_is_per_state():
    """A per-state `air_density` in the wind resource sets abl.rho, which
    linearly scales the reported turbine power (the APM velocity solution
    does not depend on it)."""
    from wifa.wayve_api import flow_io_abl

    resource = _scalar_resource(0.08)
    resource["air_density"] = {"data": [1.29, 1.11]}
    resource["wind_speed"]["data"].append(9.0)
    resource["wind_direction"]["data"].append(270.0)
    resource["turbulence_intensity"]["data"].append(0.08)
    resource["z0"]["data"].append(0.03)
    assert flow_io_abl(resource, 0, zh=78.0, h1=156.0)[0].rho == pytest.approx(1.29)
    assert flow_io_abl(resource, 1, zh=78.0, h1=156.0)[0].rho == pytest.approx(1.11)


def test_scalar_resource_air_density_defaults_when_absent():
    from wifa.wayve_api import flow_io_abl

    abl, _ = flow_io_abl(_scalar_resource(0.08), 0, zh=78.0, h1=156.0)
    assert abl.rho == pytest.approx(1.225)


def _timeseries_system(wind_speeds, subset):
    """Single-turbine system whose wind speed differs at every timestep."""
    n = len(wind_speeds)
    return {
        "site": {
            "energy_resource": {
                "wind_resource": {
                    "time": list(range(n)),
                    "wind_speed": {"data": list(wind_speeds)},
                    "wind_direction": {"data": [270.0] * n},
                    "turbulence_intensity": {"data": [0.08] * n},
                    "z0": {"data": [0.03] * n},
                    "fc": {"data": [1.0e-4] * n},
                }
            }
        },
        "wind_farm": {
            "layouts": [{"coordinates": {"x": [0.0], "y": [0.0]}}],
            "turbines": {
                "performance": {
                    "power_curve": {
                        "power_values": [0, 1e6, 3e6, 3e6],
                        "power_wind_speeds": [3, 8, 12, 25],
                    },
                    "Ct_curve": {
                        "Ct_values": [0.8, 0.8, 0.4, 0.2],
                        "Ct_wind_speeds": [3, 8, 12, 25],
                    },
                },
                "hub_height": 80.0,
                "rotor_diameter": 80.0,
            },
        },
        "attributes": {
            "flow_model": {"name": "wayve"},
            "analysis": {
                "wind_deficit_model": {
                    "wake_expansion_coefficient": {"k_a": 0.04, "k_b": 0.0},
                    "ceps": 0.2,
                },
                "superposition_model": {"ws_superposition": "Product"},
                "wm_coupling": {"method": "PB"},
            },
            "model_outputs_specification": {
                "output_folder": "output",
                "run_configuration": {
                    "times_run": {"all_occurences": False, "subset": subset}
                },
                "turbine_outputs": {
                    "turbine_nc_filename": "turbine_data.nc",
                    "output_variables": ["rotor_effective_velocity"],
                },
            },
        },
    }


def test_times_run_subset_selects_the_requested_rows(tmp_path):
    """Regression: the subsetted timestamps were re-enumerated from 0, so the
    ABL was built from rows 0..n-1 of the wind resource while being *labelled*
    with the requested timestamps — silently simulating the wrong states."""
    import xarray as xr

    from wifa.wayve_api import run_wayve

    wind_speeds = [6.0, 9.0, 12.0]
    subset = [2]
    # debug_mode skips the APM solve: the wake model still sees the full ABL.
    run_wayve(
        _timeseries_system(wind_speeds, subset),
        output_dir=tmp_path,
        debug_mode=True,
    )

    ds = xr.load_dataset(tmp_path / "turbine_data.nc")
    assert ds["states"].values.tolist() == subset
    # A lone turbine is unwaked, so its rotor-effective velocity is the
    # background speed of the *requested* state (12 m/s), not of row 0 (6 m/s).
    rews = float(ds["rotor_effective_velocity"].isel(states=0, turbine=0))
    assert rews == pytest.approx(wind_speeds[subset[0]], rel=0.02)


if __name__ == "__main__":
    test_wayve_4wts()


# ---------------------------------------------------------------------------
# Solver-frame rotation (westerly normalization) and profile-resolving
# free atmosphere
# ---------------------------------------------------------------------------

_CT_CURVE = {
    "Ct_values": [0.8, 0.8, 0.4, 0.2],
    "Ct_wind_speeds": [3, 8, 12, 25],
}
_POWER_CURVE = {
    "power_values": [0, 1e6, 3e6, 3e6],
    "power_wind_speeds": [3, 8, 12, 25],
}


def _solver_system(x, y, wind_resource, flow_field=None, layers_description=None):
    """Multi-turbine system for full (non-debug) APM solves on a small grid."""
    analysis = {
        "wind_deficit_model": {
            "wake_expansion_coefficient": {"k_a": 0.04, "k_b": 0.0},
            "ceps": 0.2,
        },
        "superposition_model": {"ws_superposition": "Product"},
        "wm_coupling": {"method": "PB"},
        "apm_grid": {"Lx": 1.0e5, "Ly": 1.0e5, "dx": 1.0e3},
    }
    if layers_description is not None:
        analysis["layers_description"] = layers_description
    outputs = {
        "output_folder": "output",
        "run_configuration": {"times_run": {"all_occurences": True}},
        "turbine_outputs": {
            "turbine_nc_filename": "turbine_data.nc",
            "output_variables": ["power", "rotor_effective_velocity"],
        },
    }
    if flow_field is not None:
        outputs["flow_field"] = flow_field
    return {
        "site": {"energy_resource": {"wind_resource": wind_resource}},
        "wind_farm": {
            "layouts": [{"coordinates": {"x": list(x), "y": list(y)}}],
            "turbines": {
                "performance": {"power_curve": _POWER_CURVE, "Ct_curve": _CT_CURVE},
                "hub_height": 80.0,
                "rotor_diameter": 80.0,
            },
        },
        "attributes": {
            "flow_model": {"name": "wayve"},
            "analysis": analysis,
            "model_outputs_specification": outputs,
        },
    }


def _scalar_wind_resource(wd, n=1):
    return {
        "time": list(range(n)),
        "wind_speed": {"data": [9.0] * n},
        "wind_direction": {"data": [float(wd)] * n},
        "turbulence_intensity": {"data": [0.08] * n},
        "z0": {"data": [0.03] * n},
        "fc": {"data": [1.0e-4] * n},
    }


def _flow_field_spec(x_bounds, y_bounds, n):
    return {
        "report": True,
        "flow_nc_filename": "flow_field.nc",
        "output_variables": ["wind_speed", "wind_direction"],
        "z_planes": {
            "xy_sampling": "grid",
            "x_bounds": list(x_bounds),
            "y_bounds": list(y_bounds),
            "Nx": n[0],
            "Ny": n[1],
            "z_sampling": "hub_heights",
        },
    }


def test_scalar_resource_is_solver_frame_aligned():
    """The ABL comes out with the hub-height wind along +x (wayve's westerly
    convention) for any input direction, and the undone rotation is
    reported."""
    import numpy as np

    from wifa.wayve_api import flow_io_abl

    resource = _scalar_resource(0.08)
    resource["wind_direction"]["data"] = [225.0]
    abl, rotation = flow_io_abl(resource, 0, zh=78.0, h1=156.0)
    assert rotation == pytest.approx(np.deg2rad(45.0))
    u_hub = np.interp(78.0, abl.zs, abl.us)
    v_hub = np.interp(78.0, abl.zs, abl.vs)
    assert u_hub > 0.0
    assert v_hub == pytest.approx(0.0, abs=1.0e-8 * u_hub)


@pytest.fixture(scope="module")
def _rotation_pair(tmp_path_factory):
    """Physically identical cases under rotation: a westerly-aligned turbine
    row (a); the same row and grids rotated 90 degrees with a southerly wind,
    run as two identical states (b); and the row rotated by a generic 37
    degrees (c). 90 degrees maps the square FFT grid onto itself, so only the
    generic angle pins the solver-frame rotation itself (see the tests)."""
    import numpy as np
    import xarray as xr

    from wifa.wayve_api import run_wayve

    span = 400.0
    xa = [0.0, span, 2 * span]
    ff_a = _flow_field_spec((-1000.0, 1800.0), (-1200.0, 1200.0), (15, 13))
    ff_b = _flow_field_spec((-1200.0, 1200.0), (-1000.0, 1800.0), (13, 15))
    system_a = _solver_system(xa, [0.0, 0.0, 0.0], _scalar_wind_resource(270.0), ff_a)
    system_b = _solver_system(
        [0.0, 0.0, 0.0], xa, _scalar_wind_resource(180.0, n=2), ff_b
    )
    # Case c: rotate the westerly row into the frame of a 233 deg (270 - 37)
    # wind: earth layout = R(+rotation) applied to the westerly layout.
    rot_c = np.deg2rad(270.0 - 233.0)
    c, s = np.cos(rot_c), np.sin(rot_c)
    xc = [c * x for x in xa]
    yc = [s * x for x in xa]
    system_c = _solver_system(xc, yc, _scalar_wind_resource(233.0))
    out_a = tmp_path_factory.mktemp("wayve_rot_a")
    out_b = tmp_path_factory.mktemp("wayve_rot_b")
    out_c = tmp_path_factory.mktemp("wayve_rot_c")
    run_wayve(system_a, output_dir=out_a)
    run_wayve(system_b, output_dir=out_b)
    run_wayve(system_c, output_dir=out_c)
    return {
        "turbines_a": xr.load_dataset(out_a / "turbine_data.nc"),
        "turbines_b": xr.load_dataset(out_b / "turbine_data.nc"),
        "turbines_c": xr.load_dataset(out_c / "turbine_data.nc"),
        "flow_a": xr.load_dataset(out_a / "flow_field.nc"),
        "flow_b": xr.load_dataset(out_b / "flow_field.nc"),
    }


def test_rotational_invariance_of_turbine_outputs(_rotation_pair):
    """A southerly case with the layout turned 90 degrees is the same physical
    system as the westerly case, so per-turbine outputs must match. Before the
    solver-frame rotation, the gravity-wave solve saw the farm footprint (and
    any anisotropic grid) at the wrong angle for non-westerly states."""
    import numpy as np

    power_a = _rotation_pair["turbines_a"]["power"].values
    power_b = _rotation_pair["turbines_b"]["power"].values
    assert np.all(np.isfinite(power_a)) and np.all(power_a > 0.0)
    np.testing.assert_allclose(power_b[0], power_a[0], rtol=1.0e-6)
    rews_a = _rotation_pair["turbines_a"]["rotor_effective_velocity"].values
    rews_b = _rotation_pair["turbines_b"]["rotor_effective_velocity"].values
    np.testing.assert_allclose(rews_b[0], rews_a[0], rtol=1.0e-6)
    # Downstream turbines are waked in both frames
    assert power_a[0, 0] > power_a[0, -1]


def test_rotational_invariance_at_generic_angle(_rotation_pair):
    """A generic 37-degree rotation is what actually pins the solver-frame
    normalization: 90-degree rotations map the square FFT grid onto itself,
    so the pre-rotation code already passed the 90-degree pair exactly, while
    it fails this case at ~4e-4 relative. At a generic angle both cos and sin
    terms of every rotation formula are exercised (90 deg blinds the cos
    terms, 45 deg hides cos/sin swaps)."""
    import numpy as np

    power_a = _rotation_pair["turbines_a"]["power"].values
    power_c = _rotation_pair["turbines_c"]["power"].values
    np.testing.assert_allclose(power_c, power_a, rtol=1.0e-6)


def test_repeated_states_give_identical_output(_rotation_pair):
    """Two identical southerly states must give identical per-turbine output.
    Guards the per-state wind-farm rebuild: a cumulative in-place rotation of
    farm_dat would rotate state 1 twice yet still produce finite, positive,
    plausibly-waked power."""
    import numpy as np

    power = _rotation_pair["turbines_b"]["power"]
    np.testing.assert_allclose(
        power.isel(states=1).values, power.isel(states=0).values, rtol=1.0e-12
    )


def test_flow_field_is_reported_in_the_earth_frame(_rotation_pair):
    """flow_field.nc stays on the requested east/north grid: the southerly
    case reports ~180 deg background flow, with the wake north (downstream)
    of the farm."""
    import numpy as np

    ds = _rotation_pair["flow_b"].isel(states=0, z=0)
    # Upstream (south) edge, away from the farm axis: background direction
    wd_upstream = float(
        ds["wind_direction"].sel(x=-1200.0, y=-1000.0, method="nearest")
    )
    assert wd_upstream % 360.0 == pytest.approx(180.0, abs=3.0)
    # The deepest deficit sits downstream of the last turbine (north)
    ws = ds["wind_speed"]
    loc = ws.argmin(dim=["x", "y"])
    assert float(ws.y[int(loc["y"])]) >= 800.0
    assert abs(float(ws.x[int(loc["x"])])) < 400.0


def test_flow_fields_map_onto_each_other_under_rotation(_rotation_pair):
    """The southerly flow field equals the westerly one rotated by 90 deg:
    ws_b(x, y) == ws_a(y, -x) on the matching grids."""
    import numpy as np

    ws_a = _rotation_pair["flow_a"]["wind_speed"].isel(states=0, z=0).values
    ws_b = _rotation_pair["flow_b"]["wind_speed"].isel(states=0, z=0).values
    # grids: a is (15 x, 13 y), b is (13 x, 15 y); earth point (x, y) in b
    # maps to solver/earth point (y, -x) in a; x_b symmetric -> index flip
    mapped = np.transpose(ws_a)[::-1, :]
    np.testing.assert_allclose(ws_b, mapped, rtol=1.0e-6)
    wd_a = _rotation_pair["flow_a"]["wind_direction"].isel(states=0, z=0).values
    wd_b = _rotation_pair["flow_b"]["wind_direction"].isel(states=0, z=0).values
    mapped_wd = (np.transpose(wd_a)[::-1, :] - 90.0) % 360.0
    np.testing.assert_allclose(wd_b % 360.0, mapped_wd, rtol=0, atol=1.0e-6)


def test_xy_points_matches_xy_plane():
    """_xy_plane_points on an axis-aligned meshgrid reproduces wayve's
    xy_plane exactly (it mirrors the same code with the meshgrid hoisted
    out) — the guard against drift between the fork helper and wayve."""
    import numpy as np
    from wayve.apm import APM
    from wayve.grid.grid import Stat2Dgrid
    from wayve.momentum_flux_parametrizations import FrictionCoefficients
    from wayve.solvers import FixedPointIteration

    from wifa.wayve_api import (
        _pressure_for_state,
        _xy_plane_points,
        flow_io_abl,
        wf_setup,
    )

    system = _solver_system(
        [0.0, 400.0, 800.0], [0.0, 0.0, 0.0], _scalar_wind_resource(270.0)
    )
    farm_dat = system["wind_farm"]
    analysis_dat = system["attributes"]["analysis"]
    resource = system["site"]["energy_resource"]["wind_resource"]
    abl, rotation = flow_io_abl(resource, 0, zh=80.0, h1=160.0)
    assert rotation == pytest.approx(0.0)
    wind_farm, forcing, _, _ = wf_setup(
        farm_dat, analysis_dat, 1.0e3, rotation=rotation
    )
    grid = Stat2Dgrid(1.0e5, 100, 1.0e5, 100)
    model = APM(grid, forcing, abl, FrictionCoefficients(), _pressure_for_state(abl, 1))
    model.solve(method=FixedPointIteration(tol=5.0e-3, relax=0.7))
    coupling = wind_farm.coupling
    wake_model = coupling.wake_model
    u_bg_evaluator = coupling.set_up_u_bg_evaluator(abl)
    apm_evaluator = coupling.apm_evaluator
    xs = np.linspace(-2000.0, 4000.0, 25)
    ys = np.linspace(-1500.0, 1500.0, 21)
    ref = wake_model.xy_plane(
        wind_farm, abl, u_bg_evaluator, apm_evaluator, xs, ys, 80.0
    )
    Xs, Ys = np.meshgrid(xs, ys, indexing="ij")
    new = _xy_plane_points(
        wake_model, wind_farm, abl, u_bg_evaluator, apm_evaluator, Xs, Ys, 80.0
    )
    for ref_field, new_field in zip(ref, new):
        np.testing.assert_array_equal(new_field, ref_field)


def _profile_wind_resource(wd_hub=225.0, veer=10.0, n=1, spelling="flat"):
    """Truncated (1.7 km) mesoscale-style profile resource: winds/theta
    profiles plus ERA5-style surface scalars and an explicit CI block.

    ``spelling`` selects how the capping inversion is written: ``"flat"``
    (the windIO schema keys, which is what a netCDF ``!include`` can carry)
    or ``"nested"`` (the older thermal_stratification block). The two must
    produce identical ABL states — see
    ``test_flat_and_nested_capping_inversion_agree``."""
    import numpy as np

    zs = np.linspace(10.0, 1700.0, 18)
    ws = 8.0 + 3.0 * zs / 1700.0
    # Small linear veer across the profile around the hub-height direction
    hub_frac = np.interp(80.0, zs, np.linspace(0.0, 1.0, len(zs)))
    wd = wd_hub + veer * (np.linspace(0.0, 1.0, len(zs)) - hub_frac)
    theta = np.where(
        zs < 1000.0,
        290.0,
        np.where(
            zs <= 1200.0,
            290.0 + 3.0 * (zs - 1000.0) / 200.0,
            293.0 + 3.5e-3 * (zs - 1200.0),
        ),
    )
    per_state = lambda arr: {"data": [list(arr)] * n}  # noqa: E731
    scalar = lambda val: {"data": [val] * n}  # noqa: E731
    if spelling == "flat":
        ci = {
            "ABL_height": scalar(1100.0),
            "capping_inversion_thickness": scalar(200.0),
            "capping_inversion_strength": scalar(3.0),
            "lapse_rate": scalar(3.5e-3),
        }
    elif spelling == "nested":
        ci = {
            "thermal_stratification": {
                "capping_inversion": {
                    "ABL_height": scalar(1100.0),
                    "dH": scalar(200.0),
                    "dtheta": scalar(3.0),
                    "lapse_rate": scalar(3.5e-3),
                }
            }
        }
    else:
        raise ValueError(f"unknown spelling {spelling!r}")
    return {
        "time": list(range(n)),
        "height": list(zs),
        "wind_speed": per_state(ws),
        "wind_direction": per_state(wd),
        "potential_temperature": per_state(theta),
        "turbulence_intensity": scalar(0.08),
        "air_density": scalar(1.2),
        "z0": scalar(0.03),
        "fc": scalar(1.0e-4),
        "friction_velocity": scalar(0.35),
        "boundary_layer_height": scalar(900.0),
        "LMO": scalar(5.0e3),
        **ci,
    }


def test_flat_and_nested_capping_inversion_agree():
    """The two windIO spellings of the capping inversion must be the same input.

    The flat keys are the windIO schema spelling (and what code_saturne's
    adapter already reads); the nested thermal_stratification block is the
    older one, kept so existing files keep working. A netCDF ``!include``
    can only carry the flat one, so the whole point is that switching a file
    over changes nothing about the atmosphere it describes.
    """
    import numpy as np

    from wifa.wayve_api import flow_io_abl

    flat, _ = flow_io_abl(_profile_wind_resource(spelling="flat"), 0, zh=80.0, h1=160.0)
    nested, _ = flow_io_abl(
        _profile_wind_resource(spelling="nested"), 0, zh=80.0, h1=160.0
    )
    # The four scalars drive exactly these: the upper-layer depth, the
    # reduced gravity, the free-atmosphere stratification and the inversion
    # bounds. Exact equality — the values travel through unchanged.
    for field in ("H", "H2", "gprime", "N", "inv_bottom", "inv_top", "U1", "U2", "U3"):
        assert getattr(flat, field) == getattr(nested, field), field
    for field in ("zs", "us", "vs", "tauxs", "tauys"):
        np.testing.assert_array_equal(
            getattr(flat, field), getattr(nested, field), err_msg=field
        )


def test_capping_inversion_rejects_both_spellings_at_once():
    """A file carrying both spellings is half-migrated; picking one silently
    would pick a physics answer for the user."""
    from wifa.wayve_api import capping_inversion_spelling

    resource = _profile_wind_resource(spelling="flat")
    resource |= {
        k: v
        for k, v in _profile_wind_resource(spelling="nested").items()
        if k == "thermal_stratification"
    }
    with pytest.raises(ValueError, match="twice"):
        capping_inversion_spelling(resource)


@pytest.mark.parametrize(
    "spelling, dropped",
    [
        ("flat", "capping_inversion_strength"),
        ("flat", "ABL_height"),
        ("nested", "dtheta"),
    ],
)
def test_capping_inversion_rejects_a_partial_block(spelling, dropped):
    """A partial block is not completed from defaults.

    dtheta alone sets the amplitude of the gravity-wave forcing, so pairing a
    real ABL_height with an invented one yields a result that looks fitted and
    is not. (WIFA's code_saturne adapter does fill per-key defaults; that
    behaviour is not copied here deliberately.)
    """
    from wifa.wayve_api import capping_inversion_spelling

    resource = _profile_wind_resource(spelling=spelling)
    if spelling == "flat":
        del resource[dropped]
    else:
        del resource["thermal_stratification"]["capping_inversion"][dropped]
    with pytest.raises(ValueError, match="Incomplete"):
        capping_inversion_spelling(resource)


def test_no_capping_inversion_falls_back_to_the_profile_fit():
    """Omitting all four scalars is a valid, distinct request: fit the
    inversion from the potential-temperature profile."""
    from wifa.wayve_api import capping_inversion_spelling, flow_io_abl

    resource = _profile_wind_resource(spelling="flat")
    for key in (
        "ABL_height",
        "capping_inversion_thickness",
        "capping_inversion_strength",
        "lapse_rate",
    ):
        del resource[key]
    assert capping_inversion_spelling(resource) is None

    fitted, _ = flow_io_abl(resource, 0, zh=80.0, h1=160.0)
    explicit, _ = flow_io_abl(
        _profile_wind_resource(spelling="flat"), 0, zh=80.0, h1=160.0
    )
    # The fit finds the 1000-1200 m jump in the test profile, not the
    # explicit block's 1100 +/- 100 m; if these agreed the test would be moot.
    assert fitted.H != explicit.H


def test_truncated_profile_free_atmosphere_capped_at_data_top():
    """The free atmosphere of a truncated profile is capped at the top of the
    data: NonUniform's layers then span actual profile points instead of the
    default 10 km tropopause (constant-fill splines above ~1.7 km)."""
    import numpy as np

    from wifa.wayve_api import flow_io_abl

    abl, rotation = flow_io_abl(_profile_wind_resource(), 0, zh=80.0, h1=160.0)
    assert abl.h_strat == pytest.approx(1700.0)
    assert abl.Uinf == pytest.approx(abl.us[-1])
    assert abl.Vinf == pytest.approx(abl.vs[-1])
    assert abl.Ninf == pytest.approx(abl.N)
    # Solver frame: hub-height wind along +x, veer preserved aloft.
    # (Interpolating components vs. directions differs by ~1e-3 m/s here.)
    u_hub = np.interp(80.0, abl.zs, abl.us)
    v_hub = np.interp(80.0, abl.zs, abl.vs)
    assert u_hub > 0.0
    assert v_hub == pytest.approx(0.0, abs=1.0e-3)
    assert rotation == pytest.approx(np.deg2rad(270.0 - 225.0), abs=0.02)
    assert np.abs(abl.vs).max() > 0.1  # veer survives the rotation


def test_pressure_for_state_selects_nonuniform_only_with_data():
    """NonUniform needs profile points between the inversion top and the
    (capped) free-atmosphere top; scalar states fall back to Uniform."""
    from wifa.wayve_api import _pressure_for_state, flow_io_abl

    abl_profile, _ = flow_io_abl(_profile_wind_resource(), 0, zh=80.0, h1=160.0)
    assert type(_pressure_for_state(abl_profile, 6)).__name__ == "NonUniform"
    assert type(_pressure_for_state(abl_profile, 1)).__name__ == "Uniform"
    # Scalar branch: the Nieuwstadt profile stops at the inversion -> no
    # data-backed free atmosphere -> warn and fall back
    abl_scalar, _ = flow_io_abl(_scalar_resource(0.08), 0, zh=78.0, h1=156.0)
    with pytest.warns(UserWarning, match="falling back to the Uniform"):
        assert type(_pressure_for_state(abl_scalar, 6)).__name__ == "Uniform"


def test_run_wayve_profile_nonuniform_end_to_end(tmp_path):
    """Non-westerly mesoscale profiles with a profile-resolving free
    atmosphere run end-to-end and produce finite, waked turbine output —
    and NonUniform is genuinely active (no silent Uniform fallback), so the
    number_of_fa_layers wiring through run_wayve is pinned, not just the
    _pressure_for_state unit behavior."""
    import warnings as _warnings

    import numpy as np
    import xarray as xr

    from wifa.wayve_api import run_wayve

    # Row aligned with the 225 deg flow (SW): downstream is to the NE
    span = 400.0 / np.sqrt(2.0)
    system = _solver_system(
        [0.0, span, 2 * span],
        [0.0, span, 2 * span],
        _profile_wind_resource(wd_hub=225.0),
        layers_description={"farm_layer_height": 160.0, "number_of_fa_layers": 6},
    )
    with _warnings.catch_warnings(record=True) as caught:
        _warnings.simplefilter("always")
        run_wayve(system, output_dir=tmp_path)
    fallbacks = [w for w in caught if "falling back to the Uniform" in str(w.message)]
    assert not fallbacks, "NonUniform silently fell back to Uniform"
    ds = xr.load_dataset(tmp_path / "turbine_data.nc")
    power = ds["power"].values
    assert power.shape == (1, 3)
    assert np.all(np.isfinite(power)) and np.all(power > 0.0)
    # SW flow: the NE-most turbine is waked
    assert power[0, 0] > power[0, -1]


# Full windIO-shaped analysis block for the foxes wake tool (mirrors what the
# flow_model_chain Experiment serializer emits; foxes' windIO reader needs the
# named model blocks, which the Lanzilao path simply ignores).
_FOXES_ANALYSIS = {
    "wind_deficit_model": {
        "name": "Bastankhah2014",
        "wake_expansion_coefficient": {"k_a": 0.04, "k_b": 0.0},
        "use_effective_ws": True,
        "ceps": 0.2,
    },
    "axial_induction_model": "Madsen",
    "deflection_model": {"name": "None"},
    "turbulence_model": {"name": "None"},
    "superposition_model": {
        "ws_superposition": "Product",
        "ti_superposition": "Squared",
    },
    "rotor_averaging": {
        "name": "none",
        "background_averaging": "center",
        "wake_averaging": "center",
        "wind_speed_exponent_for_power": 3,
        "wind_speed_exponent_for_ct": 2,
    },
    "blockage_model": {"name": "None"},
    "wm_coupling": {"method": "PB", "wake_tool": "foxes"},
    "apm_grid": {"Lx": 1.0e5, "Ly": 1.0e5, "dx": 1.0e3},
}


def test_xy_points_matches_xy_plane_foxes():
    """Foxes twin of the _xy_plane_points drift guard. The solve itself also
    exercises the background_flow_direction shim: without it, wayve@e87780a's
    broken e_spanwise import crashes every foxes-coupled solve."""
    pytest.importorskip("foxes")
    import numpy as np
    from wayve.apm import APM
    from wayve.grid.grid import Stat2Dgrid
    from wayve.momentum_flux_parametrizations import FrictionCoefficients
    from wayve.solvers import FixedPointIteration

    from wifa.wayve_api import (
        _pressure_for_state,
        _xy_plane_points,
        flow_io_abl,
        wf_setup,
    )

    system = _solver_system(
        [0.0, 400.0, 800.0], [0.0, 0.0, 0.0], _scalar_wind_resource(270.0)
    )
    farm_dat = system["wind_farm"]
    resource = system["site"]["energy_resource"]["wind_resource"]
    abl, rotation = flow_io_abl(resource, 0, zh=80.0, h1=160.0)
    wind_farm, forcing, _, _ = wf_setup(
        farm_dat, _FOXES_ANALYSIS, 1.0e3, rotation=rotation
    )
    wake_model = wind_farm.coupling.wake_model
    assert type(wake_model).__name__ == "FoxesWakeModel"
    grid = Stat2Dgrid(1.0e5, 100, 1.0e5, 100)
    model = APM(grid, forcing, abl, FrictionCoefficients(), _pressure_for_state(abl, 1))
    model.solve(method=FixedPointIteration(tol=5.0e-3, relax=0.7))
    coupling = wind_farm.coupling
    u_bg_evaluator = coupling.set_up_u_bg_evaluator(abl)
    apm_evaluator = coupling.apm_evaluator
    xs = np.linspace(-2000.0, 4000.0, 13)
    ys = np.linspace(-1500.0, 1500.0, 11)
    ref = wake_model.xy_plane(
        wind_farm, abl, u_bg_evaluator, apm_evaluator, xs, ys, 80.0
    )
    Xs, Ys = np.meshgrid(xs, ys, indexing="ij")
    new = _xy_plane_points(
        wake_model, wind_farm, abl, u_bg_evaluator, apm_evaluator, Xs, Ys, 80.0
    )
    for ref_field, new_field in zip(ref, new):
        np.testing.assert_array_equal(new_field, ref_field)


def test_turbulence_profile_stress_rotation():
    """Explicit tau_x/tau_y stress profiles are earth-frame vectors: they must
    rotate into the solver frame with the velocities (magnitude preserved)."""
    import numpy as np

    from wifa.wayve_api import flow_io_abl

    zs = np.linspace(10.0, 1700.0, 18)
    # Jet-like profile: wayve's ABL estimates the eddy viscosity from where
    # the shear turns negative, so the speed must peak below the profile top
    ws = 8.0 + 3.0 * np.minimum(zs, 1100.0) / 1100.0
    ws -= 0.5 * np.maximum(0.0, zs - 1100.0) / 600.0
    wd = np.full(zs.size, 200.0)
    theta = 290.0 + 3.5e-3 * zs
    # Stress magnitude decaying to zero at 800 m, at a fixed earth angle
    tau_mag = 0.12 * np.maximum(0.0, 1.0 - zs / 800.0) ** 1.5
    phi_tau = np.deg2rad(30.0)
    tau_x_in = tau_mag * np.cos(phi_tau)
    tau_y_in = tau_mag * np.sin(phi_tau)
    scalar = lambda val: {"data": [val]}  # noqa: E731
    resource = {
        "time": [0],
        "height": list(zs),
        "wind_speed": {"data": [list(ws)]},
        "wind_direction": {"data": [list(wd)]},
        "potential_temperature": {"data": [list(theta)]},
        "tau_x": {"data": [list(tau_x_in)]},
        "tau_y": {"data": [list(tau_y_in)]},
        "turbulence_intensity": scalar(0.08),
        "z0": scalar(0.03),
        "fc": scalar(1.0e-4),
        "thermal_stratification": {
            "capping_inversion": {
                "ABL_height": scalar(1100.0),
                "dH": scalar(200.0),
                "dtheta": scalar(3.0),
                "lapse_rate": scalar(3.5e-3),
            }
        },
    }
    abl, rotation = flow_io_abl(resource, 0, zh=80.0, h1=160.0)
    assert rotation == pytest.approx(np.deg2rad(270.0 - 200.0))
    c, s = np.cos(rotation), np.sin(rotation)
    np.testing.assert_allclose(
        abl.tauxs, c * tau_x_in + s * tau_y_in, rtol=1e-12, atol=1e-15
    )
    np.testing.assert_allclose(
        abl.tauys, c * tau_y_in - s * tau_x_in, rtol=1e-12, atol=1e-15
    )
    np.testing.assert_allclose(
        np.hypot(abl.tauxs, abl.tauys), tau_mag, rtol=1e-12, atol=1e-15
    )
