"""Tests for the Eddy Viscosity (EV) bundle dispatch in wifa/pywake_api.py.

The EV model exists only on pyWake's unmerged ``cj_add_eddy_viscosity_model``
branch (mirrored at eu-flow/wp5/PyWake@flow-model-chain), so this file lives
apart from test_pywake_submodels.py and every test is skipped on a released
pyWake.
"""

import pytest

eddy_viscosity = pytest.importorskip(
    "py_wake.deficit_models.eddy_viscosity",
    reason="EV requires pyWake's cj_add_eddy_viscosity_model branch",
)

from py_wake.rotor_avg_models.simplified_gaussian_rotor_average_model import (  # noqa: E402
    SimplifiedGaussianRotorAverageModel,
)
from py_wake.turbulence_models.quarton_and_ainslie import (  # noqa: E402
    ModifiedQuartonAndAinslieTurbulenceModel,
)

from wifa.pywake_api import (  # noqa: E402
    _configure_deficit_model,
    _configure_rotor_averaging,
    _configure_turbulence_model,
)

_RD = 126.0
_HH = 90.0


def _call_deficit(name, analysis_extra=None):
    wind_deficit_model = {"name": name, **(analysis_extra or {})}
    analysis = {"wind_deficit_model": wind_deficit_model}
    cls, args, _post = _configure_deficit_model({"name": name}, analysis, _RD, _HH)
    return cls, args


@pytest.mark.parametrize(
    "name", ["EddyViscosity", "eddyviscosity", "EDDYVISCOSITY", "eddy_viscosity"]
)
def test_configure_deficit_model_eddy_viscosity(name):
    cls, args = _call_deficit(name)
    assert cls is eddy_viscosity.EddyViscosityDeficitModel
    # EV pairs with MaxSum: deficits are referenced to the free stream, so the
    # usual use_effective_ws=True default flips to False for this model.
    assert args["use_effective_ws"] is False


def test_eddy_viscosity_use_effective_ws_override():
    _cls, args = _call_deficit("EddyViscosity", {"use_effective_ws": True})
    assert args["use_effective_ws"] is True


def test_eddy_viscosity_free_stream_ti_mapped():
    # EV is TI-capable: windIO's free_stream_ti flag maps to use_effective_ti.
    _cls, args = _call_deficit(
        "EddyViscosity", {"wake_expansion_coefficient": {"free_stream_ti": True}}
    )
    assert args["use_effective_ti"] is False


def test_eddy_viscosity_no_ct2a_injected():
    # EddyViscosityDeficitModel has no ct2a parameter; the axial-induction
    # mapping must skip it rather than pass an unexpected kwarg.
    wind_deficit_model = {"name": "EddyViscosity"}
    analysis = {
        "wind_deficit_model": wind_deficit_model,
        "axial_induction_model": "1D",
    }
    _cls, args, _post = _configure_deficit_model(
        {"name": "EddyViscosity"}, analysis, _RD, _HH
    )
    assert "ct2a" not in args


@pytest.mark.parametrize(
    "name",
    [
        "QuartonAndAinslie",
        "quarton_and_ainslie",
        "ModifiedQuartonAndAinslie",
        "modifiedquartonandainslie",
    ],
)
def test_configure_turbulence_model_quarton_and_ainslie(name):
    model = _configure_turbulence_model({"name": name})
    assert isinstance(model, ModifiedQuartonAndAinslieTurbulenceModel)


@pytest.mark.parametrize("name", ["simplified_gaussian", "SimplifiedGaussian"])
def test_configure_rotor_averaging_simplified_gaussian(name):
    model = _configure_rotor_averaging({"name": name})
    assert isinstance(model, SimplifiedGaussianRotorAverageModel)


def test_eddy_viscosity_deficit_constructs():
    # The dispatched class + args must actually build (triggers the ~45 ms
    # LUT generation) with the EV bundle defaults.
    cls, args = _call_deficit("EddyViscosity")
    deficit = cls(**args)
    assert deficit.use_effective_ws is False
    assert deficit.use_effective_ti is True
