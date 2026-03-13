"""Parametrized unit tests for PyWake submodel configuration functions.

Tests each _configure_*() function in wifa/pywake_api.py to verify:
- Correct PyWake class is returned for each model name
- Parameters are passed through correctly
- NotImplementedError raised for unsupported names
- Case-insensitive matching works
"""

import pytest
from py_wake.deficit_models import (
    HybridInduction,
    RankineHalfBody,
    SelfSimilarityDeficit,
    SelfSimilarityDeficit2020,
    VortexCylinder,
    VortexDipole,
)
from py_wake.deficit_models.gaussian import (
    BastankhahGaussianDeficit,
    BlondelSuperGaussianDeficit2020,
    BlondelSuperGaussianDeficit2023,
    CarbajofuertesGaussianDeficit,
    NiayifarGaussianDeficit,
    TurboGaussianDeficit,
    ZongGaussianDeficit,
)
from py_wake.deficit_models.gcl import GCLDeficit
from py_wake.deficit_models.noj import NOJLocalDeficit, TurboNOJDeficit
from py_wake.deficit_models.rathmann import Rathmann
from py_wake.deflection_models import JimenezWakeDeflection
from py_wake.deflection_models.gcl_hill_vortex import GCLHillDeflection
from py_wake.rotor_avg_models import (
    CGIRotorAvg,
    EqGridRotorAvg,
    GQGridRotorAvg,
    GridRotorAvg,
    PolarGridRotorAvg,
    RotorCenter,
)
from py_wake.superposition_models import (
    CumulativeWakeSum,
    LinearSum,
    MaxSum,
    SquaredSum,
    WeightedSum,
)
from py_wake.turbulence_models import (
    CrespoHernandez,
    STF2005TurbulenceModel,
    STF2017TurbulenceModel,
)
from py_wake.turbulence_models.gcl_turb import GCLTurbulence

from wifa.pywake_api import (
    DEFAULTS,
    _configure_blockage_model,
    _configure_deficit_model,
    _configure_deflection_model,
    _configure_rotor_averaging,
    _configure_superposition_model,
    _configure_turbulence_model,
    configure_wake_model,
    get_with_default,
)

# Default rotor diameter and hub height for deficit model tests
_RD = 126.0
_HH = 90.0


def _call_deficit(name, analysis_extra=None):
    """Helper to call _configure_deficit_model with minimal boilerplate."""
    wind_deficit_model = {"name": name, **(analysis_extra or {})}
    analysis = {"wind_deficit_model": wind_deficit_model}
    return _configure_deficit_model({"name": name}, analysis, _RD, _HH)


# ---------------------------------------------------------------------------
# Deficit model tests
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "name,expected_class",
    [
        ("Jensen", NOJLocalDeficit),
        ("jensen", NOJLocalDeficit),
        ("JENSEN", NOJLocalDeficit),
        ("Bastankhah2014", BastankhahGaussianDeficit),
        ("bastankhah2014", BastankhahGaussianDeficit),
        ("BASTANKHAH2014", BastankhahGaussianDeficit),
        ("SuperGaussian", BlondelSuperGaussianDeficit2020),
        ("supergaussian", BlondelSuperGaussianDeficit2020),
        ("SuperGaussian2023", BlondelSuperGaussianDeficit2023),
        ("TurboPark", TurboGaussianDeficit),
        ("TurbOPark", TurboGaussianDeficit),
        ("turbopark", TurboGaussianDeficit),
        ("Niayifar2016", NiayifarGaussianDeficit),
        ("niayifar2016", NiayifarGaussianDeficit),
        ("Zong2020", ZongGaussianDeficit),
        ("zong2020", ZongGaussianDeficit),
        ("Carbajofuertes2018", CarbajofuertesGaussianDeficit),
        ("carbajofuertes2018", CarbajofuertesGaussianDeficit),
        ("TurboNOJ", TurboNOJDeficit),
        ("turbonoj", TurboNOJDeficit),
        ("GCL", GCLDeficit),
        ("gcl", GCLDeficit),
        ("NOJLocalDeficit", NOJLocalDeficit),
        ("nojlocaldeficit", NOJLocalDeficit),
        ("NOJLOCALDEFICIT", NOJLocalDeficit),
    ],
)
def test_configure_deficit_model(name, expected_class):
    cls, args = _call_deficit(name)
    assert cls is expected_class


@pytest.mark.parametrize(
    "name,expected_class",
    [
        ("Jensen", NOJLocalDeficit),
        ("Bastankhah2014", BastankhahGaussianDeficit),
        ("SuperGaussian", BlondelSuperGaussianDeficit2020),
        ("SuperGaussian2023", BlondelSuperGaussianDeficit2023),
        ("TurboPark", TurboGaussianDeficit),
        ("Niayifar2016", NiayifarGaussianDeficit),
        ("Zong2020", ZongGaussianDeficit),
        ("Carbajofuertes2018", CarbajofuertesGaussianDeficit),
        ("TurboNOJ", TurboNOJDeficit),
        ("GCL", GCLDeficit),
        ("NOJLocalDeficit", NOJLocalDeficit),
    ],
)
def test_configure_deficit_model_instantiation(name, expected_class):
    """Verify returned kwargs can actually instantiate the model without TypeError."""
    cls, args = _call_deficit(name)
    instance = cls(**args)
    assert isinstance(instance, expected_class)


def test_configure_deficit_model_bastankhah2014_params():
    """Verify wake expansion and ceps params are passed for Bastankhah2014."""
    cls, args = _call_deficit(
        "Bastankhah2014",
        {"wake_expansion_coefficient": {"k": 0.04}, "ceps": 0.2},
    )
    assert cls is BastankhahGaussianDeficit
    assert args["k"] == 0.04
    assert args["ceps"] == 0.2


def test_configure_deficit_model_bastankhah2014_k_b():
    """Verify k_b wake expansion is used when present."""
    _, args = _call_deficit(
        "Bastankhah2014",
        {"wake_expansion_coefficient": {"k_a": 0.38, "k_b": 0.004}},
    )
    assert args["k"] == 0.004


def test_configure_deficit_model_jensen_k_b():
    """Verify Jensen k_a/k_b expansion params."""
    _, args = _call_deficit(
        "Jensen",
        {"wake_expansion_coefficient": {"k_a": 0.38, "k_b": 0.004}},
    )
    assert args["a"] == [0.38, 0.004]


def test_configure_deficit_model_nojlocaldeficit_k_b():
    """Verify NOJLocalDeficit k_a/k_b expansion params."""
    _, args = _call_deficit(
        "NOJLocalDeficit",
        {"wake_expansion_coefficient": {"k_a": 0.38, "k_b": 0.004}},
    )
    assert args["a"] == [0.38, 0.004]


def test_configure_deficit_model_gaussian_params_niayifar():
    """Verify Gaussian params pass through for Niayifar2016."""
    cls, args = _call_deficit(
        "Niayifar2016",
        {
            "wake_expansion_coefficient": {"k_a": 0.38, "k_b": 0.004},
            "ceps": 0.3,
            "use_effective_ti": True,
        },
    )
    assert cls is NiayifarGaussianDeficit
    assert args["a"] == [0.38, 0.004]
    assert args["ceps"] == 0.3
    assert args["use_effective_ti"] is True


def test_configure_deficit_model_zong_no_ceps():
    """Verify Zong2020 does not pass ceps (unsupported)."""
    _, args = _call_deficit(
        "Zong2020",
        {"ceps": 0.3, "use_effective_ti": False},
    )
    assert "ceps" not in args
    assert args["use_effective_ti"] is False


def test_configure_deficit_model_bastankhah2014_no_effective_ti():
    """Verify Bastankhah2014 does not pass use_effective_ti (unsupported)."""
    _, args = _call_deficit(
        "Bastankhah2014",
        {"wake_expansion_coefficient": {"k": 0.04}, "use_effective_ti": True},
    )
    assert "use_effective_ti" not in args
    assert args["k"] == 0.04


def test_configure_deficit_model_a_param_warns_on_scalar_k():
    """Verify warning when scalar k is provided for a=[k_a, k_b] models."""
    with pytest.warns(UserWarning, match="uses a="):
        _, args = _call_deficit(
            "Niayifar2016", {"wake_expansion_coefficient": {"k": 0.05}}
        )
    assert "k" not in args
    assert "a" not in args


def test_configure_deficit_model_a_param_warns_on_missing_k_a():
    """Verify warning when k_b is provided without k_a."""
    with pytest.warns(UserWarning, match="k_a not specified"):
        _, args = _call_deficit(
            "Zong2020", {"wake_expansion_coefficient": {"k_b": 0.004}}
        )
    assert args["a"] == [0, 0.004]


@pytest.mark.parametrize(
    "name,extra,expected_class",
    [
        (
            "Bastankhah2014",
            {"wake_expansion_coefficient": {"k": 0.04}, "ceps": 0.2},
            BastankhahGaussianDeficit,
        ),
        (
            "Niayifar2016",
            {
                "wake_expansion_coefficient": {"k_a": 0.38, "k_b": 0.004},
                "ceps": 0.3,
                "use_effective_ti": True,
            },
            NiayifarGaussianDeficit,
        ),
        ("Zong2020", {"use_effective_ti": False}, ZongGaussianDeficit),
        (
            "Jensen",
            {"wake_expansion_coefficient": {"k_a": 0.38, "k_b": 0.004}},
            NOJLocalDeficit,
        ),
        (
            "NOJLocalDeficit",
            {"wake_expansion_coefficient": {"k_a": 0.38, "k_b": 0.004}},
            NOJLocalDeficit,
        ),
    ],
)
def test_configure_deficit_model_instantiation_with_params(name, extra, expected_class):
    """Verify models with user-specified params can be instantiated."""
    cls, args = _call_deficit(name, extra)
    assert isinstance(cls(**args), expected_class)


def test_configure_deficit_model_turbonoj_A_param():
    """Verify TurboNOJ passes through the A parameter and can be instantiated."""
    cls, args = _call_deficit("TurboNOJ", {"A": 0.6})
    assert cls is TurboNOJDeficit
    assert args["A"] == 0.6
    instance = cls(**args)
    assert isinstance(instance, TurboNOJDeficit)


@pytest.mark.parametrize("name", ["Bastankhah2016", "bastankhah2016"])
def test_configure_deficit_model_bastankhah2016_not_implemented(name):
    with pytest.raises(NotImplementedError, match="Bastankhah2016"):
        _call_deficit(name)


def test_configure_deficit_model_unknown():
    with pytest.raises(NotImplementedError, match="NonexistentModel"):
        _call_deficit("NonexistentModel")


# ---------------------------------------------------------------------------
# Deflection model tests
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "name,expected_class",
    [
        ("Jimenez", JimenezWakeDeflection),
        ("jimenez", JimenezWakeDeflection),
        ("JIMENEZ", JimenezWakeDeflection),
        ("GCLHill", GCLHillDeflection),
        ("gclhill", GCLHillDeflection),
        ("GCLhill", GCLHillDeflection),
    ],
)
def test_configure_deflection_model(name, expected_class):
    model = _configure_deflection_model({"name": name, "beta": 0.1})
    assert isinstance(model, expected_class)


@pytest.mark.parametrize("name", [None, "None", "none", "NONE"])
def test_configure_deflection_model_none(name):
    assert _configure_deflection_model({"name": name, "beta": 0.1}) is None


def test_configure_deflection_model_bastankhah2016():
    with pytest.raises(NotImplementedError, match="Bastankhah2016"):
        _configure_deflection_model({"name": "Bastankhah2016", "beta": 0.1})


def test_configure_deflection_model_unknown():
    with pytest.raises(NotImplementedError, match="UnknownDeflection"):
        _configure_deflection_model({"name": "UnknownDeflection", "beta": 0.1})


# ---------------------------------------------------------------------------
# Turbulence model tests
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "name,expected_class",
    [
        ("STF2005", STF2005TurbulenceModel),
        ("stf2005", STF2005TurbulenceModel),
        ("STF2017", STF2017TurbulenceModel),
        ("stf2017", STF2017TurbulenceModel),
        ("CrespoHernandez", CrespoHernandez),
        ("crespohernandez", CrespoHernandez),
        ("CRESPOHERNANDEZ", CrespoHernandez),
        ("IEC-TI-2019", STF2017TurbulenceModel),
        ("iec-ti-2019", STF2017TurbulenceModel),
        ("GCL", GCLTurbulence),
        ("gcl", GCLTurbulence),
    ],
)
def test_configure_turbulence_model(name, expected_class):
    data = {"name": name, "c1": 1.0, "c2": 1.0}
    assert isinstance(_configure_turbulence_model(data), expected_class)


@pytest.mark.parametrize("name", [None, "None", "none", "NONE"])
def test_configure_turbulence_model_none(name):
    assert _configure_turbulence_model({"name": name, "c1": 1.0, "c2": 1.0}) is None


def test_configure_turbulence_model_unknown():
    with pytest.raises(NotImplementedError, match="UnknownTurb"):
        _configure_turbulence_model({"name": "UnknownTurb", "c1": 1.0, "c2": 1.0})


# ---------------------------------------------------------------------------
# Superposition model tests
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "name,expected_class",
    [
        ("Linear", LinearSum),
        ("linear", LinearSum),
        ("LINEAR", LinearSum),
        ("Squared", SquaredSum),
        ("squared", SquaredSum),
        ("Max", MaxSum),
        ("max", MaxSum),
        ("Weighted", WeightedSum),
        ("weighted", WeightedSum),
        ("Cumulative", CumulativeWakeSum),
        ("cumulative", CumulativeWakeSum),
    ],
)
def test_configure_superposition_model(name, expected_class):
    assert isinstance(
        _configure_superposition_model({"ws_superposition": name}), expected_class
    )


def test_configure_superposition_model_product_not_implemented():
    with pytest.raises(NotImplementedError, match="Product"):
        _configure_superposition_model({"ws_superposition": "Product"})


def test_configure_superposition_model_unknown():
    with pytest.raises(NotImplementedError, match="UnknownSuper"):
        _configure_superposition_model({"ws_superposition": "UnknownSuper"})


# ---------------------------------------------------------------------------
# Rotor averaging tests
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "name,expected_class",
    [
        ("Center", RotorCenter),
        ("center", RotorCenter),
        ("CENTER", RotorCenter),
        ("avg_deficit", GridRotorAvg),
        ("Avg_Deficit", GridRotorAvg),
        ("EqGrid", EqGridRotorAvg),
        ("eqgrid", EqGridRotorAvg),
        ("GQGrid", GQGridRotorAvg),
        ("gqgrid", GQGridRotorAvg),
        ("PolarGrid", PolarGridRotorAvg),
        ("polargrid", PolarGridRotorAvg),
        ("CGI", CGIRotorAvg),
        ("cgi", CGIRotorAvg),
    ],
)
def test_configure_rotor_averaging(name, expected_class):
    assert isinstance(_configure_rotor_averaging({"name": name}), expected_class)


def test_configure_rotor_averaging_eqgrid_n():
    assert isinstance(
        _configure_rotor_averaging({"name": "EqGrid", "n": 9}), EqGridRotorAvg
    )


def test_configure_rotor_averaging_gqgrid_params():
    data = {"name": "GQGrid", "n_x_grid_points": 3, "n_y_grid_points": 5}
    model = _configure_rotor_averaging(data)
    assert isinstance(model, GQGridRotorAvg)
    # Verify custom params produced different nodes than defaults
    default = GQGridRotorAvg()
    assert len(model.nodes_x) != len(default.nodes_x)


def test_configure_rotor_averaging_cgi_n():
    assert isinstance(_configure_rotor_averaging({"name": "CGI", "n": 7}), CGIRotorAvg)


def test_configure_rotor_averaging_unknown():
    with pytest.raises(NotImplementedError, match="UnknownRotor"):
        _configure_rotor_averaging({"name": "UnknownRotor"})


# ---------------------------------------------------------------------------
# Blockage model tests
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "name,expected_class",
    [
        ("SelfSimilarityDeficit2020", SelfSimilarityDeficit2020),
        ("selfsimilaritydeficit2020", SelfSimilarityDeficit2020),
        ("SelfSimilarityDeficit", SelfSimilarityDeficit),
        ("selfsimilaritydeficit", SelfSimilarityDeficit),
        ("RankineHalfBody", RankineHalfBody),
        ("rankinehalfbody", RankineHalfBody),
        ("Rathmann", Rathmann),
        ("rathmann", Rathmann),
        ("VortexCylinder", VortexCylinder),
        ("vortexcylinder", VortexCylinder),
        ("VortexDipole", VortexDipole),
        ("vortexdipole", VortexDipole),
        ("HybridInduction", HybridInduction),
        ("hybridinduction", HybridInduction),
    ],
)
def test_configure_blockage_model(name, expected_class):
    model = _configure_blockage_model({"name": name, "ss_alpha": 0.888}, {})
    assert isinstance(model, expected_class)


@pytest.mark.parametrize("name", [None, "None", "none", "NONE"])
def test_configure_blockage_model_none(name):
    assert _configure_blockage_model({"name": name}, {}) is None


def test_configure_blockage_model_unknown():
    with pytest.raises(NotImplementedError, match="UnknownBlockage"):
        _configure_blockage_model({"name": "UnknownBlockage"}, {})


# ---------------------------------------------------------------------------
# get_with_default preserves extra user keys
# ---------------------------------------------------------------------------


def test_get_with_default_preserves_extra_keys():
    """Verify that get_with_default merges defaults without dropping user keys."""
    analysis = {
        "rotor_averaging": {
            "name": "GQGrid",
            "n_x_grid_points": 3,
            "n_y_grid_points": 5,
        },
    }
    result = get_with_default(analysis, "rotor_averaging", DEFAULTS)
    assert result["name"] == "GQGrid"
    assert result["n_x_grid_points"] == 3
    assert result["n_y_grid_points"] == 5


def test_get_with_default_rotor_avg_eqgrid_n():
    """Verify EqGrid 'n' param survives through get_with_default."""
    analysis = {"rotor_averaging": {"name": "EqGrid", "n": 9}}
    result = get_with_default(analysis, "rotor_averaging", DEFAULTS)
    model = _configure_rotor_averaging(result)
    assert isinstance(model, EqGridRotorAvg)
    assert result["n"] == 9


def test_get_with_default_fills_missing_keys():
    """Verify that missing keys are filled from defaults."""
    # deflection_model defaults have beta=0.1; user only provides name
    analysis = {"deflection_model": {"name": "Jimenez"}}
    result = get_with_default(analysis, "deflection_model", DEFAULTS)
    assert result["name"] == "Jimenez"
    assert result["beta"] == 0.1


def test_get_with_default_recursive_nested_dicts():
    """Verify recursive merge fills deep missing keys while preserving user extras."""
    nested_defaults = {
        "model": {
            "params": {"a": 1, "b": 2},
            "name": "default",
        }
    }
    # User provides partial nested dict (missing key "b") plus an extra key "c"
    data = {
        "model": {
            "params": {"a": 10, "c": 99},
            "name": "custom",
        }
    }
    result = get_with_default(data, "model", nested_defaults)
    assert result["name"] == "custom"
    assert result["params"]["a"] == 10  # user value preserved
    assert result["params"]["b"] == 2  # missing key filled from defaults
    assert result["params"]["c"] == 99  # extra user key preserved


# ---------------------------------------------------------------------------
# configure_wake_model return contract
# ---------------------------------------------------------------------------


def test_configure_wake_model_returns_wake_deficit_key():
    """Verify configure_wake_model returns wake_deficit_key for API compat."""
    system_dat = {
        "attributes": {
            "analysis": {
                "wind_deficit_model": {"name": "Jensen"},
            }
        }
    }
    config = configure_wake_model(system_dat, rotor_diameter=126.0, hub_height=90.0)
    assert "wake_deficit_key" in config
    assert config["wake_deficit_key"] is None
