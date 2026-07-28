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
from py_wake.deficit_models.noj import NOJDeficit, NOJLocalDeficit, TurboNOJDeficit
from py_wake.deficit_models.rathmann import Rathmann
from py_wake.deficit_models.utils import ct2a_madsen, ct2a_mom1d
from py_wake.deflection_models import JimenezWakeDeflection
from py_wake.deflection_models.gcl_hill_vortex import GCLHillDeflection
from py_wake.ground_models.ground_models import Mirror
from py_wake.rotor_avg_models import (
    AreaOverlapAvgModel,
    CGIRotorAvg,
    EqGridRotorAvg,
    GaussianOverlapAvgModel,
    GQGridRotorAvg,
    GridRotorAvg,
    PolarGridRotorAvg,
    RotorCenter,
)
from py_wake.superposition_models import (
    CumulativeWakeSum,
    LinearSum,
    MaxSum,
    SqrMaxSum,
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
    _fuga_atmosphere,
    _fuga_z0_sweep,
    configure_wake_model,
    get_with_default,
)

# Default rotor diameter and hub height for deficit model tests
_RD = 126.0
_HH = 90.0


def _call_deficit(name, analysis_extra=None, analysis_top=None):
    """Helper to call _configure_deficit_model with minimal boilerplate.

    Returns ``(wake_model_class, deficit_args)``; the third element
    (``deficit_post_attrs``) is dropped for the common case.  ``analysis_top``
    merges keys at the ``analysis`` level (e.g. ``axial_induction_model``).
    """
    wind_deficit_model = {"name": name, **(analysis_extra or {})}
    analysis = {"wind_deficit_model": wind_deficit_model, **(analysis_top or {})}
    cls, args, _post = _configure_deficit_model({"name": name}, analysis, _RD, _HH)
    return cls, args


def _call_deficit_full(name, analysis_extra=None, analysis_top=None):
    """Like ``_call_deficit`` but returns the full 3-tuple including post-attrs."""
    wind_deficit_model = {"name": name, **(analysis_extra or {})}
    analysis = {"wind_deficit_model": wind_deficit_model, **(analysis_top or {})}
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
        ("Jensen_1983", NOJDeficit),
        ("jensen_1983", NOJDeficit),
        ("JENSEN_1983", NOJDeficit),
        ("NOJDeficit", NOJDeficit),
        ("nojdeficit", NOJDeficit),
        ("NOJDEFICIT", NOJDeficit),
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
        ("Jensen_1983", NOJDeficit),
        ("NOJDeficit", NOJDeficit),
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


def test_configure_deficit_model_bastankhah2014_k_a():
    """The windIO constant k_a is Bastankhah2014's scalar k; a nonzero TI
    coefficient k_b cannot be represented and warns."""
    with pytest.warns(UserWarning, match="k_b=0.38 is ignored"):
        _, args = _call_deficit(
            "Bastankhah2014",
            {"wake_expansion_coefficient": {"k_a": 0.004, "k_b": 0.38}},
        )
    assert args["k"] == 0.004


def test_configure_deficit_model_jensen_k_a_k_b():
    """Jensen: windIO k = k_a + k_b*TI maps to PyWake a = [k_b, k_a]."""
    _, args = _call_deficit(
        "Jensen",
        {"wake_expansion_coefficient": {"k_a": 0.004, "k_b": 0.38}},
    )
    assert args["a"] == [0.38, 0.004]


def test_configure_deficit_model_nojlocaldeficit_k_a_k_b():
    """NOJLocalDeficit: windIO k = k_a + k_b*TI maps to PyWake a = [k_b, k_a]."""
    _, args = _call_deficit(
        "NOJLocalDeficit",
        {"wake_expansion_coefficient": {"k_a": 0.004, "k_b": 0.38}},
    )
    assert args["a"] == [0.38, 0.004]


def test_configure_deficit_model_jensen_1983_k():
    """Verify Jensen_1983 passes scalar k and does not pass use_effective_ws."""
    cls, args = _call_deficit(
        "Jensen_1983",
        {"wake_expansion_coefficient": {"k": 0.04}},
    )
    assert cls is NOJDeficit
    assert args["k"] == 0.04
    assert "use_effective_ws" not in args


def test_configure_deficit_model_jensen_1983_k_a():
    """Jensen_1983 (NOJDeficit) takes its scalar k from windIO's constant k_a,
    since the wake_expansion_coefficient schema has no scalar k field."""
    cls, args = _call_deficit(
        "Jensen_1983",
        {"wake_expansion_coefficient": {"k_a": 0.1, "k_b": 0.0}},
    )
    assert cls is NOJDeficit
    assert args["k"] == 0.1


def test_configure_deficit_model_gaussian_params_niayifar():
    """Verify Gaussian params pass through for Niayifar2016."""
    cls, args = _call_deficit(
        "Niayifar2016",
        {
            "wake_expansion_coefficient": {
                "k_a": 0.004,
                "k_b": 0.38,
                "free_stream_ti": False,
            },
            "ceps": 0.3,
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
        {"ceps": 0.3, "wake_expansion_coefficient": {"free_stream_ti": True}},
    )
    assert "ceps" not in args
    assert args["use_effective_ti"] is False


def test_configure_deficit_model_bastankhah2014_no_effective_ti():
    """Verify Bastankhah2014 does not pass use_effective_ti (unsupported)."""
    _, args = _call_deficit(
        "Bastankhah2014",
        {"wake_expansion_coefficient": {"k": 0.04, "free_stream_ti": False}},
    )
    assert "use_effective_ti" not in args
    assert args["k"] == 0.04


def test_configure_deficit_model_a_param_warns_on_scalar_k():
    """Verify warning when scalar k is provided for k_a/k_b models."""
    with pytest.warns(UserWarning, match="uses k_a/k_b"):
        _, args = _call_deficit(
            "Niayifar2016", {"wake_expansion_coefficient": {"k": 0.05}}
        )
    assert "k" not in args
    assert "a" not in args


def test_configure_deficit_model_a_param_warns_on_missing_k_b():
    """Verify warning when only the constant k_a is given to a TI-dependent
    expansion model (k_b defaults to 0)."""
    with pytest.warns(UserWarning, match="k_b not specified"):
        _, args = _call_deficit(
            "Zong2020", {"wake_expansion_coefficient": {"k_a": 0.004}}
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
                "wake_expansion_coefficient": {
                    "k_a": 0.004,
                    "k_b": 0.38,
                    "free_stream_ti": False,
                },
                "ceps": 0.3,
            },
            NiayifarGaussianDeficit,
        ),
        (
            "Zong2020",
            {
                "wake_expansion_coefficient": {
                    "k_a": 0.004,
                    "k_b": 0.38,
                    "free_stream_ti": True,
                }
            },
            ZongGaussianDeficit,
        ),
        (
            "Jensen",
            {"wake_expansion_coefficient": {"k_a": 0.004, "k_b": 0.38}},
            NOJLocalDeficit,
        ),
        (
            "NOJLocalDeficit",
            {"wake_expansion_coefficient": {"k_a": 0.004, "k_b": 0.38}},
            NOJLocalDeficit,
        ),
        (
            "Jensen_1983",
            {"wake_expansion_coefficient": {"k": 0.04}},
            NOJDeficit,
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
# TI reference flag (windIO free_stream_ti -> PyWake use_effective_ti)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "name",
    [
        "Jensen",
        "Niayifar2016",
        "Carbajofuertes2018",
        "Zong2020",
        "TurbOPark",
        "SuperGaussian",
        "SuperGaussian2023",
    ],
)
@pytest.mark.parametrize("free_stream_ti,expected", [(False, True), (True, False)])
def test_free_stream_ti_inverts_to_use_effective_ti(name, free_stream_ti, expected):
    """free_stream_ti maps to use_effective_ti with inverted polarity, for every
    TI-capable deficit (including SuperGaussian and TurbOPark, which the previous
    narrow handling missed)."""
    _, args = _call_deficit(
        name, {"wake_expansion_coefficient": {"free_stream_ti": free_stream_ti}}
    )
    assert args["use_effective_ti"] is expected
    # kwargs must actually instantiate the model
    cls, _ = _call_deficit(name)
    cls(**args)


@pytest.mark.parametrize("name", ["Bastankhah2014"])
def test_free_stream_ti_ignored_for_non_ti_capable(name):
    """Deficits without a use_effective_ti param must not receive it, even if
    free_stream_ti is present (would raise TypeError on instantiation).

    (GCL was moved to TI_CAPABLE — GCLDeficit accepts use_effective_ti, as
    GCLLocal demonstrates.)"""
    _, args = _call_deficit(
        name, {"wake_expansion_coefficient": {"k_b": 0.04, "free_stream_ti": True}}
    )
    assert "use_effective_ti" not in args


def test_use_effective_ti_key_is_ignored():
    """The deprecated top-level use_effective_ti key is no longer read."""
    _, args = _call_deficit(
        "Niayifar2016",
        {
            "wake_expansion_coefficient": {"k_a": 0.38, "k_b": 0.004},
            "use_effective_ti": False,
        },
    )
    assert "use_effective_ti" not in args


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


def test_configure_superposition_model_vector_not_implemented():
    """Vector superposition is foxes-only; the pyWake path rejects it."""
    with pytest.raises(NotImplementedError, match="Vector"):
        _configure_superposition_model({"ws_superposition": "Vector"})


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
        ("grid", GridRotorAvg),
        ("avg_deficit", GridRotorAvg),
        ("Avg_Deficit", GridRotorAvg),
        ("gaussian_overlap", GaussianOverlapAvgModel),
        ("gaussianoverlap", GaussianOverlapAvgModel),
        ("area_overlap", AreaOverlapAvgModel),
        ("areaoverlap", AreaOverlapAvgModel),
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


@pytest.mark.parametrize(
    "name,expected_class",
    [
        ("SelfSimilarityDeficit2020", SelfSimilarityDeficit2020),
        ("SelfSimilarityDeficit", SelfSimilarityDeficit),
        ("RankineHalfBody", RankineHalfBody),
        ("Rathmann", Rathmann),
        ("VortexCylinder", VortexCylinder),
        ("VortexDipole", VortexDipole),
        ("HybridInduction", HybridInduction),
    ],
)
def test_configure_blockage_model_ground_mirror(name, expected_class):
    model = _configure_blockage_model({"name": name, "ground_mirror": True}, {})
    assert isinstance(model, expected_class)
    assert isinstance(model.groundModel, Mirror)


def test_configure_blockage_model_no_ground_mirror_by_default():
    model = _configure_blockage_model({"name": "SelfSimilarityDeficit2020"}, {})
    assert model.groundModel is None

    model = _configure_blockage_model(
        {"name": "SelfSimilarityDeficit2020", "ground_mirror": False}, {}
    )
    assert model.groundModel is None


def test_configure_blockage_model_ground_mirror_fuga_warns(tmp_path):
    # FUGA LUTs already include the ground; the flag must warn and be ignored.
    # A real LUT is not needed to hit the warning path, but FugaDeficit
    # requires a valid path, so only the warning is asserted before the
    # constructor fails on the dummy path.
    with pytest.warns(UserWarning, match="ground_mirror is ignored for FUGA"):
        try:
            _configure_blockage_model(
                {"name": "FUGA", "ground_mirror": True},
                {"LUT_path": str(tmp_path / "missing.nc")},
            )
        except Exception:
            pass


def test_configure_blockage_model_ground_mirror_increases_deficit():
    # Physics sanity check: the image rotor must increase the upstream
    # blockage deficit relative to the unmirrored model.
    import numpy as np
    from py_wake.deficit_models.no_wake import NoWakeDeficit
    from py_wake.examples.data.hornsrev1 import V80
    from py_wake.site._site import UniformSite
    from py_wake.wind_farm_models.engineering_models import All2AllIterative

    site = UniformSite([1], ti=0.06)
    wt = V80()

    def upstream_ws(cfg):
        wfm = All2AllIterative(
            site,
            wt,
            NoWakeDeficit(),
            blockage_deficitModel=_configure_blockage_model(cfg, {}),
        )
        # second turbine acts as a probe 3D upstream (wd=270 blows along +x)
        sim_res = wfm([0, -240], [0, 0], wd=270, ws=10)
        return sim_res.WS_eff.sel(wt=1).item()

    ws_plain = upstream_ws({"name": "SelfSimilarityDeficit2020"})
    ws_mirror = upstream_ws(
        {"name": "SelfSimilarityDeficit2020", "ground_mirror": True}
    )
    assert ws_mirror < ws_plain < 10.0
    assert np.isclose(ws_plain, ws_mirror, atol=0.1)  # small correction, same order


def test_ground_mirror_works_with_squared_wake_superposition():
    # Regression: Mirror without an explicit superposition falls back to the
    # wind farm model's superposition. With SquaredSum (e.g. the TurbOPark
    # recipe) that asserts on the speed-up (negative deficit) regions every
    # blockage model produces ("SquaredSum only works for deficit - not
    # speedups"). The mirror must therefore sum real+image linearly.
    from py_wake.deficit_models.utils import ct2a_madsen
    from py_wake.examples.data.hornsrev1 import V80
    from py_wake.site._site import UniformSite
    from py_wake.superposition_models import SquaredSum
    from py_wake.turbulence_models.stf import STF2017TurbulenceModel
    from py_wake.wind_farm_models.engineering_models import All2AllIterative
    from py_wake.deficit_models.gaussian import TurboGaussianDeficit

    site = UniformSite([1], ti=0.06)
    wt = V80()
    wfm = All2AllIterative(
        site,
        wt,
        TurboGaussianDeficit(ct2a=ct2a_madsen),
        superpositionModel=SquaredSum(),
        turbulenceModel=STF2017TurbulenceModel(),
        blockage_deficitModel=_configure_blockage_model(
            {"name": "VortexCylinder", "ground_mirror": True}, {}
        ),
    )
    sim_res = wfm([0, 400], [0, 0], wd=270, ws=10)
    assert (sim_res.WS_eff < 10.0).all()


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


def _weighted_system(rotor_name, ws_superposition="Weighted"):
    return {
        "attributes": {
            "analysis": {
                "wind_deficit_model": {
                    "name": "Zong2020",
                    "wake_expansion_coefficient": {"k_a": 0.38, "k_b": 0.004},
                },
                "superposition_model": {"ws_superposition": ws_superposition},
                "rotor_averaging": {"name": rotor_name},
                "deflection_model": {"name": "None"},
                "turbulence_model": {"name": "None"},
                "blockage_model": {"name": None},
            }
        }
    }


@pytest.mark.parametrize("rotor_name", ["gaussian_overlap", "area_overlap", "center"])
@pytest.mark.parametrize("ws_superposition", ["Weighted", "Cumulative"])
def test_weighted_superposition_requires_node_rotor_avg(rotor_name, ws_superposition):
    """WeightedSum/CumulativeWakeSum with a non-node rotor-averaging model raises
    a clear ValueError instead of PyWake's deep AssertionError."""
    with pytest.raises(ValueError, match="node"):
        configure_wake_model(
            _weighted_system(rotor_name, ws_superposition),
            rotor_diameter=126.0,
            hub_height=90.0,
        )


@pytest.mark.parametrize("rotor_name", ["grid", "eq_grid", "gq_grid", "cgi"])
def test_weighted_superposition_allows_node_rotor_avg(rotor_name):
    """A node rotor-averaging model is accepted with Weighted superposition."""
    config = configure_wake_model(
        _weighted_system(rotor_name), rotor_diameter=126.0, hub_height=90.0
    )
    assert isinstance(config["superposition_model"], WeightedSum)


def _weighted_deficit_system(deficit_name):
    """Weighted superposition + node rotor-avg, varying only the deficit."""
    return {
        "attributes": {
            "analysis": {
                "wind_deficit_model": {
                    "name": deficit_name,
                    "wake_expansion_coefficient": {"free_stream_ti": True},
                },
                "superposition_model": {"ws_superposition": "Weighted"},
                "rotor_averaging": {"name": "grid"},
                "deflection_model": {"name": "None"},
                "turbulence_model": {"name": "None"},
                "blockage_model": {"name": None},
            }
        }
    }


@pytest.mark.parametrize("deficit_name", ["SuperGaussian", "SuperGaussian2023", "GCL"])
def test_weighted_superposition_requires_convection_deficit(deficit_name):
    """WeightedSum/CumulativeWakeSum with a non-ConvectionDeficitModel deficit
    (super-Gaussian, GCL) raises a clear ValueError instead of PyWake's deep
    AssertionError, even when the rotor-averaging model is a node model."""
    with pytest.raises(ValueError, match="ConvectionDeficitModel"):
        configure_wake_model(
            _weighted_deficit_system(deficit_name),
            rotor_diameter=126.0,
            hub_height=90.0,
        )


@pytest.mark.parametrize("deficit_name", ["Zong2020", "Niayifar2016", "Bastankhah2014"])
def test_weighted_superposition_allows_convection_deficit(deficit_name):
    """Convection-based deficits are accepted with Weighted superposition."""
    config = configure_wake_model(
        _weighted_deficit_system(deficit_name), rotor_diameter=126.0, hub_height=90.0
    )
    assert isinstance(config["superposition_model"], WeightedSum)


# ---------------------------------------------------------------------------
# CrespoHernandez calibration coefficients (Phase 2 / Fix C)
# ---------------------------------------------------------------------------


def test_crespo_default_without_c():
    """No c -> the PyWake-default CrespoHernandez (ct2a_madsen)."""
    tm = _configure_turbulence_model({"name": "CrespoHernandez"})
    assert isinstance(tm, CrespoHernandez)
    assert tm.ct2a is ct2a_madsen


def test_crespo_with_c_uses_literature_recipe():
    """c -> CrespoHernandez with those coefficients, 1D induction and SqrMaxSum."""
    c = [0.73, 0.83, 0.03, -0.32]
    tm = _configure_turbulence_model({"name": "CrespoHernandez", "c": c})
    assert isinstance(tm, CrespoHernandez)
    assert list(tm.c) == c
    assert tm.ct2a is ct2a_mom1d
    assert isinstance(tm.addedTurbulenceSuperpositionModel, SqrMaxSum)


# ---------------------------------------------------------------------------
# 'none' rotor averaging + WeightedSum (Phase 2 / Fix D)
# ---------------------------------------------------------------------------


def test_rotor_averaging_none():
    assert _configure_rotor_averaging({"name": "none"}) is None


def test_weighted_superposition_allows_none_rotor():
    """WeightedSum accepts rotorAvgModel=None (rotor centre), as Zong (2020) uses."""
    sd = _weighted_deficit_system("Zong2020")
    sd["attributes"]["analysis"]["rotor_averaging"] = {"name": "none"}
    config = configure_wake_model(sd, rotor_diameter=126.0, hub_height=90.0)
    assert config["rotor_averaging"] is None
    assert isinstance(config["superposition_model"], WeightedSum)


def test_weighted_superposition_rejects_center_rotor():
    """A non-node, non-None rotor model is still rejected for WeightedSum."""
    sd = _weighted_deficit_system("Zong2020")
    sd["attributes"]["analysis"]["rotor_averaging"] = {"name": "center"}
    with pytest.raises(ValueError, match="node"):
        configure_wake_model(sd, rotor_diameter=126.0, hub_height=90.0)


# ---------------------------------------------------------------------------
# Zong ceps -> eps_coeff (Phase 2)
# ---------------------------------------------------------------------------


def test_zong_ceps_maps_to_eps_coeff():
    """Zong's near-wake epsilon is named eps_coeff (not ceps) in PyWake."""
    cls, args = _call_deficit(
        "Zong2020",
        {"wake_expansion_coefficient": {"k_a": 0.38, "k_b": 0.004}, "ceps": 0.35},
    )
    assert cls is ZongGaussianDeficit
    assert args["eps_coeff"] == 0.35
    assert "ceps" not in args


# ---------------------------------------------------------------------------
# Axial induction -> ct2a (honoring axial_induction_model)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "axial, expected",
    [("1D", ct2a_mom1d), ("Madsen", ct2a_madsen), ("madsen", ct2a_madsen)],
)
def test_axial_induction_sets_ct2a(axial, expected):
    """axial_induction_model maps to the deficit's ct2a on a ct2a-capable model."""
    _, args = _call_deficit(
        "Bastankhah2014", analysis_top={"axial_induction_model": axial}
    )
    assert args["ct2a"] is expected


def test_axial_induction_default_absent_keeps_deficit_default():
    """No axial_induction_model -> no ct2a override (deficit keeps its default)."""
    _, args = _call_deficit("Bastankhah2014")
    assert "ct2a" not in args


def test_axial_induction_skipped_when_deficit_has_no_ct2a():
    """Blondel2020 has no ct2a parameter, so it is left untouched."""
    _, args = _call_deficit(
        "SuperGaussian", analysis_top={"axial_induction_model": "1D"}
    )
    assert "ct2a" not in args


def test_axial_induction_ct2a_instantiates():
    """The injected ct2a is a valid constructor argument."""
    cls, args = _call_deficit(
        "Niayifar2016", analysis_top={"axial_induction_model": "1D"}
    )
    assert isinstance(cls(**args), NiayifarGaussianDeficit)


# ---------------------------------------------------------------------------
# TurbOPark canonical recipe (Nygaard 2022)
# ---------------------------------------------------------------------------


def test_turbopark_recipe_ground_and_ctlim():
    """TurbOPark gets a Mirror ground model and ctlim=0.96 as constructor args."""
    cls, args, post = _call_deficit_full("TurbOPark")
    assert cls is TurboGaussianDeficit
    assert isinstance(args["groundModel"], Mirror)
    assert args["ctlim"] == 0.96


def test_turbopark_recipe_ws_key_post_attr():
    """TurbOPark scales by the downstream ambient WS via WS_key='WS_jlk'."""
    _, _, post = _call_deficit_full("TurbOPark")
    assert post == {"WS_key": "WS_jlk"}


def test_turbopark_recipe_instantiates_and_applies_ws_key():
    """The recipe builds a valid deficit and WS_key is set post-construction."""
    cls, args, post = _call_deficit_full("TurbOPark")
    deficit = cls(**args)
    for attr, value in post.items():
        setattr(deficit, attr, value)
    assert deficit.WS_key == "WS_jlk"


def test_non_turbopark_has_no_post_attrs():
    """Only TurbOPark carries post-construction attributes."""
    _, _, post = _call_deficit_full("Bastankhah2014")
    assert post == {}


# ---------------------------------------------------------------------------
# use_effective_ws / use_effective_ti for GCL (Phase 3)
# ---------------------------------------------------------------------------


def test_use_effective_ws_honored():
    """The windIO use_effective_ws flag is passed through (not hardcoded)."""
    _, args = _call_deficit("Bastankhah2014", {"use_effective_ws": False})
    assert args["use_effective_ws"] is False


def test_use_effective_ws_defaults_true():
    _, args = _call_deficit("Bastankhah2014")
    assert args["use_effective_ws"] is True


def test_gcl_honors_free_stream_ti():
    """GCLDeficit accepts use_effective_ti (GCLLocal); free_stream_ti is honored."""
    _, args = _call_deficit(
        "GCL", {"wake_expansion_coefficient": {"free_stream_ti": False}}
    )
    assert args["use_effective_ti"] is True


# ---------------------------------------------------------------------------
# Fuga LUT-atmosphere helpers (pure functions; no LUT generation)
# ---------------------------------------------------------------------------
import numpy as np  # noqa: E402


def _resource(ti=None, abl=None, z0=None):
    wr = {}
    if ti is not None:
        wr["turbulence_intensity"] = {"data": np.asarray(ti, dtype=float)}
    if abl is not None:
        wr["ABL_height"] = {"data": np.asarray(abl, dtype=float)}
    if z0 is not None:
        wr["z0"] = {"data": np.asarray(z0, dtype=float)}
    return {"wind_resource": wr}


def test_fuga_atmosphere_derives_z0_from_ti():
    # Neutral inversion: z0 = zhub * exp(-1/TI). Mean TI 0.10 at hub 78.
    z0, zi, zeta0, ti = _fuga_atmosphere(
        _resource(ti=np.full(50, 0.10), abl=np.full(50, 600.0)), {}, 78.0
    )
    assert zeta0 == 0.0
    assert zi == 600.0  # from ABL_height
    assert ti == pytest.approx(0.10)
    assert z0 == pytest.approx(78.0 * np.exp(-1 / 0.10), rel=1e-3)


def test_fuga_atmosphere_defaults_without_resource():
    z0, zi, zeta0, ti = _fuga_atmosphere(None, {}, 78.0)
    assert z0 == 0.03 and zi == 500.0 and ti is None


def test_fuga_atmosphere_config_overrides_site():
    z0, zi, _, _ = _fuga_atmosphere(
        _resource(ti=np.full(10, 0.10)), {"z0": 0.001, "zi": 800.0}, 78.0
    )
    assert z0 == 0.001 and zi == 800.0


def test_fuga_z0_sweep_spans_distribution_and_is_sorted():
    rng = np.random.default_rng(0)
    ti = np.clip(rng.normal(0.10, 0.03, 5000), 0.03, 0.25)
    z0s = _fuga_z0_sweep(_resource(ti=ti), {"n_z0": 5}, 78.0, 0.0, 0.03)
    assert len(z0s) >= 2
    assert z0s == sorted(z0s)
    # all physical: clamp keeps z0 within ~[1e-5, 0.3] m
    assert z0s[0] >= 1e-5 and z0s[-1] <= 0.31


def test_fuga_z0_sweep_clamps_high_ti_to_physical_z0():
    # All-high TI must not produce absurd roughness (TI 0.30 -> z0 ~2.8 m).
    z0s = _fuga_z0_sweep(_resource(ti=np.full(200, 0.30)), {"n_z0": 5}, 78.0, 0.0, 0.5)
    assert max(z0s) <= 0.31


def test_fuga_z0_sweep_single_when_n1_or_no_ti():
    assert _fuga_z0_sweep(
        _resource(ti=np.full(10, 0.1)), {"n_z0": 1}, 78.0, 0.0, 0.04
    ) == [0.04]
    assert _fuga_z0_sweep(None, {"n_z0": 5}, 78.0, 0.0, 0.04) == [0.04]


def test_fuga_z0_sweep_explicit_z0_overrides():
    z0s = _fuga_z0_sweep(
        _resource(ti=np.full(10, 0.1)), {"z0": [0.01, 0.05]}, 78.0, 0.0, 0.04
    )
    assert z0s == [0.01, 0.05]


def test_fuga_atmosphere_ignores_list_z0():
    # An explicit z0 list is a sweep (handled by _fuga_z0_sweep); _fuga_atmosphere
    # must not choke on it when computing the scalar fallback.
    z0, zi, _, _ = _fuga_atmosphere(
        _resource(ti=np.full(10, 0.10)), {"z0": [1e-4, 1e-2]}, 78.0
    )
    assert isinstance(z0, float)  # derived from TI, not the list


# ---------------------------------------------------------------------------
# Fuga mixed turbine geometries (LUT generation intercepted; no pyfuga)
#
# Regression for mixed-rotor-diameter/hub-height layouts (neighbor farms):
# single-hub-height LUTs at different hub heights turn FugaMultiLUTDeficit's
# merged table all-NaN (xarray cannot interpolate a size-1 z axis), which
# surfaced as zero production at every waked turbine. Mixed layouts must get
# LUTs spanning all hub heights, one shared x/y grid, and a z_lst pinning the
# merged z axis to the hub heights.
# ---------------------------------------------------------------------------
import wifa.pywake_api as _pywake_api  # noqa: E402


def _call_fuga(monkeypatch, turbine_geometries=None):
    """Call the FUGA deficit branch with _ensure_fuga_luts intercepted."""
    captured = {}

    def fake_ensure_luts(**kwargs):
        captured.update(kwargs)
        n = len(kwargs["geometries"]) * len(kwargs["z0_list"])
        return [f"/fake/lut_{i}.nc" for i in range(n)]

    monkeypatch.setattr(_pywake_api, "_ensure_fuga_luts", fake_ensure_luts)
    analysis = {
        "wind_deficit_model": {"name": "FUGA", "fuga": {"z0": 0.03, "n_z0": 1}}
    }
    _, args, _ = _configure_deficit_model(
        {"name": "FUGA"}, analysis, _RD, _HH, turbine_geometries=turbine_geometries
    )
    return args, captured


def test_fuga_mixed_geometries_span_hub_heights_and_share_grid(monkeypatch):
    args, luts = _call_fuga(monkeypatch, [(80.0, 80.0), (100.0, 90.0)])
    # Every LUT spans all hub heights so the merged z axis has no NaN slices.
    assert luts["zlow"] == 80.0 and luts["zhigh"] == 90.0
    # One shared x/y grid (finest natural resolution: min diameter / 4, / 16).
    assert luts["dx"] == 20.0 and luts["dy"] == 5.0
    # The merged z axis is pinned to exactly the hub heights.
    assert args["z_lst"] == [80.0, 90.0]
    assert isinstance(args["LUT_path"], list) and len(args["LUT_path"]) == 2


def test_fuga_mixed_diameters_same_hub_share_grid_only(monkeypatch):
    args, luts = _call_fuga(monkeypatch, [(80.0, 80.0), (100.0, 80.0)])
    assert luts["dx"] == 20.0 and luts["dy"] == 5.0
    # Same hub height everywhere: single-level LUTs stay valid, no z axis.
    assert "zlow" not in luts and "z_lst" not in args


def test_fuga_single_geometry_keeps_hub_level_lut(monkeypatch):
    args, luts = _call_fuga(monkeypatch, None)
    # Unchanged cheap path: per-geometry defaults, single plain LUT path.
    assert "zlow" not in luts and "dx" not in luts
    assert "z_lst" not in args
    assert isinstance(args["LUT_path"], str)


def test_fuga_duplicate_geometries_deduped(monkeypatch):
    # Two turbine types sharing one geometry must not produce duplicate d_h
    # coordinates (they break FugaMultiLUTDeficit's merge) nor a second LUT.
    args, luts = _call_fuga(monkeypatch, [(80.0, 80.0), (80.0, 80.0)])
    assert luts["geometries"] == [(80.0, 80.0)]
    assert isinstance(args["LUT_path"], str)
