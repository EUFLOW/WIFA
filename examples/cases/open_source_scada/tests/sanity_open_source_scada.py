"""Sanity check for examples/cases/open_source_scada/outputs/observedPower.nc.

Verifies that the cleaned per-turbine power timeseries is consistent with
the Vestas V82 1.65 MW OEM power curve packaged in the same case
directory. Catches regressions of the historical "multiply by rated MW"
bug fixed in the cleanup PR.

Run from the WIFA repository root:

    python examples/cases/open_source_scada/tests/sanity_open_source_scada.py
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import xarray as xr
import yaml


CASE_DIR = Path(__file__).resolve().parents[1]
OBSERVED = CASE_DIR / "outputs" / "observedPower.nc"
TURBINE = CASE_DIR / "plant_energy_turbine" / "turbine.yaml"

RATED_W = 1.65e6
TOL_RATED = 0.02 * RATED_W            # 2% above rated tolerated (overshoot)
TOL_BIN_FRAC = 0.20                   # 20% bin-median tolerance vs OEM
WS_TEST_BINS = list(range(5, 14))     # exclude very-low (start-up scatter)
                                      # and very-high (storm shutdown scatter)


def load_oem_pc(turbine_yaml: Path):
    d = yaml.safe_load(turbine_yaml.read_text())
    pc = d["performance"]["power_curve"]
    ws = np.array(pc["power_wind_speeds"], dtype=float)
    p_kw = np.array(pc["power_values"], dtype=float)
    return ws, p_kw * 1000.0          # convert to W


def main():
    print(f"Checking {OBSERVED} ...")
    ds = xr.open_dataset(OBSERVED)
    P = ds["power"].values             # (time, turbine), W
    WS = ds["rotor_effective_velocity"].values
    oem_ws, oem_p = load_oem_pc(TURBINE)

    valid = np.isfinite(P) & np.isfinite(WS)
    p_max = float(np.nanmax(P))
    print(f"  max(power) = {p_max:.1f} W (rated = {RATED_W:.1f} W)")
    assert p_max <= RATED_W + TOL_RATED, (
        f"max(power)={p_max} exceeds rated+tol={RATED_W + TOL_RATED}; "
        "is the historical 'multiply by rated MW' scaling re-introduced?"
    )

    fails = []
    for lo in WS_TEST_BINS:
        m = valid & (WS >= lo) & (WS < lo + 1)
        if m.sum() < 30:
            continue
        oem_at = float(np.interp(lo + 0.5, oem_ws, oem_p))
        med = float(np.median(P[m]))
        rel = abs(med - oem_at) / max(oem_at, 1.0)
        flag = "OK" if rel <= TOL_BIN_FRAC else "FAIL"
        print(f"  WS={lo:2d}-{lo+1:2d}: n={m.sum():6d}  "
              f"OEM={oem_at:7.0f} W  median={med:7.0f} W  rel-err={rel:.2f} [{flag}]")
        if rel > TOL_BIN_FRAC:
            fails.append((lo, oem_at, med, rel))

    if fails:
        msg = ", ".join(f"WS={lo}: rel-err={r:.2f}" for lo, _, _, r in fails)
        raise AssertionError(
            f"Power-curve agreement worse than {TOL_BIN_FRAC:.0%} on bins: {msg}"
        )

    nan_frac = float(np.isnan(P).mean())
    print(f"  NaN fraction (operational gaps): {nan_frac:.2f}  (expect ~0.20)")
    print("Sanity check passed.")


if __name__ == "__main__":
    main()
