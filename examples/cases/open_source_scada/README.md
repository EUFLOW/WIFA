# open_source_scada — Twin Groves Wind Farm

This example case packages one year (2010-09-01 to 2011-08-31, hourly) of
per-turbine SCADA power and wind-speed observations from the **Twin
Groves Wind Farm**, owned and operated by Invenergy in McLean County,
Illinois, USA. The farm has approximately 200 Vestas V82 1.65 MW
turbines (82 m rotor, 80 m hub height).

The packaged data is suitable for benchmarking flow models against
operational observations: it provides a public reference case with a
known turbine type, a publicly available manufacturer power curve, and
permissive licensing.

## Files

| path | content |
|--|--|
| `plant_wind_farm/FLOW_toy_study_wind_farm.yaml` | Turbine x/y coordinates (UTM zone 16N, EPSG:32616). |
| `plant_energy_turbine/turbine.yaml` | Vestas V82 1.65 MW power and thrust curves (OEM values, also at https://github.com/NatLabRockies/turbine-models). |
| `plant_energy_turbine/VestasV82_1.65MW_82.csv` | Same OEM curve in CSV form (WS, P, Cp, Ct). |
| `plant_energy_resource/` | ERA5 reanalysis wind timeseries at hub height + Weibull resource + utility scripts. |
| `plant_energy_site/FLOW_toy_study_energy_site.yaml` | Site polygon and energy-resource include. |
| `wind_energy_system/system.yaml` | Top-level WindIO definition tying farm, site, and analysis together. |
| `wind_energy_system/analysis.yaml` | Example Bastankhah Gaussian wake-model analysis configuration. |
| `outputs/observedPower.nc` | **Observed** per-turbine power (W) and wind speed (m/s) per hour, shape `(time=8760, turbine=200)`. |

## Units and conventions

- `observedPower.nc.power`: per-turbine electrical power output in **watts (W)**.
  Maximum across all (time, turbine) entries is the rated power
  $1.65\times10^{6}$ W, consistent with the V82 OEM curve packaged in
  `plant_energy_turbine/turbine.yaml`.
- `observedPower.nc.rotor_effective_velocity`: per-turbine wind speed in
  **m/s**. Currently the raw nacelle-anemometer reading; this is **not** a
  freestream value (no Nacelle Transfer Function correction has been
  applied, and no upstream met mast is provided).
- Time axis is hourly UTC. Approximately 20% of (time, turbine) entries
  are NaN due to operational gaps, maintenance, and turbine outages — this
  is real operational SCADA, not a synthetic dataset.

## Historical note on the power scaling

Prior to the cleanup PR that accompanies this README, the
`observedPower.nc.power` values were stored with a **multiplicative factor
of 1.65 applied** (i.e., the rated-power numerical value in MW). Maximum
values were near 2.72 MW, which is inconsistent with the V82 1.65 MW
rated power and made the file unusable for any benchmark relying on the
packaged OEM power curve.

If you have a copy of the older file, you can recover the corrected
values by dividing `power` by `1.65`:

```python
import xarray as xr
ds = xr.open_dataset("outputs/observedPower.nc")
ds["power"] = ds["power"] / 1.65   # only if you have the pre-PR file
```

The current file in this repository has the correction applied. A
self-consistency check is included in `tests/sanity_open_source_scada.py`
(see PR description).

## Suggested usage

```python
import xarray as xr
ds = xr.open_dataset("outputs/observedPower.nc")
print(ds)                       # (time: 8760, turbine: 200)
print(ds.power.max().item())    # ~1.65e6 W
print(ds.power.attrs)           # units, source
```

## License

MIT (inherited from the WIFA repository).

## Provenance and citation

If you use this dataset in published work, please cite the WIFA
repository and acknowledge Twin Groves Wind Farm / Invenergy as the
source of the SCADA observations. The Vestas V82 1.65 MW OEM power and
thrust curves are reproduced from public manufacturer data (also
available at the NatLabRockies turbine-models repository:
https://github.com/NatLabRockies/turbine-models/blob/main/turbine_models/data/Onshore/VestasV82_1.65MW_82.csv).
