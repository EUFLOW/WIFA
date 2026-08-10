# Multiple wind farms example

windIO wind energy system with three wind farms sharing one site and one
energy resource, used to exercise multi-farm support in the WIFA runners
(`wind_farm` given as a list).

Runner support differs:

- **pywake** (`run_pywake`) returns one AEP per farm (a list, in input order),
  with wake interaction between the farms accounted for.
- **foxes** (`run_foxes`, foxes >= 1.8.4) reads the farm list but merges all
  farms into a single combined farm result; per-farm AEPs are not currently
  reported.
- **wayve**, **floris**, and **code_saturne** raise `NotImplementedError` for
  multi-farm input.

## Attribution

The plant data (layouts, turbine, site, and energy resource) are derived from
the IEA Wind 2.2-GW, 22-MW Reference Offshore Wind Plant:

- Source: [IEAWindSystems/IEA-Wind-2200-22-ROWP](https://github.com/IEAWindSystems/IEA-Wind-2200-22-ROWP)
- License: Apache License 2.0 (see the source repository for the full text)

The data have been adapted for use as a WIFA example. WIFA itself is
MIT-licensed; the Apache-2.0 license of the upstream dataset applies to the
files in this directory that are derived from it.
