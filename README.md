# WIFA (Wind Farm API)

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![Documentation](https://img.shields.io/badge/docs-online-green.svg)](https://euflow.github.io/WIFA/)

WIFA is an open-source multi-fidelity wind farm simulation framework that provides a unified interface to multiple flow modeling tools.

## Features

- **Multi-fidelity support**: Engineering wake models (PyWake, foxes), atmospheric perturbation model (wayve), and CFD (code_saturne)
- **Unified interface**: Same input format works with all tools via WindIO schema
- **CLI and Python API**: Run simulations from command line or integrate into workflows
- **Flexible inputs**: Time series, Weibull wind roses, turbine-specific conditions
- **Comprehensive outputs**: Turbine power, AEP, flow fields in NetCDF format

## Supported Tools

| Tool | Type | Speed | Use Case |
|------|------|-------|----------|
| **PyWake** | Engineering wake model | Fast | AEP estimation, layout optimization |
| **foxes** | Engineering wake model | Fast | Large farms, long time series |
| **wayve** | Atmospheric perturbation model | Medium | Gravity waves, farm blockage |
| **code_saturne** | CFD (RANS) | Slow (HPC) | Detailed flow analysis |

## Quick Start

### Installation

We recommend using [uv](https://docs.astral.sh/uv/) for fast, reliable Python environment management:

```bash
# Install uv (if not already installed)
curl -LsSf https://astral.sh/uv/install.sh | sh

# Create environment and install WIFA
uv venv --python 3.11
source .venv/bin/activate
uv pip install wifa
```

Or with pip:

```bash
pip install wifa
```

### Run a Simulation

```bash
# Using CLI
wifa examples/cases/windio_4turbines/wind_energy_system/system.yaml

# Tool-specific
wifa_foxes system.yaml
wifa_pywake system.yaml
```

### Python API

```python
from wifa.main_api import run_api

# Run simulation (tool selected from YAML)
run_api("path/to/system.yaml")

# Or use tool-specific API
from wifa.foxes_api import run_foxes
results = run_foxes("system.yaml", engine="process", n_procs=4)
```

## Input Format

WIFA uses WindIO-formatted YAML files:

```yaml
# system.yaml
name: My Wind Farm Study
site: !include ../site/site.yaml
wind_farm: !include ../farm/farm.yaml
attributes:
  flow_model:
    name: foxes  # or pywake, wayve, codesaturne
  analysis: !include analysis.yaml
  model_outputs_specification:
    output_folder: "results"
    turbine_outputs:
      turbine_nc_filename: 'turbine_data.nc'
```

See the [documentation](https://euflow.github.io/WIFA/) for complete schema details.

## Examples

The `examples/cases/` directory contains ready-to-run cases:

- `windio_4turbines` - Basic 4-turbine example
- `simple_wind_rose` - Weibull wind rose input
- `heterogeneous_wind_rose_at_turbines` - Per-turbine wind conditions
- `timeseries_with_operating_flag` - Turbine availability modeling
- `windio_4turbines_ABL` - Atmospheric boundary layer case
- And more...

## Architecture

![WIFA Architecture](docs/img/wifa_diagram.png)

## Documentation

Full documentation: https://euflow.github.io/WIFA/

- [Getting Started](https://euflow.github.io/WIFA/getting_started.html)
- [Installation Guide](https://euflow.github.io/WIFA/installation.html)
- [Examples](https://euflow.github.io/WIFA/examples.html)
- [API Reference](https://euflow.github.io/WIFA/main_api.html)
- [Troubleshooting](https://euflow.github.io/WIFA/troubleshooting.html)

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch
3. Set up development environment:
   ```bash
   uv venv --python 3.11
   source .venv/bin/activate
   uv pip install -e ".[dev,test]"
   ```
4. Make your changes
5. Run tests: `pytest tests/`
6. Submit a pull request

## Citation

If you use WIFA in your research, please cite:

```bibtex
@software{wifa2024,
  title = {WIFA: Wind Farm API},
  author = {EU-FLOW Consortium},
  year = {2024},
  url = {https://github.com/EUFLOW/WIFA}
}
```

## License

MIT License - see [LICENSE](LICENSE) for details.

## Acknowledgments

WIFA is developed within the EU-FLOW project. The embedded tools are developed by:

- **PyWake**: DTU Wind Energy
- **foxes**: Fraunhofer IWES
- **wayve**: KU Leuven
- **code_saturne**: EDF
