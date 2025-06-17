# WIFA (Wind Farm API)

[![Tests](https://github.com/EUFLOW/WIFA/workflows/Run%20tests/badge.svg)](https://github.com/EUFLOW/WIFA/actions)
[![codecov](https://codecov.io/gh/WIFA/workflows/branch/main/graph/badge.svg)](https://codecov.io/gh/EUFLOW/WIFA)


WIFA is an open-source wind farm prediction tool that interfaces models with various degrees of fidelity. It includes engineering wake models (PyWake and foxes), an Atmospheric Perturbation Model or APM (wayve), and a 3D CFD code (code_saturne).

The different tools embedded in WIFA are interfaced with a Python wrapping through an API. All the models available in WIFA rely on a common standardized input data structure (WindIO) to return power production estimations.

A common set of example cases can be run by the different tools interfaced in WIFA. The example cases cover the range of expected types of data provided as input (wind speed and wind direction at hub height, vertical profiles, etc.) and simulation/post-processing options.

## Documentation

https://euflow.github.io/WIFA/
