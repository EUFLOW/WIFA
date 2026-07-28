# General packages
import argparse
import warnings
from pathlib import Path

import numpy as np
import xarray as xr
from scipy.interpolate import interp1d
from scipy.special import gamma as scipy_gamma
from windIO import load_yaml

from wifa._optional import require


def run_wayve(yamlFile, output_dir="output", debug_mode=False):
    require("wayve")

    # General APM setup
    from wayve.apm import APM
    from wayve.grid.grid import Stat2Dgrid
    from wayve.momentum_flux_parametrizations import FrictionCoefficients
    from wayve.solvers import FixedPointIteration

    #####################
    # Read out yaml file
    #####################
    # Yaml loading
    if isinstance(yamlFile, dict):
        system_dat = yamlFile
    else:
        system_dat = load_yaml(yamlFile)
    # WindIO components
    farm_dat = system_dat["wind_farm"]
    resource_dat = system_dat["site"]["energy_resource"]
    analysis_dat = system_dat["attributes"]["analysis"]

    ######################
    # construct APM grid
    ######################
    # Default numerical parameters (values from Allaerts and Meyers, 2019)
    Lx = 1.0e6  # grid size in x-direction [m]
    Ly = 1.0e6  # grid size in y-direction [m]
    dx = 500.0  # Grid spacing
    L_filter = 1.0e3
    # Read out numerical parameters
    if "apm_grid" in analysis_dat:
        grid_dat = analysis_dat["apm_grid"]
        if "Lx" in grid_dat:
            Lx = float(grid_dat["Lx"])
        if "Ly" in grid_dat:
            Ly = float(grid_dat["Ly"])
        if "dx" in grid_dat:
            dx = float(grid_dat["dx"])
        if "L_filter" in grid_dat:
            L_filter = float(grid_dat["L_filter"])
    # Grid points
    Nx = int(Lx / dx)  # grid points in x-direction
    Ny = int(Ly / dx)  # grid points in y-direction
    # Generate 2D grid object
    grid = Stat2Dgrid(Lx, Nx, Ly, Ny)

    ####################
    # Set up WindFarm and Forcing objects
    ####################
    wind_farm, forcing, wf_offset_x, wf_offset_y = wf_setup(
        farm_dat, analysis_dat, L_filter, debug_mode
    )
    coupling = wind_farm.coupling
    wake_model = coupling.wake_model
    Nt = wind_farm.Nturb
    hh = np.mean([turb.zh for turb in wind_farm.turbines])
    h1_min = np.max([turb.zh + turb.D / 2 for turb in wind_farm.turbines])

    # Determine H1
    h1 = 2.0 * hh  # Default
    if "layers_description" in analysis_dat:
        if "farm_layer_height" in analysis_dat["layers_description"]:
            h1 = analysis_dat["layers_description"]["farm_layer_height"]
    if h1 < h1_min:
        raise UserWarning("Lower layer height too low, please specify a higher value")

    # Optional ABL-setup knobs (analysis.abl_setup), forwarded to flow_io_abl:
    # geostrophic-wind mode and capping-inversion fit settings.
    abl_setup = analysis_dat.get("abl_setup", {})
    abl_kwargs = {}
    if "Gmode" in abl_setup:
        abl_kwargs["gmode"] = abl_setup["Gmode"]
    if "dh_max" in abl_setup:
        abl_kwargs["dh_max"] = abl_setup["dh_max"]
    if "serz" in abl_setup:
        abl_kwargs["serz"] = bool(abl_setup["serz"])

    ##################
    # Other APM components
    ##################
    # Momentum flux parametrization
    mfp = FrictionCoefficients()
    # Pressure feedback parametrization: any value but 1 selects the
    # profile-resolving NonUniform gravity-wave closure (<= 0 means "one
    # sublayer per profile level"). The object is built per state
    # (_pressure_for_state), because whether the profile actually supports
    # NonUniform depends on the per-state inversion fit.
    n_fa_layers = 1
    if "layers_description" in analysis_dat:
        n_fa_layers = analysis_dat["layers_description"].get("number_of_fa_layers", 1)

    ######################
    # Read output settings
    ######################
    # Select timestamps. `time_indices` are positions into the *full* wind
    # resource arrays; they must be carried alongside the timestamp labels,
    # because flow_io_abl() indexes the full arrays. Enumerating the subsetted
    # timestamps instead would simulate rows 0..n-1 while labelling them with
    # the requested timestamps.
    all_times = resource_dat["wind_resource"]["time"]
    time_indices = list(range(len(all_times)))
    run_config = system_dat["attributes"]["model_outputs_specification"][
        "run_configuration"
    ]
    if "times_run" in run_config and not run_config["times_run"].get(
        "all_occurences", True
    ):
        if "subset" in run_config["times_run"]:
            time_indices = list(run_config["times_run"]["subset"])
    times = [all_times[i] for i in time_indices]
    # Get turbine variables to output
    turbine_nc_filename = "turbine_data.nc"
    turbine_output_variables = ["power", "rotor_effective_velocity"]
    if "turbine_outputs" in system_dat["attributes"]["model_outputs_specification"]:
        turb_out_dat = system_dat["attributes"]["model_outputs_specification"][
            "turbine_outputs"
        ]
        if "turbine_nc_filename" in turb_out_dat:
            turbine_nc_filename = turb_out_dat["turbine_nc_filename"]
        # The schema (and every other WIFA engine) calls this `output_variables`;
        # `turbine_output_variables` is kept as a fallback for legacy yamls.
        if "output_variables" in turb_out_dat:
            turbine_output_variables = turb_out_dat["output_variables"]
        elif "turbine_output_variables" in turb_out_dat:
            turbine_output_variables = turb_out_dat["turbine_output_variables"]
    # Check flow field output specification
    flow_nc_filename = "flow_field.nc"
    flow_output_variables = ["wind_speed", "wind_direction"]
    report_flow = False
    x_ff = []
    y_ff = []
    z_ff = []
    if "flow_field" in system_dat["attributes"]["model_outputs_specification"]:
        if (
            "report"
            in system_dat["attributes"]["model_outputs_specification"]["flow_field"]
        ):
            report_flow = system_dat["attributes"]["model_outputs_specification"][
                "flow_field"
            ]["report"]
    if report_flow:
        # Output settings
        flow_out_dat = system_dat["attributes"]["model_outputs_specification"][
            "flow_field"
        ]
        if "flow_nc_filename" in flow_out_dat:
            flow_nc_filename = flow_out_dat["flow_nc_filename"]
        if "output_variables" in flow_out_dat:
            flow_output_variables = flow_out_dat["output_variables"]
        # Output grid
        if flow_out_dat["z_planes"]["xy_sampling"] != "grid":
            report_flow = False
            raise UserWarning("xy_sampling not supported")
        x_bounds = flow_out_dat["z_planes"]["x_bounds"]
        y_bounds = flow_out_dat["z_planes"]["y_bounds"]
        if "Nx" in flow_out_dat["z_planes"]:
            Nx = flow_out_dat["z_planes"]["Nx"]
        else:
            dx = flow_out_dat["z_planes"]["dx"]
            Nx = int((x_bounds[1] - x_bounds[0]) / dx)
        if "Ny" in flow_out_dat["z_planes"]:
            Ny = flow_out_dat["z_planes"]["Ny"]
        else:
            dy = flow_out_dat["z_planes"]["dy"]
            Ny = int((y_bounds[1] - y_bounds[0]) / dy)
        x_ff = np.linspace(x_bounds[0], x_bounds[1], Nx)
        y_ff = np.linspace(y_bounds[0], y_bounds[1], Ny)
        if flow_out_dat["z_planes"]["z_sampling"] == "hub_heights":
            z_ff = np.unique([turb.zh for turb in wind_farm.turbines])
        else:
            z_ff = flow_out_dat["z_planes"]["z_list"]

    #####################
    # Perform model runs
    #####################
    # Validate the capping-inversion spelling once, before the loop: an
    # ambiguous or incomplete block is a defect of the file, and the loop's
    # ``except Exception`` would otherwise report it as every state crashing.
    capping_inversion_spelling(resource_dat["wind_resource"])
    # Initialize crash counter
    crashes = 0
    # NonUniform free-atmosphere tally (per-state fallback is warned once)
    nonuniform_states = 0
    # List of datasets
    ds_list = []
    ds_ff_list = []
    # Loop over timeseries
    for run_index, time_index in enumerate(time_indices):
        time = times[run_index]
        if debug_mode:
            # Print timestep
            print(f"time {run_index + 1}/{len(times)}")
        try:
            # Set up ABL (solver frame: hub-height wind along +x)
            abl, rotation = flow_io_abl(
                resource_dat["wind_resource"], time_index, hh, h1, **abl_kwargs
            )
            # Rebuild the wind farm in the solver frame: wayve's gravity-wave
            # solve assumes westerly flow, so the layout turns with the wind.
            wind_farm, forcing, wf_offset_x, wf_offset_y = wf_setup(
                farm_dat, analysis_dat, L_filter, debug_mode, rotation=rotation
            )
            coupling = wind_farm.coupling
            wake_model = coupling.wake_model
            # Pressure feedback for this state
            pressure = _pressure_for_state(abl, n_fa_layers)
            if type(pressure).__name__ == "NonUniform":
                nonuniform_states += 1
            # Set up APM from components
            model = APM(grid, forcing, abl, mfp, pressure)
            # Use a fixed-point iteration solver with a relaxation factor of 0.7
            solver = FixedPointIteration(tol=5.0e-3, relax=0.7)
            # Solve model
            if not debug_mode:
                _ = model.solve(method=solver)  # APM run
            else:
                wind_farm.preprocess(model)  # Wake model run
            # Turbine level outputs #
            # Turbine output dictionary
            turb_out_dict = {}
            if "power" in turbine_output_variables:
                turb_out_dict["power"] = (
                    ["turbine"],
                    wind_farm.power_turbines(abl.rho),
                )
            if "rotor_effective_velocity" in turbine_output_variables:
                turb_out_dict["rotor_effective_velocity"] = (
                    ["turbine"],
                    wind_farm.coupling.St,
                )
            # State-level ABL diagnostics. The capping inversion and the
            # geostrophic wind are often *fitted* here rather than given (see
            # flow_io_abl), so without these the atmosphere a run actually used
            # is unrecoverable from its outputs. Scalars, no turbine dim.
            diagnostics = _abl_diagnostics(abl, pressure)
            assert not set(diagnostics) & set(turb_out_dict)
            turb_out_dict |= diagnostics
            # NC setup
            ds = xr.Dataset(
                turb_out_dict,
                coords={"states": time, "turbine": range(Nt)},
            )
            # Flow field outputs #
            ds_ff = None
            if report_flow and not debug_mode:
                # Callables for flow evaluation
                u_bg_evaluator = coupling.set_up_u_bg_evaluator(
                    abl
                )  # Background velocity callable
                apm_evaluator = coupling.apm_evaluator  # APM lower layer state callable
                # Output arrays
                wind_speed = np.zeros([len(x_ff), len(y_ff), len(z_ff)])
                wind_dir = np.zeros([len(x_ff), len(y_ff), len(z_ff)])
                # Query points: the requested earth-frame grid, rotated into
                # the solver frame (where the farm sits and the wind blows
                # along +x). Evaluating the rotated points directly keeps the
                # output on the requested grid with no regridding.
                c, s = np.cos(rotation), np.sin(rotation)
                x_e, y_e = np.meshgrid(
                    x_ff - wf_offset_x, y_ff - wf_offset_y, indexing="ij"
                )
                x_q = c * x_e + s * y_e
                y_q = c * y_e - s * x_e
                # Loop over z-planes
                for k, z_k in enumerate(z_ff):
                    # Get velocities (solver frame)
                    u_bg, v_bg, u_wm, v_wm = _xy_plane_points(
                        wake_model,
                        wind_farm,
                        abl,
                        u_bg_evaluator,
                        apm_evaluator,
                        x_q,
                        y_q,
                        z_k,
                    )
                    # Rotate velocity vectors back to the earth frame
                    u_wm, v_wm = c * u_wm - s * v_wm, s * u_wm + c * v_wm
                    # Convert to speed and direction (wrapped to [0, 360))
                    wind_speed[:, :, k] = np.sqrt(np.square(u_wm) + np.square(v_wm))
                    wind_dir[:, :, k] = (
                        np.rad2deg(np.pi / 2 - (np.arctan2(v_wm, u_wm) + np.pi)) % 360.0
                    )
                # Flow output dictionary
                flow_out_dict = {}
                if "wind_speed" in flow_output_variables:
                    flow_out_dict["wind_speed"] = (["x", "y", "z"], wind_speed)
                if "wind_direction" in flow_output_variables:
                    flow_out_dict["wind_direction"] = (["x", "y", "z"], wind_dir)
                # NC setup
                ds_ff = xr.Dataset(
                    flow_out_dict,
                    coords={"states": time, "x": x_ff, "y": y_ff, "z": z_ff},
                )
            # Append outputs together, so a flow-field failure cannot leave
            # turbine_data.nc and flow_field.nc with different states axes.
            ds_list.append(ds)
            if ds_ff is not None:
                ds_ff_list.append(ds_ff)

        except Exception as exc:
            print(exc)
            # Update crash counter
            crashes += 1
            continue
    if debug_mode:
        print(f"crashes: {crashes}/{len(times)}")
    if n_fa_layers != 1:
        # Make the per-state Uniform/NonUniform closure mixture visible
        print(
            f"NonUniform free atmosphere: {nonuniform_states}/{len(times)} "
            "states (the rest fell back to Uniform)"
        )

    # Combine into total dataset
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    ds_full = xr.concat(ds_list, dim="states")
    output_fn = Path(output_dir) / turbine_nc_filename
    ds_full.to_netcdf(output_fn)
    if report_flow and ds_ff_list:
        ds_ff_full = xr.concat(ds_ff_list, dim="states")
        output_fn = Path(output_dir) / flow_nc_filename
        ds_ff_full.to_netcdf(output_fn)


def _abl_diagnostics(abl, pressure):
    """Per-state record of the atmosphere the solver actually ran on.

    ``flow_io_abl`` fits the capping inversion when the wind resource does not
    state one, derives the geostrophic wind through ``Gmode``, and closes the
    stress profile from surface scalars; ``_pressure_for_state`` then picks the
    free-atmosphere closure per state. All of that is invisible in the turbine
    outputs, so record it alongside them: it is what makes a run reproducible
    and its inversion checkable against an independent boundary-layer height.

    Names are prefixed ``abl_`` so they cannot collide with a turbine output:
    ``air_density`` and ``turbulence_intensity`` are per-turbine fields in other
    engines, and these share one dataset with them.

    Returned as ``{name: ([], value)}`` — scalars on the state axis, no turbine
    dimension.
    """
    gravity = 9.80665  # [m s-2], as in flow_io_abl
    # gprime = g * dtheta / theta0, with theta0 the mixed-layer potential
    # temperature at the inversion; invert it to report dtheta itself.
    dtheta = abl.gprime * np.interp(abl.H, abl.zs, abl.ths) / gravity
    values = {
        "abl_height": (abl.H, "m"),
        "abl_capping_inversion_thickness": (abl.inv_top - abl.inv_bottom, "m"),
        "abl_capping_inversion_strength": (dtheta, "K"),
        "abl_free_atmosphere_N": (abl.N, "s-1"),
        "abl_geostrophic_wind_speed": (abl.S3, "m s-1"),
        # Angle of the geostrophic wind in the solver frame, where the
        # hub-height wind lies along +x: the cross-isobar angle, not a
        # compass direction.
        "abl_geostrophic_veer": (abl.WD3, "degree"),
        # In the turbulence-profile branch wayve's ``utau`` carries the surface
        # stress itself (m2 s-2), not its square root; everywhere else it is a
        # velocity. Reported as wayve holds it.
        "abl_friction_velocity": (abl.utau, "m s-1"),
        "abl_air_density": (abl.rho, "kg m-3"),
        "abl_turbulence_intensity": (abl.TI, "1"),
        # 1 when this state resolved the free atmosphere on the profile,
        # 0 when it fell back to the bulk closure.
        "abl_nonuniform_free_atmosphere": (
            float(type(pressure).__name__ == "NonUniform"),
            "1",
        ),
    }
    # Reporting must not be able to cost a state its results: anything
    # unexpected here (a None, a non-scalar) becomes NaN rather than an
    # exception that the caller's `except Exception` would turn into a drop.
    out = {}
    for name, (value, units) in values.items():
        try:
            out[name] = ([], float(value))
        except (TypeError, ValueError):
            out[name] = ([], float("nan"))
        out[name] = (*out[name], {"units": units})
    return out


def _pressure_for_state(abl, n_fa_layers):
    """Build the gravity-wave pressure closure for one state.

    Anything other than exactly 1 requests wayve's NonUniform closure, which
    resolves the free-atmosphere N(z) and wind shear by slicing the profile
    between the inversion top and ``abl.h_strat`` into layers. A value <= 0 is
    wayve's own "use the ABL's vertical grid" setting: the sublayers follow the
    profile levels themselves, so a reanalysis column resolves the free
    atmosphere at the resolution it actually has instead of an arbitrary count.

    Either way that needs actual profile points in that range: truncated or
    synthetic profiles (e.g. the scalar branch's Nieuwstadt profile, which
    stops at the inversion) fall back to the bulk Uniform closure with a
    warning. The warning text is deliberately state-independent so Python's
    warning dedup collapses it on long time series; run_wayve prints a per-run
    NonUniform/Uniform tally.
    """
    from wayve.pressure.gravity_waves.gravity_waves import NonUniform, Uniform

    if n_fa_layers != 1:
        h_min = max(abl.H, abl.inv_top if abl.inv_top is not None else 0.0)
        n_pts = int(np.sum((abl.zs > h_min) & (abl.zs < abl.h_strat)))
        if n_pts >= 2:
            return NonUniform(n_layers=n_fa_layers, order=1)
        warnings.warn(
            "number_of_fa_layers requested but too few profile points lie "
            "between the inversion top and the free-atmosphere top; falling "
            "back to the Uniform free atmosphere for such states"
        )
    return Uniform(dynamic=True, rotating=False)


def _xy_plane_points(
    wake_model, wind_farm, abl, u_bg_evaluator, apm_evaluator, Xs, Ys, z
):
    """Evaluate the coupled flow at arbitrary (x, y) coordinate arrays.

    Mirror of wayve 2.0.0's ``xy_plane`` methods (Lanzilao and foxes wake
    models) with the axis-aligned ``meshgrid(xs, ys)`` replaced by
    caller-provided 2-D coordinate arrays, so a rotated (solver-frame) grid
    can be evaluated directly — every operation downstream of the meshgrid is
    pointwise in the coordinates. Returns ``(u_bg, v_bg, u_wm, v_wm)`` with
    the shape of ``Xs``. Equivalence with ``xy_plane`` on axis-aligned grids
    is pinned by regression tests for both wake-model paths.
    """
    Nx, Ny = Xs.shape
    if hasattr(wake_model, "_algo"):  # foxes coupling
        import foxes.variables as FV
        from foxes import Engine
        from foxes.utils import wd2uv

        locations = np.stack(
            [np.ravel(Xs), np.ravel(Ys), np.full(Xs.size, float(z))], axis=1
        )
        with Engine.new(**wake_model._engine_pars):
            point_results = wake_model._algo.calc_points(
                wake_model._farm_results,
                locations[None],
                outputs=[FV.AMB_WS, FV.WS, FV.WD],
            )
        amb_uv = wd2uv(
            point_results[FV.WD].to_numpy()[0], point_results[FV.AMB_WS].to_numpy()[0]
        )
        uv = wd2uv(
            point_results[FV.WD].to_numpy()[0], point_results[FV.WS].to_numpy()[0]
        )
        amb_uv = amb_uv.reshape(Nx, Ny, 2)
        uv = uv.reshape(Nx, Ny, 2)
        return amb_uv[..., 0], amb_uv[..., 1], uv[..., 0], uv[..., 1]

    # Lanzilao (UniDirectionalSelfSimilar) path
    from wayve.forcing.wind_farms.wake_model_coupling.wake_models.lanzilao_merging import (
        array_of_matrices,
        dot_matrix_vec_arrays,
    )
    from wayve.forcing.wind_farms.wake_model_coupling.wake_models.wake_model_tools import (
        evaluate_TI,
    )

    Nz = 1
    zs = np.array([z])
    # Get the turbine thrust coefficients
    _, Ct, _ = wake_model.get_St_Ct_et(wind_farm, abl, u_bg_evaluator, apm_evaluator)
    # Get wind farm information
    turbines = wind_farm.turbines
    Nt = wind_farm.Nturb
    xloc = np.array([turbines[k].x for k in range(Nt)])
    yloc = np.array([turbines[k].y for k in range(Nt)])
    # Ambient TI
    TI_inf = abl.TI
    # Get turbine direction (streamwise)
    e_str, e_span = wake_model.background_flow_direction(wind_farm, abl)
    theta_str = np.arctan2(e_str[1], e_str[0])
    # Sort turbines along wind direction
    order = wake_model.sort_turbines(wind_farm, e_str)
    # Turbine direction evaluation
    if wake_model.wake_deflection:  # Base turbine direction on APM velocity
        u_1, v_1, h_1 = apm_evaluator(xloc, yloc)
        theta_turb = np.array([np.arctan2(v_1[i], u_1[i]) for i in range(Nt)])
    else:
        theta_turb = np.array([theta_str for _ in range(Nt)])
    # Vector normal (t) and parallel (p) to rotor
    ets = np.array(
        [np.array([np.cos(theta_turb[i]), np.sin(theta_turb[i])]) for i in order]
    )
    eps = np.array(
        [np.array([-np.sin(theta_turb[i]), np.cos(theta_turb[i])]) for i in order]
    )
    # Sort turbines along wind direction - for TI
    xloc_sort = np.array([xloc[i] for i in order])
    yloc_sort = np.array([yloc[i] for i in order])
    D_sort = np.array([turbines[i].D for i in order])
    Ct_sort = Ct[order]
    theta_turb_sort = np.array([theta_turb[i] for i in order])
    et_sort = np.array([ets[i, :] for i in order])
    ep_sort = np.array([eps[i, :] for i in order])
    # Evaluate TI at turbine locations
    TI = evaluate_TI(
        Nt,
        et_sort,
        ep_sort,
        order,
        xloc_sort,
        yloc_sort,
        D_sort,
        Ct_sort,
        TI_inf,
        theta_turb_sort,
        wake_model.ka,
        wake_model.kb,
    )
    # Set up output
    if wake_model.wake_deflection:
        W_comb = np.zeros((Nx, Ny, Nz, 2)) + e_str[None, None, None, :]
    else:
        W_comb = np.ones((Nx, Ny, Nz))
    for turb in range(Nt):
        if Ct[turb] != 0.0:
            # Evaluate W_t
            W_t = type(wake_model).wake_function(
                Nx,
                Ny,
                Nz,
                Xs,
                Ys,
                zs,
                TI[turb],
                Ct[turb],
                xloc[turb],
                yloc[turb],
                turbines[turb].D,
                turbines[turb].zh,
                theta_turb[turb],
                ka=wake_model.ka,
                kb=wake_model.kb,
                ind=wake_model.induction,
                mirr=wake_model.mirrored,
            )
            if wake_model.wake_deflection:
                # Matrix of current turbine
                et = ets[turb]
                ep = eps[turb]
                A_t = array_of_matrices(1 - W_t, np.outer(et, et)) + np.outer(ep, ep)
                # Update wake deficit field
                W_comb = dot_matrix_vec_arrays(A_t, W_comb)
            else:
                W_comb *= 1 - W_t
    # Wake deficit in main flow direction
    if wake_model.wake_deflection:
        wake_deficit = np.inner(e_str, W_comb)
    else:
        wake_deficit = W_comb
    # Evaluate background velocities
    locations = np.stack(
        [np.ravel(Xs), np.ravel(Ys), np.full(Xs.size, float(z))], axis=1
    )
    vel_bg = u_bg_evaluator(locations)
    u_bg = np.reshape(vel_bg[:, 0], (Nx, Ny, Nz))
    v_bg = np.reshape(vel_bg[:, 1], (Nx, Ny, Nz))
    # Split up into streamwise and spanwise components
    str_bg = u_bg * e_str[0] + v_bg * e_str[1]
    span_bg = u_bg * e_span[0] + v_bg * e_span[1]
    # Add wakes
    str_wm = np.multiply(str_bg, wake_deficit)
    # Convert to u and v components
    u_wm = str_wm * e_str[0] + span_bg * e_span[0]
    v_wm = str_wm * e_str[1] + span_bg * e_span[1]
    return u_bg[:, :, 0], v_bg[:, :, 0], u_wm[:, :, 0], v_wm[:, :, 0]


def nieuwstadt83_profiles(zh, v, wd, z0=1.0e-1, h=1.5e3, fc=1.0e-4, ust=0.666):
    """Set up the cubic analytical profile from Nieuwstadt (1983), based on hub height velocity information"""
    import mpmath

    # Atmospheric state setup
    from wayve.abl.abl_tools import Cg_cubic, alpha_cubic

    # Constants #
    kappa = 0.41  # Von Karman constant
    # # We iterate until we find a profile that has the requested speed at hub height, by varying ust # #
    # Iteration settings
    ust_i = ust
    error = np.inf
    attempt = 0
    max_attempts = 30
    tolerance = 1.0e-3
    # Iteration
    while error > tolerance and attempt < max_attempts:
        # # # Nieuwstadt solution # # #
        # Vertical grid
        Nz = 100
        zs = np.linspace(z0, h, Nz)
        # Dimensionless groups
        hstar = h * fc / ust_i
        z0_h = z0 / h
        # Nieuwstadt relations
        Cg = Cg_cubic(hstar, z0_h, kappa)  # Geostrophic drag Cg = utau/G
        geo_angle = alpha_cubic(hstar, z0_h, kappa)  # Geostrophic wind angle
        # Nieuwstadt solution #
        C = h * fc / kappa / ust_i
        alpha = 0.5 + 0.5 * np.sqrt(1 + 4j * C)
        sigma_s = np.zeros(len(zs), dtype=np.complex128)
        wd_s = np.zeros(len(zs), dtype=np.complex128)
        # Pre-compute gamma values once per iteration (not per z-point)
        gamma_alpha = complex(scipy_gamma(alpha))
        gamma_2alpha = complex(scipy_gamma(2 * alpha))
        prefactor_sigma = alpha * gamma_alpha**2 / gamma_2alpha
        prefactor_wd = (1j * alpha**2 * gamma_alpha**2) / (kappa * C * gamma_2alpha)
        with np.errstate(
            invalid="ignore"
        ):  # z>=h will result in Nan. This is set to 0 below.
            for k in range(len(zs)):
                z_ratio = 1.0 - zs[k] / h
                sigma_s[k] = (
                    prefactor_sigma
                    * np.power(z_ratio, alpha)
                    * complex(mpmath.hyp2f1(alpha - 1, alpha, 2 * alpha, z_ratio))
                )
                wd_s[k] = (
                    prefactor_wd
                    * np.power(z_ratio, alpha - 1)
                    * complex(mpmath.hyp2f1(alpha + 1, alpha - 1, 2 * alpha, z_ratio))
                )
        # Set Nan to 0
        sigma_s[np.isnan(sigma_s)] = np.complex128(0.0)
        wd_s[np.isnan(wd_s)] = np.complex128(0.0)
        # Velocity arrays
        us = ((Cg**-1) * np.cos(geo_angle) + np.real(wd_s)) * ust_i
        vs = ((Cg**-1) * np.sin(geo_angle) + np.imag(wd_s)) * ust_i
        # Error
        u_hh = np.interp(zh, zs, np.sqrt(np.square(us) + np.square(vs)))
        error = np.abs(u_hh - v) / v
        ust_i *= v / u_hh
        attempt += 1
    # Velocity arrays
    us = ((Cg**-1) * np.cos(geo_angle) + np.real(wd_s)) * ust_i
    vs = ((Cg**-1) * np.sin(geo_angle) + np.imag(wd_s)) * ust_i
    # Momentum flux arrays
    tauxs = np.real(sigma_s) * ust_i**2
    tauys = np.imag(sigma_s) * ust_i**2
    nus = (
        kappa
        * ust
        * np.multiply(zs, (1 - zs / h) ** 2, out=np.zeros_like(zs), where=(zs <= h))
    )
    # # Rotate to match wind direction at hub height # #
    # Current wind direction (angle w.r.t. x-axis)
    wd_hh_0 = np.arctan2(np.interp(zh, zs, vs), np.interp(zh, zs, us))
    # Rotation angle
    rotation_angle = -(
        wd_hh_0 + np.deg2rad(wd) + np.pi / 2.0
    )  # +pi/2 for wd convention
    # Velocity components
    us, vs = rotate_xy_arrays(us, vs, rotation_angle)
    tauxs, tauys = rotate_xy_arrays(tauxs, tauys, rotation_angle)
    # Upper atmosphere
    U3 = us[-1]
    V3 = vs[-1]
    return zs, us, vs, U3, V3, tauxs, tauys, nus


def rotate_xy_arrays(xs, ys, angle):
    """
    Rotate the given vectors around the given angle.

    Parameters
    ----------
    xs : array_like
        x-components of the vectors
    ys : array_like
        y-components of the vectors
    angle : float
        angle over which to rotate the vectors (in radians)
    """
    # Angle cosine and sine
    c, s = np.cos(angle), np.sin(angle)
    # Rotation matrix
    R = np.array(((c, -s), (s, c)))
    # Output arrays
    xs_rot, ys_rot = 0.0 * xs, 0.0 * ys
    # Loop over vectors
    Ns = len(xs)
    for i in range(Ns):
        # Vector i
        vec = np.array([xs[i], ys[i]])
        # Multiply with rotation matrix
        rotated_vec = np.matmul(R, vec)
        # Store rotated vector
        xs_rot[i] = rotated_vec[0]
        ys_rot[i] = rotated_vec[1]
    return xs_rot, ys_rot


def ci_fitting(
    zs, ths, l_mo=5.0e3, blh=1.0e3, dh_max=300.0, serz=True, plot_fits=False
):
    import matplotlib.pyplot as plt

    # Atmospheric state setup
    from wayve.abl import ci_methods

    # Stable or unstable atmosphere
    stable = 0.0 < l_mo < 100
    # Estimate inversion parameters with RZ fit #
    # Relevant part of the vertical profiles
    max_z_fit = 5.0e3
    z_ci = zs[zs <= max_z_fit]
    th_ci = ths[zs <= max_z_fit]
    # Surface-Extended RZ or regular RZ
    if serz:
        # Stable or unstable profile determines the initial guess for the CI height
        if stable:
            l_p0 = 1.0e3
        else:
            l_p0 = blh
        # Initial estimate for MBL temperature in fit
        th_mbl = np.interp(l_p0, z_ci, th_ci)
        # Fitting procedure
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            ci_estimate = ci_methods.SERZ_fit(
                z_ci,
                th_ci,
                p0=[0.9, 0.1, 0.0, th_mbl, l_p0, 100, 0.05],
                initialGuess="RZ",
                dh_max=dh_max,
            )
    else:
        # Stable or unstable profile
        if stable:
            # Ignore temperature decrease inside SBL
            # (we are trying to identify the mixing layer that preceded this SBL)
            # We want to capture the mixed layer (or residual layer since we are in
            # SBL) that preceded the SBL. Therefore, we take the potential temperature
            # at the top of the ABL where we have the mixed layer and we extrapolate
            # till the bottom. We use this constant value in the ABL.
            # p0 are the initial guess for [a,b,thm,l,dh] used in Ramp&Zar model
            th_ci[z_ci < blh] = np.interp(blh, z_ci, th_ci)
            l_p0 = 1.0e3
        else:
            # Ignore temperature increase in CBL surface layer, therefore we take
            # the lowest value of potential temperature. We are able to capture the
            # mixed layer in this way, where the temperature is constant and equal
            # to the lowest theta. We use this constant value in the ABL.
            # p0 are the initial guess for [a,b,thm,l,dh] used in Ramp&Zar model
            th_ci[0 : np.argmin(th_ci)] = np.min(th_ci)
            l_p0 = blh
        # Initial estimate for MBL temperature in fit
        th_mbl = np.interp(l_p0, z_ci, th_ci)
        # RZ fit
        ci_estimate = ci_methods.RZfit(
            z_ci, th_ci, p0=[0.9, 0.1, th_mbl, l_p0, 100.0], dh_max=dh_max
        )
    # Plot fitted potential temperature profile
    if plot_fits:
        fig, ax = plt.subplots()
        ax.plot(ths[zs <= max_z_fit], z_ci / 1.0e3, "b", label="Data")
        ax.plot(ci_estimate["thfit"], z_ci / 1.0e3, "--k", label="RZ fit", zorder=-1)
        ax.set_xlim([285.0, 312.0])
        ax.set_ylim([0.0, 4.0])
        ax.set_ylabel("$z$ [km]")
        ax.set_xlabel("$\\theta$ [K]")
        plt.legend()
        plt.tight_layout()
        plt.show()
    # CI altitudes
    inv_bottom = ci_estimate["h0"]
    H = ci_estimate["h1"]
    inv_top = ci_estimate["h2"]
    # Determine reference potential temperature
    th0 = np.interp(H, zs, ths)
    # Inversion strength
    if ci_estimate["a"] <= 0.2 or ci_estimate["a"] <= 2 * ci_estimate["b"]:
        # No inversion strength in the following cases:
        # a<=0.2: encroachment (No inversion layer, so the entire profile is given by g and a=0
        #           (considered a,0.2 as in paper))
        # a<=2*b: inversion lapse rate is equal to or smaller than free lapse rate
        dth = 0.0
    else:
        dth = ci_estimate["dth"]
    # Lapse rate
    dthdz = ci_estimate["gamma"]
    return inv_bottom, H, inv_top, th0, dth, dthdz


def read_turbine_type(turb_dat):
    # Turbine geometry
    hh = turb_dat["hub_height"]
    rd = turb_dat["rotor_diameter"]
    # Ct curve data
    ct = turb_dat["performance"]["Ct_curve"]["Ct_values"]
    ct_ws = turb_dat["performance"]["Ct_curve"]["Ct_wind_speeds"]
    # Cp curve data
    air_density = 1.225  # Hard-coded for now
    if "Cp_curve" in turb_dat["performance"]:
        # Read out Cp curve
        cp = turb_dat["performance"]["Cp_curve"]["Cp_values"]
        cp_ws = turb_dat["performance"]["Cp_curve"]["Cp_wind_speeds"]
    elif "power_curve" in turb_dat["performance"]:
        # Power curve data
        cp_ws = np.array(turb_dat["performance"]["power_curve"]["power_wind_speeds"])
        pows = np.array(turb_dat["performance"]["power_curve"]["power_values"])
        # Filter out Nan values and zero wind speeds
        selection = np.logical_and(
            np.greater(cp_ws, 0.0),
            np.logical_not(np.logical_or(np.isnan(cp_ws), np.isnan(pows))),
        )
        cp_ws = cp_ws[selection]
        pows = pows[selection]
        # Convert power curve to Cp curve
        rotor_area = np.pi * (rd / 2) ** 2
        cp = np.divide(
            np.array(pows), 0.5 * air_density * np.array(cp_ws) ** 3 * rotor_area
        )
    else:
        raise Exception("Bad Power Curve")
    # Ct and Cp curves
    ct_curve = interp1d(ct_ws, ct, fill_value="extrapolate")
    cp_curve = interp1d(cp_ws, cp, fill_value="extrapolate")
    return hh, rd, ct_curve, cp_curve


def wf_setup(farm_dat, analysis_dat, L_filter=1.0e3, debug_mode=False, rotation=0.0):
    # WAYVE imports
    from wayve.forcing.apm_forcing import ForcingComposite
    from wayve.forcing.wind_farms.dispersive_stresses import DispersiveStresses
    from wayve.forcing.wind_farms.entrainment import ConstantFlux
    from wayve.forcing.wind_farms.wind_farm import Turbine, WindFarm

    ####################
    # Set up WindFarm object
    ####################
    # Get x and y positions. Copy: this is called once per state now (the
    # solver-frame rotation differs per state), so farm_dat must stay pristine.
    x = np.asarray(farm_dat["layouts"][0]["coordinates"]["x"], dtype=float)
    y = np.asarray(farm_dat["layouts"][0]["coordinates"]["y"], dtype=float)
    # Reposition to be at grid center
    wf_offset_x = np.mean(x)
    wf_offset_y = np.mean(y)
    x = x - wf_offset_x
    y = y - wf_offset_y
    # Rotate the layout into the solver frame (see flow_io_abl: the ABL is
    # built with the hub-height wind along +x, wayve's westerly convention,
    # so the layout must turn with it to keep the wind-relative geometry).
    if rotation != 0.0:
        c, s = np.cos(rotation), np.sin(rotation)
        x, y = c * x + s * y, c * y - s * x
    # Number of turbines
    Nt = len(x)
    # Get turbine types
    turb_types = {}
    if "turbines" in farm_dat:
        type_inds = [0 for _ in range(Nt)]
        hh, rd, ct_curve, cp_curve = read_turbine_type(farm_dat["turbines"])
        turb_types[0] = [hh, rd, ct_curve, cp_curve]
    else:
        type_inds = farm_dat["layouts"][0]["turbine_types"]
        for i in np.unique(type_inds):
            hh, rd, ct_curve, cp_curve = read_turbine_type(farm_dat["turbine_types"][i])
            turb_types[i] = [hh, rd, ct_curve, cp_curve]
    # Turbine setup
    turbines = []
    for t in range(Nt):
        hh, rd, ct_curve, cp_curve = turb_types[type_inds[t]]
        turbine = Turbine(x[t], y[t], rd, hh, ct_curve, cp_curve)
        turbines.append(turbine)
    # Set up wake model
    wake_model = wake_model_setup(analysis_dat, debug_mode)
    # Set up coupling object
    coupling = wm_coupling_setup(analysis_dat, wake_model)
    # Generate wind farm object
    wind_farm = WindFarm(turbines, L_filter, coupling)
    # Combined forcing object
    forcing = ForcingComposite([wind_farm])
    # Additional forcing components
    if "APM_additional_terms" in analysis_dat:
        apm_terms_dat = analysis_dat["APM_additional_terms"]
        if "apm_disp_stresses" in apm_terms_dat:
            if apm_terms_dat["apm_disp_stresses"]["ds_type"] == "subgrid":
                if wind_farm.coupling.wm_velocity_handler is None:
                    raise ValueError(
                        "Subgrid dispersive stresses parametrization requires a subgrid to be included"
                    )
                disp_stresses = DispersiveStresses(wind_farm)
                forcing.add_child(disp_stresses)
        if (
            "momentum_entrainment" in apm_terms_dat
            and apm_terms_dat["momentum_entrainment"]["mfp_type"] == "constant_flux"
            and wind_farm.area > 0.0
        ):
            mfp_dat = apm_terms_dat["momentum_entrainment"]
            a_tau = 0.12
            if "a_mfp" in mfp_dat["apm_mfp_settings"]:
                a_tau = mfp_dat["apm_mfp_settings"]["a_mfp"]
            d_tau = 27.8
            if "d_mfp" in mfp_dat["apm_mfp_settings"]:
                d_tau = mfp_dat["apm_mfp_settings"]["d_mfp"]
            mfp = ConstantFlux(wind_farm, a=a_tau, d=d_tau)
            forcing.add_child(mfp)
    return wind_farm, forcing, wf_offset_x, wf_offset_y


def wm_coupling_setup(analysis_dat, wake_model):
    # WAYVE imports
    from wayve.forcing.wind_farms.wake_model_coupling.coupling_methods.pressure_based import (
        PressureBased,
    )
    from wayve.forcing.wind_farms.wake_model_coupling.coupling_methods.upstream import (
        Upstream,
    )
    from wayve.forcing.wind_farms.wake_model_coupling.coupling_methods.varying_background import (
        SelfSimilarWMVH,
        WakeModelVelocityHandler,
    )
    from wayve.forcing.wind_farms.wake_model_coupling.coupling_methods.velocity_matching import (
        VelocityMatching,
    )

    # Read inputs
    wmc_dat = analysis_dat["wm_coupling"]
    # Subgrid settings
    wm_velocity_handler = None
    if "subgrid" in wmc_dat and wmc_dat["subgrid"]["include_subgrid"]:
        sg_ratio = wmc_dat["subgrid"]["D_to_dx"]
        if analysis_dat["superposition_model"]["ws_superposition"] == "Product":
            wm_velocity_handler = SelfSimilarWMVH(sg_ratio)
        else:
            wm_velocity_handler = WakeModelVelocityHandler(sg_ratio)
    # Wake model coupling settings
    if "method" not in wmc_dat or wmc_dat["method"] == "PB":
        # Use pressure-based method
        coupling = PressureBased(wake_model, wm_velocity_handler)
    elif wmc_dat["method"] == "VM":
        if analysis_dat["superposition_model"]["ws_superposition"] != "Product":
            raise ValueError("VM method requires product-based superposition")
        # Read settings
        alpha = wmc_dat["settings"]["alpha"]
        # Use velocity matching method
        coupling = VelocityMatching(wake_model, wm_velocity_handler, alpha)
    elif wmc_dat["method"] == "US":
        # Read settings
        distance = wmc_dat["settings"]["distance"]
        # Use velocity matching method
        coupling = Upstream(wake_model, wm_velocity_handler, distance)
    else:
        raise ValueError("Wake model coupling not implemented!")
    return coupling


def wake_model_setup(analysis_dat, debug_mode=False):
    from wayve.forcing.wind_farms.wake_model_coupling.wake_models.lanzilao_merging import (
        Lanzilao,
    )

    # WM tool. The windIO schema places `wake_tool` inside `wm_coupling`; a
    # top-level `analysis.wake_tool` is rejected by validation (the schema sets
    # additionalProperties: false). Older hand-written yamls put it at the
    # analysis level, so keep reading that as a fallback.
    wake_tool = analysis_dat.get("wm_coupling", {}).get(
        "wake_tool", analysis_dat.get("wake_tool", "wayve")
    )
    if wake_tool == "wayve":
        # Read wake model settings #
        wm_dat = analysis_dat["wind_deficit_model"]
        k_dat = wm_dat["wake_expansion_coefficient"]
        # windIO convention: k = k_a + k_b * TI. wayve's Lanzilao computes
        # kwake = ka * TI + kb, so the pair maps swapped: ka=k_b, kb=k_a.
        # A scalar k is a constant expansion -> kb, with no TI term.
        if "k_a" in k_dat and "k_b" in k_dat and "ceps" in wm_dat:
            ti_coef = k_dat["k_b"]
            k_const = k_dat["k_a"]
            ceps = wm_dat["ceps"]
        elif "k" in k_dat and "ceps" in wm_dat:
            ti_coef = 0.0
            k_const = k_dat["k"]
            ceps = wm_dat["ceps"]
        else:
            raise ValueError("Wake spreading parameter not specified!")
        # Use wake merging method of Lanzilao and Meyers (2021)
        wake_model = Lanzilao(ka=ti_coef, kb=k_const, eps_beta=ceps)
    elif wake_tool == "foxes":
        require("foxes")
        from foxes import ModelBook
        from foxes.input.yaml.windio.read_attributes import _read_analysis
        from foxes.utils import Dict
        from wayve.couplings.foxes_coupling import FoxesWakeModel

        verbosity = 1 if debug_mode else 0

        # foxes' Dict takes its own label as `_name`; a plain `name=` kwarg is
        # stored as a *data* key and would be forwarded into the foxes
        # Algorithm constructor (TypeError: unexpected keyword argument 'name').
        algo_dict = Dict(
            algo_type="Downwind",
            wake_models=[],
            verbosity=verbosity,
            _name="wayve.algorithm",
        )

        ana_dict = Dict(analysis_dat, _name="analysis")
        idict = Dict(algorithm=algo_dict, _name="wayve")
        mbook = ModelBook()

        _read_analysis(ana_dict, idict, mbook=mbook, verbosity=verbosity)

        wake_model = FoxesWakeModel(mbook=mbook, **idict["algorithm"])
        # wayve@e87780a overrides background_flow_direction with logic
        # identical to the base class but a broken lazy import (e_spanwise
        # from wake_models.wake_model_tools instead of forcing_tools), which
        # crashes every foxes-coupled solve. Rebind the base implementation
        # while the override still carries the broken import; probing the
        # override's own source (rather than the import target) keeps the
        # shim inert once upstream fixes or rewrites the method.
        import inspect

        try:
            bfd_src = inspect.getsource(FoxesWakeModel.background_flow_direction)
        except (OSError, TypeError):
            bfd_src = ""
        if "wake_model_tools import e_spanwise" in bfd_src:
            import types

            from wayve.forcing.wind_farms.wake_model_coupling.wake_model_interface import (
                UniDirectionalSelfSimilar,
            )

            wake_model.background_flow_direction = types.MethodType(
                UniDirectionalSelfSimilar.background_flow_direction, wake_model
            )
    else:
        raise NotImplementedError(f"Wake tool '{wake_tool}' not implemented!")
    return wake_model


# Capping inversion: the windIO energy-resource schema spells the four scalars
# flat, and WIFA's own code_saturne adapter already reads them that way
# (cs_launch_modules.py). The nested
# ``thermal_stratification.capping_inversion`` block is the older spelling and
# is still accepted, so files written against it keep working; it cannot be
# expressed in a netCDF ``!include``, whose root group is flat.
# Ordered (h, dH, dtheta, lapse_rate) — the tuple read_capping_inversion returns.
_CI_FLAT_KEYS = (
    "ABL_height",
    "capping_inversion_thickness",
    "capping_inversion_strength",
    "lapse_rate",
)
_CI_NESTED_KEYS = ("ABL_height", "dH", "dtheta", "lapse_rate")


def capping_inversion_spelling(wind_resource_dat):
    """Which spelling this wind resource uses for the capping inversion.

    Parameters
    ----------
    wind_resource_dat: dict
        Wind resource data

    Returns
    -------
    str or None
        ``"flat"``, ``"nested"``, or None when no inversion is specified (the
        caller then fits one from the potential-temperature profile, or falls
        back to defaults).

    Raises
    ------
    ValueError
        If both spellings are present, or if either is incomplete. Both are
        defects of the *file*, not of a single state, which is why this check
        is separate from the per-state value read: ``run_wayve`` calls it once
        before the state loop, whose ``except Exception`` would otherwise turn
        a malformed file into a run where every state silently "crashed".

        A partial block is rejected rather than completed from defaults: of the
        four scalars, ``capping_inversion_strength`` alone sets the amplitude of
        the gravity-wave forcing, so silently pairing a real ``ABL_height`` with
        an invented dtheta yields a result that looks fitted and is not.
    """
    flat = [key for key in _CI_FLAT_KEYS if key in wind_resource_dat]
    nested = wind_resource_dat.get("thermal_stratification", {}).get(
        "capping_inversion"
    )
    if flat and nested is not None:
        raise ValueError(
            "Wind resource specifies the capping inversion twice: flat keys "
            f"{sorted(flat)} and a nested thermal_stratification."
            "capping_inversion block. Keep one (the flat keys are the windIO "
            "schema spelling)."
        )
    if nested is not None:
        missing = [key for key in _CI_NESTED_KEYS if key not in nested]
        if missing:
            raise ValueError(
                "Incomplete thermal_stratification.capping_inversion block; "
                f"missing {missing}. All of {list(_CI_NESTED_KEYS)} are required."
            )
        return "nested"
    if not flat:
        return None
    if len(flat) != len(_CI_FLAT_KEYS):
        raise ValueError(
            f"Incomplete capping inversion: found {sorted(flat)}, missing "
            f"{sorted(set(_CI_FLAT_KEYS) - set(flat))}. All four are required "
            "(they are not completed from defaults; see "
            "capping_inversion_spelling). Omit all four to fit the inversion "
            "from the potential-temperature profile instead."
        )
    return "flat"


def read_capping_inversion(wind_resource_dat, time_index):
    """The capping inversion of one state, in either windIO spelling.

    Parameters
    ----------
    wind_resource_dat: dict
        Wind resource data
    time_index: int
        Index of the timestamp to read

    Returns
    -------
    tuple or None
        ``(h, dh, dth, dthdz)`` — inversion centre height [m], thickness [m],
        strength [K] and free-atmosphere lapse rate [K/m] — or None when the
        resource specifies no inversion.
    """
    spelling = capping_inversion_spelling(wind_resource_dat)
    if spelling is None:
        return None
    if spelling == "nested":
        block = wind_resource_dat["thermal_stratification"]["capping_inversion"]
        return tuple(float(block[key]["data"][time_index]) for key in _CI_NESTED_KEYS)
    return tuple(
        float(wind_resource_dat[key]["data"][time_index]) for key in _CI_FLAT_KEYS
    )


def flow_io_abl(
    wind_resource_dat, time_index, zh, h1, dh_max=None, serz=True, gmode="avg"
):
    """
    Method to set up an ABL object based on FLOW IO

    Parameters
    ----------
    wind_resource_dat: dict
        Wind resource data
    time_index: int
        Index of the timestamp to set up ABL for
    zh: float
        Mean turbine hub height
    h1: float
        Lower layer height
    dh_max (optional): float
        Maximum depth of the inversion layer used in the inversion curve fitting procedure (default: None)
    serz (optional): boolean
        Whether the surface-extended version of the RZ model is used (default: True)
    gmode (optional): str
        How the free-atmosphere (geostrophic) velocity is derived from
        height-resolved wind resources with mesoscale surface scalars
        (see wayve's ``abl_setup.mesoscale_based``): "h1" (at the inversion
        center), "h2" (at the inversion top), "avg" (profile average between
        the inversion and 5 km; needs profiles reaching 5 km), or "trop"
        (average to the tropopause; needs profiles reaching the stratosphere).
        Ignored by the scalar and turbulence-profile input paths.
        (default: "avg")

    Returns
    -------
    abl: wayve.abl.abl.ABL
        Atmospheric state in the *solver frame*: the hub-height wind blows
        along +x (wayve's westerly convention, which its gravity-wave
        machinery and anisotropic grids assume). Vertical veer is preserved
        as spanwise components.
    rotation: float
        Angle (radians, counterclockwise) of the earth-frame hub-height flow
        vector measured from east. Rotating earth-frame coordinates by
        ``-rotation`` maps them into the solver frame (``wf_setup`` does this
        to the layout); rotating solver-frame vectors by ``+rotation`` maps
        results back to earth.
    """
    # Atmospheric state setup
    from wayve.abl.abl import ABL

    # Constants #
    gravity = 9.80665  # [m s-2]
    kappa = 0.41  # Von Karman constant
    omega = 7.2921159e-5  # angular speed of the Earth [rad/s]
    # Basic atmospheric scalars #
    # Operating air density. The APM velocity solution does not depend on it
    # (see wayve's power_turbines docstring); it linearly scales the reported
    # turbine power — the Cp curve itself stays referenced to the standard
    # 1.225 kg/m3 of the windIO power curve (read_turbine_type) — and is
    # forwarded to foxes as FV.RHO by the foxes coupling.
    air_density = 1.225
    if "air_density" in wind_resource_dat.keys():
        air_density = wind_resource_dat["air_density"]["data"][time_index]
    # Surface roughness
    z0 = 1.0e-1
    if "z0" in wind_resource_dat.keys():
        z0 = wind_resource_dat["z0"]["data"][time_index]
    # Monin-Obhukov length
    l_mo = 5.0e3
    if "LMO" in wind_resource_dat.keys():
        l_mo = wind_resource_dat["LMO"]["data"][time_index]
    # Coriolis parameter #
    phi = 0.377  # Assume latitude location
    fc = 2 * omega * np.sin(phi)
    if "fc" in wind_resource_dat.keys():
        fc = wind_resource_dat["fc"]["data"][time_index]
    # Check if wind resource contains vertical profile
    profile_input = "height" in wind_resource_dat.keys()
    if not profile_input:
        # Wind speed and direction
        v = wind_resource_dat["wind_speed"]["data"][time_index]
        wd = wind_resource_dat["wind_direction"]["data"][time_index]
        # Solver-frame normalization: build the profile as if the wind were
        # westerly (the Ekman veer relative to the hub-height direction is
        # unchanged) and report the rotation undone by doing so.
        rotation = np.deg2rad(270.0 - wd)
        wd = 270.0
        # Friction velocity
        ust = 0.666
        if "friction_velocity" in wind_resource_dat.keys():
            ust = wind_resource_dat["friction_velocity"]["data"][time_index]
        # Turbulence intensity. windIO carries TI as a fraction (as does the
        # vertical-profile branch below, and the 0.04 default just above), so
        # it is used as-is. Guard on the variable actually being read.
        TI = 0.04
        if "turbulence_intensity" in wind_resource_dat.keys():
            TI = wind_resource_dat["turbulence_intensity"]["data"][time_index]
        # Capping inversion information
        h = 1.5e3
        dh = 100.0
        dth = 5.0
        dthdz = 2.0e-3
        th0 = 293.15
        ci = read_capping_inversion(wind_resource_dat, time_index)
        if ci is not None:
            h, dh, dth, dthdz = ci
        inv_bottom, inv_top = h - dh / 2, h + dh / 2
        # Nieuwstadt profiles for velocity and shear stress
        zs, us, vs, U3, V3, tauxs, tauys, nus = nieuwstadt83_profiles(
            zh, v, wd, z0=z0, h=h, ust=ust, fc=fc
        )
        # Potential temperature profile constant
        ths = th0 * np.ones_like(zs)
    else:
        # Read out vertical profile
        zs = np.array(wind_resource_dat["height"])
        if not np.all(np.diff(zs) > 0):
            # np.interp silently returns garbage on unsorted abscissae, and
            # the hub interpolation below now derives the whole solver frame.
            raise UserWarning("wind resource 'height' must be strictly ascending")
        vs = np.array(wind_resource_dat["wind_speed"]["data"][time_index])
        wds = np.array(wind_resource_dat["wind_direction"]["data"][time_index])
        ths = np.array(wind_resource_dat["potential_temperature"]["data"][time_index])
        # Solver-frame normalization: shift the direction profile so the
        # hub-height wind is westerly (flow along +x). Shifting the input
        # keeps the veer and makes everything derived below (velocity
        # components, Gmode geostrophic wind, momentum-flux alignment) land
        # in the solver frame consistently — without wayve's ABL.rotate(),
        # whose momentum-flux rotation is broken in wayve 2.0.0. Unwrap
        # first so interpolation across the 0/360 seam is safe.
        wds = np.rad2deg(np.unwrap(np.deg2rad(wds)))
        wd_hub = np.interp(zh, zs, wds)
        rotation = np.deg2rad(270.0 - wd_hub)
        wds = wds + (270.0 - wd_hub)
        # Turbulence intensity: height-resolved profile or one value per state
        TIs = np.atleast_1d(
            np.array(wind_resource_dat["turbulence_intensity"]["data"][time_index])
        )
        TI = np.interp(zh, zs, TIs) if TIs.size > 1 else float(TIs[0])
        # Velocity components
        us = -vs * np.sin(np.deg2rad(wds))
        vs = -vs * np.cos(np.deg2rad(wds))
        # Check available inputs
        if "k" in wind_resource_dat.keys() or "tau_x" in wind_resource_dat.keys():
            # Turbulence profiles provided directly (LES/RANS-like inputs)
            if "k" in wind_resource_dat.keys():  # RANS-like inputs
                tkes = np.array(wind_resource_dat["k"]["data"][time_index])
                eps = np.array(wind_resource_dat["epsilon"]["data"][time_index])
                # Eddy viscosity
                C_mu = 0.09  # k-epsilon model value
                nus = C_mu * np.divide(
                    np.square(tkes), eps, out=np.zeros_like(tkes), where=eps != 0
                )
                # Momentum fluxes
                dudz = np.gradient(us, zs, edge_order=2)
                dvdz = np.gradient(vs, zs, edge_order=2)
                tauxs = nus * dudz
                tauys = nus * dvdz
            else:  # Shear stress profile directly available
                tauxs = np.array(wind_resource_dat["tau_x"]["data"][time_index])
                tauys = np.array(wind_resource_dat["tau_y"]["data"][time_index])
                # Explicit stress components are earth-frame vectors: rotate
                # them into the solver frame along with the velocities.
                c_r, s_r = np.cos(rotation), np.sin(rotation)
                tauxs, tauys = c_r * tauxs + s_r * tauys, c_r * tauys - s_r * tauxs
                nus = None
            # Total momentum flux
            taus = np.sqrt(np.square(tauxs) + np.square(tauys))
            # Friction velocity
            ust = taus[0]  # Assume friction velocity is not given explicitly
            # Estimate boundary layer height based on momentum flux #
            f_tau = interp1d(taus, zs)
            blh = f_tau(0.01 * ust)
            # Capping inversion information
            ci = read_capping_inversion(wind_resource_dat, time_index)
            if ci is not None:
                h, dh, dth, dthdz = ci
                # Mixed-layer temperature from the profile the resource
                # supplies, as every other branch does. A fixed 293.15 K here
                # made the buoyancy jump g*dtheta/th0 disagree with the ABL's
                # own theta profile by up to ~4% at realistic temperatures.
                th0 = np.interp(h, zs, ths)
                inv_bottom, inv_top = h - dh / 2, h + dh / 2
            else:
                inv_bottom, h, inv_top, th0, dth, dthdz = ci_fitting(
                    zs, ths, l_mo, blh, dh_max=dh_max, serz=serz
                )
            # Geostrophic wind speed
            z = np.linspace(h, 15.0e3, 1000)
            _trapezoid = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
            U3 = _trapezoid(np.interp(z, zs, us), z) / (15.0e3 - h)
            V3 = _trapezoid(np.interp(z, zs, vs), z) / (15.0e3 - h)
        elif (
            "friction_velocity" in wind_resource_dat.keys()
            and "boundary_layer_height" in wind_resource_dat.keys()
        ):
            # Mesoscale/reanalysis-style inputs (e.g. ERA5): u/v/theta profiles
            # plus surface scalars (ust, blh, L_MO) instead of turbulence
            # profiles — wayve's abl_setup.mesoscale_based() territory.
            ust = wind_resource_dat["friction_velocity"]["data"][time_index]
            blh = wind_resource_dat["boundary_layer_height"]["data"][time_index]
            if zs[-1] >= 12.0e3:
                # Profiles reach the stratosphere: use wayve's full setup,
                # including the tropopause two-line fit. Note the radiating
                # top-layer state differs from the truncated branch below:
                # here Uinf/Vinf are means above the fitted tropopause and
                # Ninf is the fitted stratospheric N, while the truncated
                # branch uses the top-of-profile point and the inversion
                # fit's free lapse rate.
                from wayve.abl.abl_setup import mesoscale_based

                lat = np.rad2deg(np.arcsin(fc / (2.0 * omega)))
                # us/vs are already solver-frame (see above), so the ABL
                # comes out westerly-aligned without further rotation.
                return (
                    mesoscale_based(
                        zs,
                        us,
                        vs,
                        ths,
                        ust,
                        blh,
                        l_mo,
                        lat,
                        h1,
                        z0=(z0 if "z0" in wind_resource_dat.keys() else None),
                        TI=TI,
                        rho=air_density,
                        dh_max=(300.0 if dh_max is None else dh_max),
                        Gmode=gmode,
                        serz=serz,
                    ),
                    rotation,
                )
            # Truncated profiles (reanalysis subsets often stop below the
            # tropopause, e.g. ERA5 levels 96-137 end near 6 km): run the same
            # core setup but skip mesoscale_based's unconditional tropopause
            # fit, which would feed a garbage stratosphere into the ABL.  The
            # free atmosphere is capped at the profile top instead (see the
            # ABL constructor below), and Gmode is restricted to what the
            # data supports.
            stable = 0.0 < l_mo < 100
            # Capping inversion: an explicit windIO block wins over the fit
            ci = read_capping_inversion(wind_resource_dat, time_index)
            if ci is not None:
                h, dh, dth, dthdz = ci
                inv_bottom, inv_top = h - dh / 2, h + dh / 2
                th0 = np.interp(h, zs, ths)
            else:
                inv_bottom, h, inv_top, th0, dth, dthdz = ci_fitting(
                    zs,
                    ths,
                    l_mo,
                    blh,
                    dh_max=(300.0 if dh_max is None else dh_max),
                    serz=serz,
                )
            # Momentum flux and eddy viscosity profiles from surface scalars
            # (same closure as mesoscale_based: stable turbulence extends to
            # blh, convective turbulence to the inversion)
            tau = np.zeros(zs.shape)
            nus = np.zeros(zs.shape)
            if stable:
                m = zs <= blh
                tau[m] = ust**2 * (1 - zs[m] / blh) ** 1.5
                nus[m] = kappa * ust * zs[m] * (1 - zs[m] / blh) ** 2
            else:
                m = zs <= h
                tau[m] = ust**2 * (1 - zs[m] / h)
                nus[m] = kappa * ust * zs[m] * (1 - zs[m] / h) ** 2
            # Geostrophic wind from the profile, restricted to the data range
            if gmode == "trop":
                raise UserWarning(
                    "Gmode 'trop' needs profiles reaching the stratosphere; "
                    f"these stop at {zs[-1]:.0f} m"
                )
            if gmode == "avg" and zs[-1] < 5.0e3:
                warnings.warn(
                    f"Gmode 'avg' averages the wind up to 5 km but the "
                    f"profiles stop at {zs[-1]:.0f} m; falling back to 'h2' "
                    "(wind at the inversion top)"
                )
                gmode = "h2"
            if gmode == "h1":
                U3 = np.interp(h, zs, us)
                V3 = np.interp(h, zs, vs)
            elif gmode == "h2":
                U3 = np.interp(inv_top, zs, us)
                V3 = np.interp(inv_top, zs, vs)
            elif gmode == "avg":
                z = np.linspace(h, 5.0e3, 1000)
                _trapezoid = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
                U3 = _trapezoid(np.interp(z, zs, us), z) / (5.0e3 - h)
                V3 = _trapezoid(np.interp(z, zs, vs), z) / (5.0e3 - h)
            else:
                raise UserWarning(f"Gmode '{gmode}' unknown")
            # Momentum flux aligned with the geostrophic wind (as in
            # mesoscale_based)
            tau_angle = np.arctan2(V3, U3)
            tauxs = np.cos(tau_angle) * tau
            tauys = np.sin(tau_angle) * tau
            # Log-law surface roughness estimate when not provided
            if "z0" not in wind_resource_dat.keys():
                z0 = zs[0] / np.exp(kappa * np.hypot(us[0], vs[0]) / ust)
        else:
            raise UserWarning(
                "Vertical-profile wind resource needs either turbulence "
                "profiles ('k'/'epsilon' or 'tau_x'/'tau_y') or mesoscale "
                "surface scalars ('friction_velocity' and "
                "'boundary_layer_height')"
            )
    # Upper layer thickness
    h2 = h - h1
    if (
        inv_bottom <= h1 + 10.0
    ):  # H cannot be lower than H1 and the upper layer must be at least 10m
        raise RuntimeWarning(f"CI too low, CI bottom located at z={int(inv_bottom)}m")
    # CI check
    if dth == 0.0:
        raise RuntimeWarning("No CI present!")
    # gprime and N
    gprime = gravity * dth / th0
    N = np.sqrt(gravity * dthdz / th0)
    # Set up ABL object. The free atmosphere is capped at the top of the
    # profile data: NonUniform's layer setup spans [inversion top, h_strat],
    # and the default h_strat of 10 km would put most of its layers in a
    # no-data region where the profile splines just hold constants. Above
    # the cap sits the semi-infinite radiating layer with the top-of-profile
    # state (Uinf/Vinf) and the fitted free-atmosphere stratification (Ninf).
    return (
        ABL(
            zs,
            us,
            vs,
            ths,
            tauxs,
            tauys,
            h1,
            h2,
            gprime,
            N,
            U3,
            V3,
            fc,
            nus=nus,
            rho=air_density,
            TI=TI,
            z0=z0,
            ust=ust,
            inv_bottom=inv_bottom,
            inv_top=inv_top,
            h_strat=min(10.0e3, zs[-1]),
            Uinf=us[-1],
            Vinf=vs[-1],
            Ninf=N,
        ),
        rotation,
    )


def run():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_yaml", help="The input yaml file")
    args = parser.parse_args()

    run_wayve(args.input_yaml)


if __name__ == "__main__":
    run()
