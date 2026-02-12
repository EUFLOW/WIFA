# General packages
import argparse
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import mpmath
import numpy as np
import xarray as xr
from scipy.interpolate import interp1d
from scipy.special import gamma as scipy_gamma
from windIO import load_yaml


def run_wayve(yamlFile, output_dir="output", debug_mode=False):
    # General APM setup
    from wayve.apm import APM
    from wayve.grid.grid import Stat2Dgrid
    from wayve.momentum_flux_parametrizations import FrictionCoefficients
    from wayve.pressure.gravity_waves.gravity_waves import NonUniform, Uniform
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

    ##################
    # Other APM components
    ##################
    # Momentum flux parametrization
    mfp = FrictionCoefficients()
    # Pressure feedback parametrization
    pressure = Uniform(dynamic=True, rotating=False)
    if "layers_description" in analysis_dat:
        if "number_of_fa_layers" in analysis_dat["layers_description"]:
            n_layers = analysis_dat["layers_description"]["number_of_fa_layers"]
            if n_layers > 1:
                pressure = NonUniform(n_layers=n_layers, order=1)

    ######################
    # Read output settings
    ######################
    # Select timestamps
    times = resource_dat["wind_resource"]["time"]
    run_config = system_dat["attributes"]["model_outputs_specification"][
        "run_configuration"
    ]
    if "times_run" in run_config and not run_config["times_run"].get(
        "all_occurences", True
    ):
        if "subset" in run_config["times_run"]:
            subset = run_config["times_run"]["subset"]
            times = [times[i] for i in subset]
    # Get turbine variables to output
    turbine_nc_filename = "turbine_data.nc"
    turbine_output_variables = ["power", "rotor_effective_velocity"]
    if "turbine_outputs" in system_dat["attributes"]["model_outputs_specification"]:
        turb_out_dat = system_dat["attributes"]["model_outputs_specification"][
            "turbine_outputs"
        ]
        if "turbine_nc_filename" in turb_out_dat:
            turbine_nc_filename = turb_out_dat["turbine_nc_filename"]
        if "turbine_output_variables" in turb_out_dat:
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
    # Initialize crash counter
    crashes = 0
    # List of datasets
    ds_list = []
    ds_ff_list = []
    # Loop over timeseries
    for time_index, time in enumerate(times):
        if debug_mode:
            # Print timestep
            print(f"time {time_index+1}/{len(times)}")
        try:
            # Set up ABL
            abl = flow_io_abl(resource_dat["wind_resource"], time_index, hh, h1)
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
            # NC setup
            ds = xr.Dataset(
                turb_out_dict,
                coords={"states": time, "turbine": range(Nt)},
            )
            # Add to output list
            ds_list.append(ds)
            # Flow field outputs #
            if report_flow and not debug_mode:
                # Callables for flow evaluation
                u_bg_evaluator = coupling.set_up_u_bg_evaluator(
                    abl
                )  # Background velocity callable
                apm_evaluator = coupling.apm_evaluator  # APM lower layer state callable
                # Output arrays
                wind_speed = np.zeros([len(x_ff), len(y_ff), len(z_ff)])
                wind_dir = np.zeros([len(x_ff), len(y_ff), len(z_ff)])
                # Loop over z-planes
                for k, z_k in enumerate(z_ff):
                    # Get velocities
                    u_bg, v_bg, u_wm, v_wm = wake_model.xy_plane(
                        wind_farm,
                        abl,
                        u_bg_evaluator,
                        apm_evaluator,
                        x_ff - wf_offset_x,
                        y_ff - wf_offset_y,
                        z_k,
                    )
                    # Convert to speed and direction
                    wind_speed[:, :, k] = np.sqrt(np.square(u_wm) + np.square(v_wm))
                    wind_dir[:, :, k] = np.rad2deg(
                        np.pi / 2 - (np.arctan2(v_wm, u_wm) + np.pi)
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
                # Add to output list
                ds_ff_list.append(ds_ff)

        except Exception as exc:
            print(exc)
            # Update crash counter
            crashes += 1
            continue
    if debug_mode:
        print(f"crashes: {crashes}/{len(times)}")

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


def nieuwstadt83_profiles(zh, v, wd, z0=1.0e-1, h=1.5e3, fc=1.0e-4, ust=0.666):
    """Set up the cubic analytical profile from Nieuwstadt (1983), based on hub height velocity information"""
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


def wf_setup(farm_dat, analysis_dat, L_filter=1.0e3, debug_mode=False):
    # WAYVE imports
    from wayve.forcing.apm_forcing import ForcingComposite
    from wayve.forcing.wind_farms.dispersive_stresses import DispersiveStresses
    from wayve.forcing.wind_farms.entrainment import ConstantFlux
    from wayve.forcing.wind_farms.wind_farm import Turbine, WindFarm

    ####################
    # Set up WindFarm object
    ####################
    # Get x and y positions
    x = farm_dat["layouts"][0]["coordinates"]["x"]
    y = farm_dat["layouts"][0]["coordinates"]["y"]
    # Reposition to be at grid center
    wf_offset_x = np.mean(x)
    wf_offset_y = np.mean(y)
    x -= wf_offset_x
    y -= wf_offset_y
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
    # WAYVE imports
    from wayve.couplings.foxes_coupling import FoxesWakeModel
    from wayve.forcing.wind_farms.wake_model_coupling.wake_models.lanzilao_merging import (
        Lanzilao,
    )

    # WM tool
    wake_tool = analysis_dat.get(
        "wake_tool", "wayve"
    )  # updated by Jonas -TODO update this according to updated schema
    if wake_tool == "wayve":
        # Read wake model settings #
        wm_dat = analysis_dat["wind_deficit_model"]
        k_dat = wm_dat["wake_expansion_coefficient"]
        # k, k_a, k_b, ceps
        if "k_a" in k_dat and "k_b" in k_dat and "ceps" in wm_dat:
            k_a = k_dat["k_a"]
            k_b = k_dat["k_b"]
            ceps = wm_dat["ceps"]
        elif "k" in k_dat and "ceps" in wm_dat:
            k_a = k_dat["k"]
            k_b = 0.0
            ceps = wm_dat["ceps"]
        else:
            raise ValueError("Wake spreading parameter not specified!")
        # Use wake merging method of Lanzilao and Meyers (2021)
        wake_model = Lanzilao(ka=k_a, kb=k_b, eps_beta=ceps)
    elif wake_tool == "foxes":
        from foxes import ModelBook
        from foxes.input.yaml.windio.read_attributes import _read_analysis
        from foxes.utils import Dict

        verbosity = 1 if debug_mode else 0

        algo_dict = Dict(
            algo_type="Downwind",
            wake_models=[],
            verbosity=verbosity,
            name="wayve.algorithm",
        )

        ana_dict = Dict(analysis_dat, name="analysis")
        idict = Dict(algorithm=algo_dict, name="wayve")
        mbook = ModelBook()

        _read_analysis(ana_dict, idict, mbook=mbook, verbosity=verbosity)

        wake_model = FoxesWakeModel(mbook=mbook, **idict["algorithm"])
    else:
        raise NotImplementedError(f"Wake tool '{wake_tool}' not implemented!")
    return wake_model


def flow_io_abl(wind_resource_dat, time_index, zh, h1, dh_max=None, serz=True):
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
    """
    # Atmospheric state setup
    from wayve.abl.abl import ABL

    # Constants #
    gravity = 9.80665  # [m s-2]
    kappa = 0.41  # Von Karman constant
    omega = 7.2921159e-5  # angular speed of the Earth [rad/s]
    # Basic atmospheric scalars #
    air_density = 1.225  # Hard-coded for now
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
        # Friction velocity
        ust = 0.666
        if "friction_velocity" in wind_resource_dat.keys():
            ust = wind_resource_dat["friction_velocity"]["data"][time_index]
        # Turbulence intensity
        TI = 0.04
        if "z0" in wind_resource_dat.keys():
            TI = wind_resource_dat["turbulence_intensity"]["data"][time_index] / 100.0
        # Capping inversion information
        h = 1.5e3
        dh = 100.0
        dth = 5.0
        dthdz = 2.0e-3
        th0 = 293.15
        if "thermal_stratification" in wind_resource_dat.keys():
            thermal_data = wind_resource_dat["thermal_stratification"]
            if "capping_inversion" in thermal_data.keys():
                ci_data = thermal_data["capping_inversion"]
                h = ci_data["ABL_height"]["data"][time_index]
                dh = ci_data["dH"]["data"][time_index]
                dth = ci_data["dtheta"]["data"][time_index]
                dthdz = ci_data["lapse_rate"]["data"][time_index]
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
        vs = np.array(wind_resource_dat["wind_speed"]["data"][time_index])
        wds = np.array(wind_resource_dat["wind_direction"]["data"][time_index])
        ths = np.array(wind_resource_dat["potential_temperature"]["data"][time_index])
        TIs = np.array(wind_resource_dat["turbulence_intensity"]["data"][time_index])
        # Interpolate TI
        TI = np.interp(zh, zs, TIs)
        # Velocity components
        us = -vs * np.sin(np.deg2rad(wds))
        vs = -vs * np.cos(np.deg2rad(wds))
        # Check available inputs
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
            nus = None
        # Total momentum flux
        taus = np.sqrt(np.square(tauxs) + np.square(tauys))
        # Friction velocity
        ust = taus[0]  # Assume friction velocity is not given explicitly
        # Estimate boundary layer height based on momentum flux #
        f_tau = interp1d(taus, zs)
        blh = f_tau(0.01 * ust)
        # Capping inversion information
        if (
            "thermal_stratification" in wind_resource_dat.keys()
            and "capping_inversion"
            in wind_resource_dat["thermal_stratification"].keys()
        ):
            thermal_data = wind_resource_dat["thermal_stratification"]
            ci_data = thermal_data["capping_inversion"]
            th0 = 293.15
            h = ci_data["ABL_height"]["data"][time_index]
            dh = ci_data["dH"]["data"][time_index]
            dth = ci_data["dtheta"]["data"][time_index]
            dthdz = ci_data["lapse_rate"]["data"][time_index]
            inv_bottom, inv_top = h - dh / 2, h + dh / 2
        else:
            inv_bottom, h, inv_top, th0, dth, dthdz = ci_fitting(
                zs, ths, l_mo, blh, dh_max=dh_max, serz=serz
            )
        # Geostrophic wind speed
        z = np.linspace(h, 15.0e3, 1000)
        U3 = np.trapz(np.interp(z, zs, us), z) / (15.0e3 - h)
        V3 = np.trapz(np.interp(z, zs, vs), z) / (15.0e3 - h)
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
    # Set up ABL object
    return ABL(
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
    )


def run():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_yaml", help="The input yaml file")
    args = parser.parse_args()

    run_wayve(args.input_yaml)


if __name__ == "__main__":
    run()
