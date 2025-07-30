/*============================================================================
 * User functions for input of calculation parameters.
 *============================================================================*/

/* VERS */

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

  This program is free software; you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free Software
  Foundation; either version 2 of the License, or (at your option) any later
  version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
  details.

  You should have received a copy of the GNU General Public License along with
  this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
  Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/*----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>
#include <string.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_coupling.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Select physical model options, including user fields.
 *
 * This function is called at the earliest stages of the data setup,
 * so field ids are not available yet.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_model(void)
{

  cs_wall_functions_t *wf = cs_get_glob_wall_functions();
  wf->iwalfs = CS_WALL_F_S_MONIN_OBUKHOV;
  wf->iwallf = CS_WALL_F_2SCALES_SMOOTH_ROUGH;

  // cs_glob_atmo_option->open_bcs_treatment = 1; //was 0

  /* Read the meteo file (1) or LMO profile with following parameters (2) */
  cs_glob_atmo_option->meteo_profile = cs_notebook_parameter_value_by_name("meteo_profile");
  if (cs_glob_atmo_option->meteo_profile == 2)
  {
    /* Elevation for reference velocity */
    cs_glob_atmo_option->meteo_zref = cs_notebook_parameter_value_by_name("zref");
    /* Friction velocity */
    cs_glob_atmo_option->meteo_uref = cs_notebook_parameter_value_by_name("ureff");
  }

  // Needed for to compute the Coriolis force
  cs_glob_atmo_option->longitude = cs_notebook_parameter_value_by_name("long");
  cs_glob_atmo_option->latitude = cs_notebook_parameter_value_by_name("lat");
  /* Large scale roughness */
  cs_glob_atmo_option->meteo_z0 = cs_notebook_parameter_value_by_name("z0");
  /* Velocity direction */
  cs_glob_atmo_option->meteo_angle = cs_notebook_parameter_value_by_name("teta");

  /* Inverse LMO length */
  cs_glob_atmo_option->meteo_dlmo = cs_notebook_parameter_value_by_name("Lmoinv");
  /* Ground temperature */
  cs_glob_atmo_option->meteo_t0 = cs_notebook_parameter_value_by_name("t0");

  /* Post processing options */

  /* Output ground in the results file*/
  cs_glob_atmo_option->compute_z_ground = true;

  /* Added properties at the boundary */
  /* /\* To post process roughness read in tif file *\/ */
  /* cs_parameters_add_property("gdalroughness", */
  /*                            1, */
  /*                            CS_MESH_LOCATION_BOUNDARY_FACES); */

  /* cs_parameters_add_property("roughness_zone", */
  /*                            1, */
  /*                            CS_MESH_LOCATION_BOUNDARY_FACES); */

  /* To post-process k and epsilon */
  cs_parameters_add_property("tke_transport", 1, CS_MESH_LOCATION_CELLS);
  cs_parameters_add_property("eps_transport", 1, CS_MESH_LOCATION_CELLS);

  /* To post-process u* and uk */
  cs_parameters_add_property("boundary_ustar", 1, CS_MESH_LOCATION_CELLS);

  cs_parameters_add_property("boundary_uk", 1, CS_MESH_LOCATION_CELLS);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define or modify general numerical and physical user parameters.
 *
 * At the calling point of this function, most model-related most variables
 * and other fields have been defined, so specific settings related to those
 * fields may be set here.
 *
 * At this stage, the mesh is not built or read yet, so associated data
 * such as field values are not accessible yet, though pending mesh
 * operations and some fields may have been defined.
 *
 * \param[in, out]   domain    pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_parameters(cs_domain_t   *domain)
{
  CS_UNUSED(domain);

#if 0
    /* We only specify XYZ0 if we explicitely fix Dirichlet conditions
       for the pressure. */

    fp->xyzp0[0] = 0.;
    fp->xyzp0[1] = 0.;
    fp->xyzp0[2] = 0.;
#endif

  /* Example: activate mesh robustness options */
  /*-------------------------------------------*/

  cs_glob_mesh_quantities_flag |= CS_CELL_CENTER_CORRECTION;
  cs_glob_mesh_quantities_flag |= CS_CELL_FACE_CENTER_CORRECTION;
  cs_glob_mesh_quantities_flag |= CS_FACE_RECONSTRUCTION_CLIP;
  cs_glob_mesh_quantities_flag |= CS_CELL_VOLUME_RATIO_CORRECTION;
  cs_glob_mesh_quantities_flag |= CS_BAD_CELLS_WARPED_CORRECTION;
  //cs_glob_mesh_quantities_flag |= CS_BAD_CELLS_REGULARISATION;
  cs_glob_mesh_quantities_flag |= CS_FACE_DISTANCE_CLIP;

  // Uncomment following if global forcing term is applied
  cs_velocity_pressure_param_t *vp_param = cs_get_glob_velocity_pressure_param();
  vp_param->igpust = 0;

  cs_time_step_t *ts = cs_get_glob_time_step();
  ts->nt_max = cs_notebook_parameter_value_by_name("ntmax");

  /* Warning, meteo file does not overwrite reference values... */
  cs_fluid_properties_t *phys_pro = cs_get_glob_fluid_properties();
  // cs_real_t rair = phys_pro->r_pg_cnst;
  /* Reference fluid properties set from meteo values */
  phys_pro->p0 = cs_glob_atmo_option->meteo_psea;
  // phys_pro->t0 = cs_glob_atmo_option->meteo_t0; /* ref temp T0 */
  // phys_pro->ro0 = phys_pro->p0/(rair * cs_glob_atmo_option->meteo_t0); /* ref
  // density T0 */

  /* isstpc:
     0: slope test enabled
     1: slope test disabled (default)
     2: continuous limiter ensuring boundedness (beta limiter) enabled */

  // set numerical options of epsilon
  cs_equation_param_t *eqp = cs_field_get_equation_param(CS_F_(eps));
  eqp->isstpc = 2;

  int kccmin = cs_field_key_id("min_scalar");
  /* Set the Value for the Sup and Inf of the studied scalar
   * for the Gamma beta limiter for the temperature */
  cs_field_set_key_double(CS_F_(eps), kccmin, 0.);

  // set numerical options for k
  eqp = cs_field_get_equation_param(CS_F_(k));
  eqp->isstpc = 2;

  /* Set the Value for the Sup and Inf of the studied scalar
   * for the Gamma beta limiter for the temperature */
  cs_field_set_key_double(CS_F_(k), kccmin, 0.);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
