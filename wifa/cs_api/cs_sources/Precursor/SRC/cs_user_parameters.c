/*============================================================================
 * User functions for input of calculation parameters.
 *============================================================================*/

/* VERS */

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2020 EDF S.A.

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

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_parameters-base.c
 *
 * \brief User functions for input of calculation parameters.
 *
 * See \ref parameters for examples.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

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
cs_user_model(void)
{
  {
    cs_wall_functions_t *wf = cs_get_glob_wall_functions();
     wf->iwalfs = CS_WALL_F_S_MONIN_OBUKHOV;
     wf->iwallf = CS_WALL_F_2SCALES_SMOOTH_ROUGH;
  }
  /* Atmospheric module options
   */
  cs_glob_physical_model_flag[CS_ATMOSPHERIC] = CS_ATMO_DRY;

  /* Option to compute ground elevation in the domain */
  cs_glob_atmo_option->compute_z_ground = true;

  cs_glob_atmo_option->meteo_z0 = cs_notebook_parameter_value_by_name("z0");

  cs_glob_atmo_option->meteo_t0 = cs_notebook_parameter_value_by_name("t0");

  cs_glob_atmo_option->meteo_dlmo = cs_notebook_parameter_value_by_name("Lmoinv");

  cs_glob_atmo_option->meteo_uref = cs_notebook_parameter_value_by_name("ureff");

  /* Automatic open boundary conditions
   *   1: meteo mass flow rate is imposed with a constant large scale
   *      pressure gradient
   *   2: same plus velocity profile imposed at ingoing faces
   */
  //  cs_glob_atmo_option->open_bcs_treatment = 1;

  /* Geographic position
   *  longitude: longitude of the domain origin
   *  latitude: latitude of the domain origin
   */

  //Needed for to compute the Coriolis force
  //cs_glob_atmo_option->longitude = cs_notebook_parameter_value_by_name("long");
  //cs_glob_atmo_option->latitude = cs_notebook_parameter_value_by_name("lat");

  cs_parameters_add_property("tke_transport",
                             1,CS_MESH_LOCATION_CELLS);
  cs_parameters_add_property("eps_transport",
                             1,CS_MESH_LOCATION_CELLS);

  /* To post-process u* and uk */
  cs_parameters_add_property("boundary_ustar",
                             1,
                             CS_MESH_LOCATION_CELLS);

  cs_parameters_add_property("boundary_uk",
                             1,
                             CS_MESH_LOCATION_CELLS);

}


void
cs_user_parameters(cs_domain_t *domain)
{

  /*! [param_var_limiter_choice] */
  {

    /* ischcv is the type of convective scheme:
       0: second order linear upwind
       1: centered
       2: pure upwind gradient in SOLU
       3: blending SOLU and centered
       4: NVD/TVD Scheme */

    /* isstpc:
      0: slope test enabled
      1: slope test disabled (default)
      2: continuous limiter ensuring boundedness (beta limiter) enabled */

    cs_fluid_properties_t *phys_pro = cs_get_glob_fluid_properties();
    cs_velocity_pressure_param_t *vp_param = cs_get_glob_velocity_pressure_param();
    vp_param->igpust=0;

    cs_time_step_t *ts = cs_get_glob_time_step();
    ts->nt_max = cs_notebook_parameter_value_by_name("precntmax");

    /* Warning, meteo file does not overwrite reference values... */
    cs_real_t rair = phys_pro->r_pg_cnst;
    /* Reference fluid properties set from meteo values */
    phys_pro->p0 = cs_glob_atmo_option->meteo_psea;

    //phys_pro->t0 = cs_glob_atmo_option->meteo_t0; /* ref temp T0 */
    //phys_pro->ro0 = phys_pro->p0/(rair * cs_glob_atmo_option->meteo_t0); /* ref density T0 */

    cs_var_cal_opt_t vcopt;
    int key_cal_opt_id = cs_field_key_id("var_cal_opt");

    cs_field_get_key_struct(CS_F_(eps), key_cal_opt_id, &vcopt);

    vcopt.ischcv = 1;
    vcopt.isstpc = 2;
    cs_field_set_key_struct(CS_F_(eps), key_cal_opt_id, &vcopt);
    int kccmin = cs_field_key_id("min_scalar");

    /* Set the Value for the Sup and Inf of the studied scalar
     * for the Gamma beta limiter for the temperature */
    cs_field_set_key_double(CS_F_(eps), kccmin,0.);

    cs_field_get_key_struct(CS_F_(k), key_cal_opt_id, &vcopt);
    vcopt.ischcv = 1;
    vcopt.isstpc = 2;
    cs_field_set_key_struct(CS_F_(k), key_cal_opt_id, &vcopt);

    /* Set the Value for the Sup and Inf of the studied scalar
     * for the Gamma beta limiter for the temperature */
    cs_field_set_key_double(CS_F_(k), kccmin, 0.);

  }

}
/*----------------------------------------------------------------------------*/

END_C_DECLS
