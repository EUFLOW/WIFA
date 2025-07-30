/*============================================================================
 * General-purpose user-defined functions called before time stepping, at
 * the end of each time step, and after time-stepping.
 *
 * These can be used for operations which do not fit naturally in any other
 * dedicated user function.
 *============================================================================*/

/* VERS */

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2023 EDF S.A.

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

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_extra_operations.c
 *
 * \brief This function is called at the end of each time step, and has a very
 * general purpose (i.e. anything that does not have another dedicated
 * user function)
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/


/*----------------------------------------------------------------------------*/
/*!
 * \brief This function is called at the end of each time step.
 *
 * It has a very general purpose, although it is recommended to handle
 * mainly postprocessing or data-extraction type operations.
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_extra_operations
void
cs_user_extra_operations(cs_domain_t     *domain)
{
//  if (cs_glob_atmo_option->meteo_dlmo>0){
  const cs_real_t  cs_energy_solver=cs_notebook_parameter_value_by_name("energy");
  /* mesh quantities */
  const cs_mesh_t *m = domain->mesh;
  const cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_lnum_t  n_cells = m->n_cells;
  const cs_real_3_t *cell_cen = (const cs_real_3_t *)mq->cell_cen;

  cs_field_t *f_thm = cs_thermal_model_field();

  int nbmett = cs_glob_atmo_option->nbmett; //nprofz
  int nbmetm = cs_glob_atmo_option->nbmetm; //nproft, dim_u_met, dim_pot_t_met, ..

  if (cs_energy_solver==0){
    // 1: continue solving for the temperature profile
    // 0: "activate" this function, which keeps the temperature profile constant
    // Start with energy=0, to fix the temp at the start of the calculation
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      const cs_real_t z = cell_cen[c_id][2];
  	  f_thm->val[c_id] = cs_intprf(nbmett, //nprofz
                                  nbmetm, //nproft
                                  cs_glob_atmo_option->z_dyn_met, //profz
                                  cs_glob_atmo_option->time_met, //proft
                                  cs_glob_atmo_option->pot_t_met, //prof theta
                                  z, //xz
                                  cs_glob_time_step->t_cur); //t ;
    }
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
