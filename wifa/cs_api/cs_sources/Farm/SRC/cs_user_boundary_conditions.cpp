/*============================================================================
 * User definition of boundary conditions.
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
#include <stdio.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief User definition of boundary conditions
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 * \param[in, out]  bc_type  boundary face types
 *
 * The icodcl and rcodcl arrays are pre-initialized based on default
 * and GUI-defined definitions, and may be modified here.
 *
 * For a given variable field f, and a given face "face_id", these arrays
 * may be used as follows:
 *
 * - Boundary condition type code given at:
 *   f->bc_coeffs->icodcl[face_id]
 *
 * - Dirichlet value defined at:
 *   f->bc_coeffs->rcodcl1[face_id]
 *
 * - Interior exchange coefficient (infinite if no exchange) at:
 *   f->bc_coeffs->rcodcl2[face_id]
 *
 * - Flux density defined at:
 *   f->bc_coeffs->rcodcl3[face_id]
 *
 * For vector or tensor fields, these arrays are not interleaved,
 * so for a given face "face_id" and field component "comp_id", access
 * is as follows (where n_b_faces is domain->mesh->n_b_faces):
 *
 *   f->bc_coeffs->rcodcl1[n_b_faces*comp_id + face_id]\n
 *   f->bc_coeffs->rcodcl2[n_b_faces*comp_id + face_id]\n
 *   f->bc_coeffs->rcodcl3[n_b_faces*comp_id + face_id]\n\n
 *
 * Only the icodcl code values from the first component are used in the case
 * of vector or tensor fields, so the icodcl values can be defined as for
 * a scalar.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_boundary_conditions(cs_domain_t  *domain,
                            int           bc_type[])
{
  CS_UNUSED(domain);
  CS_UNUSED(bc_type);

  cs_field_t *f_roughness = cs_field_by_name("boundary_roughness");
  cs_field_t *f_thermal_roughness = cs_field_by_name("boundary_thermal_roughness");

  cs_real_t dlmo = cs_glob_atmo_option->meteo_dlmo;
  const cs_zone_t *z = cs_boundary_zone_by_name("Sol");

  for (cs_lnum_t face_count = 0; face_count < z->n_elts; face_count++) {
    const cs_lnum_t face_id = z->elt_ids[face_count];
    f_roughness->val[face_id] = cs_glob_atmo_option->meteo_z0;
    /* /\* How to treat thermal rugosity is still uncertain *\/ */
    f_thermal_roughness->val[face_id] = cs_glob_atmo_option->meteo_z0;
    if (dlmo > 0) {
      CS_F_(t)->bc_coeffs->icodcl[face_id] = 6;
      CS_F_(t)->bc_coeffs->rcodcl1[face_id] = cs_glob_atmo_option->meteo_t0;
    }
  }

}

/*----------------------------------------------------------------------------*/

END_C_DECLS
