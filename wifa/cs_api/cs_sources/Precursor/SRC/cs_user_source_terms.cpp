/*============================================================================
 * Additional user-defined source terms for variable equations.
 *============================================================================*/

/* VERS */

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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
 * \brief Function called at each time step to define physical properties.
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 * \param[in]       f_id     field id of the variable
 * \param[out]      st_exp   explicit source term
 * \param[out]      st_imp   implicit part of the source term
 */
/*----------------------------------------------------------------------------*/

void
cs_user_source_terms(cs_domain_t  *domain,
                     int           f_id,
                     cs_real_t    *st_exp,
                     cs_real_t    *st_imp)
{
  const cs_mesh_t *m = domain->mesh;
  cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;

  const cs_real_t *cell_vol = mq->cell_vol;
  const cs_real_3_t *cell_cen
    = (const cs_real_3_t *)mq->cell_cen;

  const cs_real_t  cp0 = cs_glob_fluid_properties->cp0;

  const cs_real_t *cpro_rom = CS_F_(rho)->val;

  const cs_real_3_t *vel = (const cs_real_3_t *)CS_F_(vel)->val;

  const cs_field_t *fld = cs_field_by_id(f_id);

  const cs_real_t  u_ref = cs_notebook_parameter_value_by_name("ureff");
  const cs_real_t  z_ref = cs_notebook_parameter_value_by_name("zref");

  const cs_real_t  zi = cs_notebook_parameter_value_by_name("zi");

  /* u* */
  const cs_real_t ustar = cs_notebook_parameter_value_by_name("ustar");

  /* theta* */
  const cs_real_t tstar = cs_notebook_parameter_value_by_name("tstar");

  cs_time_step_t *ts = cs_get_glob_time_step();
  const cs_real_t  cs_energy_solver=cs_notebook_parameter_value_by_name("energy");

  /* For velocity
   * ============ */

  if (fld == CS_F_(vel)) {

    const cs_equation_param_t *eqp = cs_field_get_equation_param(CS_F_(vel));
    if (eqp->verbosity == 1) {
      bft_printf(" User source terms for variable %s \n ",
                 cs_field_get_label(CS_F_(vel)));
    }

    /* relaxation with the Wangara geostrophic wind */
    cs_real_3_t *_st_exp = (cs_real_3_t *)st_exp;
    cs_real_33_t *_st_imp = (cs_real_33_t *)st_imp;

    /* 63°N ! /!\ if we change the value of fcorio, must also change
     * fcorio in cs_user_postprocess.c */
    const cs_real_t dpi = atan(1.0)*4.0;
    const cs_real_t  lat = cs_notebook_parameter_value_by_name("lat");
    cs_real_t fcorio;
    fcorio = 2.0*7.292115e-5*sin(lat*dpi/180.0);

    if ((ts->nt_cur) > 20000 && (ts->nt_cur)%1000==0 && cs_energy_solver==0){
      cs_real_t uzrefn = u_ref;
      cs_real_t uzref, vzref;
      cs_lnum_t closest_id;
      int closest_id_rank;
      cs_real_t xyz_ref[3] = {0.0, 0.0, z_ref};
      cs_geom_closest_point(m->n_cells,
  		      (const cs_real_3_t *)(mq->cell_cen),
  		      xyz_ref,
  		      &closest_id,
  		      &closest_id_rank);

      if (closest_id_rank == cs_glob_rank_id) {
  	    uzref = vel[closest_id][0];
  	    vzref = vel[closest_id][1];
      }
      cs_parall_bcast(closest_id_rank, 1, CS_REAL_TYPE, &uzref);
      cs_parall_bcast(closest_id_rank, 1, CS_REAL_TYPE, &vzref);

      uzrefn = sqrt(pow(uzref, 2.0)+pow(vzref, 2.0));
      if (fabs(u_ref-uzrefn)>0.05){
          cs_glob_atmo_option->meteo_uref += CS_MAX(CS_MIN(u_ref-uzrefn, 0.01), -0.01);
      }
      else{
        ts->nt_max = ts->nt_cur + 7000;
        //modification du param notebook pour passer a energy=1
        cs_notebook_parameter_set_value("energy", 1);
      }
      bft_printf("uzrefn = %.2f ; uref = %.2f\n ",uzrefn, cs_glob_atmo_option->meteo_uref);
    }

    if(cs_glob_atmo_option->meteo_dlmo>0 && ts->nt_cur==1)
    {
      cs_real_t utop;
      cs_lnum_t closest_id_top;
      int closest_id_rank_top;
      cs_real_t xyz_ref_top[3] = {0.0, 0.0, 20000};
      cs_geom_closest_point(m->n_cells,
  		      (const cs_real_3_t *)(mq->cell_cen),
  		      xyz_ref_top,
  		      &closest_id_top,
  		      &closest_id_rank_top);

      if (closest_id_rank_top == cs_glob_rank_id) {
  	    utop = vel[closest_id_top][0];
      }
      cs_parall_bcast(closest_id_rank_top, 1, CS_REAL_TYPE, &utop);
	  cs_glob_atmo_option->meteo_uref = utop;
	}

    /* Geostrophic wind */
    cs_real_t xuref = cs_glob_atmo_option->meteo_uref;
    cs_real_t xvref = 0;

    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      _st_exp[c_id][0]    = -cell_vol[c_id]*cpro_rom[c_id]*xvref*fcorio;
      _st_imp[c_id][1][0] = -cell_vol[c_id]*cpro_rom[c_id]*fcorio;

      _st_exp[c_id][1]    = cell_vol[c_id]*cpro_rom[c_id]*xuref*fcorio;
      _st_imp[c_id][0][1] = cell_vol[c_id]*cpro_rom[c_id]*fcorio;
    }
  }

  /* For the temperature forcing
   * =========================== */

  if (fld == CS_F_(t)) {
	if (cs_glob_atmo_option->meteo_dlmo>0)
	{
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        cs_real_t z = cell_cen[c_id][2];
        cs_real_t sterm = cp0*cpro_rom[c_id]*ustar*tstar*2.0/zi*pow(1-z/zi, 2.0-1);
        if (z < zi) {
          st_exp[c_id] = cell_vol[c_id]*sterm;
        }
      }
    }
  }

  /* Dyunkerke model */
  if (f_id == CS_F_(eps)->id &&
      cs_notebook_parameter_value_by_name("Dyunkerke") > 0) {
    const cs_lnum_t n_b_faces = m->n_b_faces;
    const cs_lnum_t n_i_faces = m->n_i_faces;
    cs_real_t *cromo = CS_F_(rho)->val;
    cs_real_t *cpro_pcvto = CS_F_(mu_t)->val;
    cs_real_t sigmak=1.0;
    const cs_equation_param_t *eqp_k
      = cs_field_get_equation_param_const(CS_F_(k));
    cs_real_t hint;
    const cs_field_bc_coeffs_t *bc_coeffs;
    cs_real_t *vol_divmugradk = nullptr;
    BFT_MALLOC(vol_divmugradk, n_cells_ext, cs_real_t);
    cs_real_t *w3 = nullptr;
    BFT_MALLOC(w3, n_cells_ext, cs_real_t);
    cs_real_t *viscf, *viscb;
    BFT_MALLOC(viscf, n_i_faces, cs_real_t);
    BFT_MALLOC(viscb, n_b_faces, cs_real_t);
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
      w3[c_id] = cpro_pcvto[c_id] / cromo[c_id] / sigmak;
    }

    cs_face_viscosity(m,
                      mq,
                      eqp_k->imvisf,
                      w3,
                      viscf,
                      viscb);

    bc_coeffs = CS_F_(k)->bc_coeffs;

    /* Compute - div(mu_T/sigmak grad (k)) time the volume of the cell */
    cs_diffusion_potential(CS_F_(k)->id,
                           m,
                           mq,
                           1, /* init */
                           0, /* inc */
                           eqp_k->imrgra,
                           eqp_k->nswrgr,
                           eqp_k->imligr,
                           0, /* iphydp */
                           eqp_k->iwgrec,
                           eqp_k->verbosity,
                           eqp_k->epsrgr,
                           eqp_k->climgr,
                           nullptr,
                           CS_F_(k)->val_pre,
                           bc_coeffs,
                           viscf,
                           viscb,
                           w3,
                           vol_divmugradk);


    cs_field_t *f_eps_transport = cs_field_by_name_try("eps_transport");
    for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id ++) {
      /* Already contains cell volume */
      st_exp[cell_id] = CS_F_(eps)->val_pre[cell_id] / CS_F_(k)->val_pre[cell_id] * cs_turb_ce1 * cromo[cell_id] * fmax( 0., -vol_divmugradk[cell_id] ) ;
      if (f_eps_transport != nullptr) {
        f_eps_transport->val[cell_id] = st_exp[cell_id]/cell_vol[cell_id]/cromo[cell_id];
      }
    }

    BFT_FREE(w3);
    BFT_FREE(vol_divmugradk);

  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
