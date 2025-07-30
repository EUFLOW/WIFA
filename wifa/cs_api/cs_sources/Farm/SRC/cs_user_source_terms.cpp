/*============================================================================
 * Additional user-defined source terms for variable equations.
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

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_coupling.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"
#include "cs_wind_farm.h"
#include "prototypes.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Additional user-defined source terms for variable equations
 *        (momentum, scalars, turbulence...).
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
  const cs_real_t *cell_vol = domain->mesh_quantities->cell_vol;
  const cs_real_3_t *cell_cen = domain->mesh_quantities->cell_cen;

  const cs_mesh_t *m = domain->mesh;
  cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_lnum_t n_cells = m->n_cells;
  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;

  const cs_real_t *cpro_rom = CS_F_(rho)->val;
  const cs_real_3_t *vel = (const cs_real_3_t *)CS_F_(vel)->val;

  const cs_field_t *fld = cs_field_by_id(f_id);

  /*physical variables for density variation */
  cs_real_3_t *source_term = (cs_real_3_t *)cs_field_by_name("source_term")->val;
  cs_real_3_t *WF_source_term = (cs_real_3_t *)cs_field_by_name("WF_source_term")->val;
  cs_real_t *nudamp_top = cs_field_by_name("nudamp_top")->val;
  cs_real_t *nudamp_bound = cs_field_by_name("nudamp_bound")->val;

  cs_time_step_t *ts = cs_get_glob_time_step();

  char name[128];

  /* For velocity
   * ============ */

  if (fld == CS_F_(vel))
  {

    const cs_equation_param_t *eqp = cs_field_get_equation_param(CS_F_(vel));
    if (eqp->verbosity == 1)
    {
      bft_printf(" User source terms for variable %s \n ",
                 cs_field_get_label(CS_F_(vel)));
    }
    const cs_real_t dpi = atan(1.0) * 4.0;
    //
    cs_real_3_t *_st_exp = (cs_real_3_t *)st_exp;
    cs_real_33_t *_st_imp = (cs_real_33_t *)st_imp;
    //
    cs_real_t ugeo = 0.0;
    cs_real_t vgeo = 0.0;
    // get domain height
    sprintf(name, "Sommet");
    const cs_zone_t *zs = cs_boundary_zone_by_name_try(name);
    // z-coordinate of the center of gravity of the boundary zone
    cs_real_t zsommet = zs->cog[2];
    //
    //
    if (cs_glob_atmo_option->meteo_profile == 1)
    {
      int nbmett = cs_glob_atmo_option->nbmett; // nprofz
      int nbmetm = cs_glob_atmo_option->nbmetm; // nproft, dim_u_met, dim_pot_t_met, ..
      if (zsommet <= cs_glob_atmo_option->z_dyn_met[cs_glob_atmo_option->nbmett - 1])
      {
        ugeo = cs_intprf(nbmett,                         // nprofz
                         nbmetm,                         // nproft
                         cs_glob_atmo_option->z_dyn_met, // profz
                         cs_glob_atmo_option->time_met,  // proft
                         cs_glob_atmo_option->u_met,     // profu
                         zsommet,                        // xz
                         cs_glob_time_step->t_cur);      // t
        vgeo = cs_intprf(nbmett,                         // nprofz
                         nbmetm,                         // nproft
                         cs_glob_atmo_option->z_dyn_met, // profz
                         cs_glob_atmo_option->time_met,  // proft
                         cs_glob_atmo_option->v_met,     // profv
                         zsommet,                        // xz
                         cs_glob_time_step->t_cur);      // t
      }
      else
      {
        ugeo = cs_glob_atmo_option->u_met[nbmett - 1];
        vgeo = cs_glob_atmo_option->v_met[nbmett - 1];
      }
      if ((ts->nt_cur) <= 1)
      {
        bft_printf("Geostrophic wind interpolated from meteo_file at z= %.2f "
                   "is (u,v)=(%.2f,%.2f)\n ",
                   zsommet,
                   ugeo,
                   vgeo);
      }
    }
    else {
      const cs_real_3_t *cpro_met_vel = (cs_real_3_t *)(cs_field_by_name("meteo_velocity")->val);
      // geostrophic wind from prescribed meteo velocity
      cs_real_t closest_x, closest_y, closest_z;
      cs_lnum_t closest_id;
      int closest_id_rank;
      cs_real_t xyz_ref[3] = {0.0, 0.0, zsommet};
      cs_geom_closest_point(m->n_cells,
                            (const cs_real_3_t *)(mq->cell_cen),
                            xyz_ref,
                            &closest_id,
                            &closest_id_rank);

      if (closest_id_rank == cs_glob_rank_id) {
        ugeo = cpro_met_vel[closest_id][0];
        vgeo = cpro_met_vel[closest_id][1];
        closest_x = cell_cen[closest_id][0];
        closest_y = cell_cen[closest_id][1];
        closest_z = cell_cen[closest_id][2];
      }
      cs_parall_bcast(closest_id_rank, 1, CS_REAL_TYPE, &ugeo);
      cs_parall_bcast(closest_id_rank, 1, CS_REAL_TYPE, &vgeo);
      cs_parall_bcast(closest_id_rank, 1, CS_REAL_TYPE, &closest_x);
      cs_parall_bcast(closest_id_rank, 1, CS_REAL_TYPE, &closest_y);
      cs_parall_bcast(closest_id_rank, 1, CS_REAL_TYPE, &closest_z);
      if ((ts->nt_cur) <= 1) {
        bft_printf(
            "Geostrophic wind interpolated from prescribed meteo velocity at "
            "point (x,y,z)=(%.2f,%.2f,%.2f) is (u,v)=(%.2f,%.2f).\n ",
            closest_x,
            closest_y,
            closest_z,
            ugeo,
            vgeo);
      }
    }
    //
    const cs_real_t lat = cs_notebook_parameter_value_by_name("lat");
    cs_real_t fcorio;
    fcorio = 2.0 * 7.292115e-5 * sin(lat * dpi / 180.0);

    /* Coriolis force */
    if (cs_notebook_parameter_value_by_name("Coriolis") > 0)
    {
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
      {
        /* Geostrophic wind */
        _st_exp[c_id][0] = -cell_vol[c_id] * cpro_rom[c_id] * vgeo * fcorio;
        _st_imp[c_id][1][0] = -cell_vol[c_id] * cpro_rom[c_id] * fcorio;

        _st_exp[c_id][1] = cell_vol[c_id] * cpro_rom[c_id] * ugeo * fcorio;
        _st_imp[c_id][0][1] = cell_vol[c_id] * cpro_rom[c_id] * fcorio;
      }
    }
    /* Damping layer */
    if (cs_notebook_parameter_value_by_name("damping") > 0)
    {
      const cs_real_t gamma = cs_notebook_parameter_value_by_name("gamma");
      const cs_real_t nura = cs_notebook_parameter_value_by_name("nura");
      const cs_real_t lra = cs_notebook_parameter_value_by_name("Lra");
      const cs_real_t sra = cs_notebook_parameter_value_by_name("Sra");
      //
      const cs_real_t start_rad = cs_notebook_parameter_value_by_name("start_rad");
      //
      const cs_real_t zra = zsommet - lra;
      cs_real_t pref = cs_glob_atmo_constants->ps;
      cs_real_t rair = cs_glob_fluid_properties->r_pg_cnst;
      cs_real_t cp0 = cs_glob_fluid_properties->cp0;
      cs_real_t rscp = rair / cp0;
      cs_real_t psea = cs_glob_atmo_option->meteo_psea;
      cs_real_t theta0 = cs_glob_atmo_option->meteo_t0 * pow(pref / psea, rscp);
      int nbmett = cs_glob_atmo_option->nbmett; // nprofz
      int nbmetm = cs_glob_atmo_option->nbmetm; // nproft, dim_u_met, dim_pot_t_met, ..
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
      {
        // TODO : deal with case where domain is not centered
        cs_real_t x = cell_cen[c_id][0];
        cs_real_t y = cell_cen[c_id][1];
        cs_real_t z = cell_cen[c_id][2];
        cs_real_t radius = sqrt(pow(x, 2) + pow(y, 2));
        cs_real_t umet, vmet;
        if (radius > start_rad)
        {
          // TODO : deal with case where no meteo file
          umet = cs_intprf(nbmett,                         // nprofz
                           nbmetm,                         // nproft
                           cs_glob_atmo_option->z_dyn_met, // profz
                           cs_glob_atmo_option->time_met,  // proft
                           cs_glob_atmo_option->u_met,     // profu
                           z,                              // xz
                           cs_glob_time_step->t_cur);      // t
          vmet = cs_intprf(nbmett,                         // nprofz
                           nbmetm,                         // nproft
                           cs_glob_atmo_option->z_dyn_met, // profz
                           cs_glob_atmo_option->time_met,  // proft
                           cs_glob_atmo_option->v_met,     // profv
                           z,                              // xz
                           cs_glob_time_step->t_cur);      // t
          nudamp_bound[c_id] = nura * sqrt(9.81 * gamma / theta0) * (1 - cos(dpi / sra * (radius - start_rad) / lra));
          _st_exp[c_id][0] += -cell_vol[c_id] * cpro_rom[c_id] * nudamp_bound[c_id] * (vel[c_id][0] - umet);
          _st_exp[c_id][1] += -cell_vol[c_id] * cpro_rom[c_id] * nudamp_bound[c_id] * (vel[c_id][1] - vmet);
          _st_exp[c_id][2] += -cell_vol[c_id] * cpro_rom[c_id] * nudamp_bound[c_id] * (vel[c_id][2]);
        }
        if (z >= zra)
        {
          nudamp_top[c_id] = nura * sqrt(9.81 * gamma / theta0) * (1 - cos(dpi / sra * (z - zra) / lra));
          _st_exp[c_id][0] += -cell_vol[c_id] * cpro_rom[c_id] * nudamp_top[c_id] * (vel[c_id][0] - ugeo);
          _st_exp[c_id][1] += -cell_vol[c_id] * cpro_rom[c_id] * nudamp_top[c_id] * (vel[c_id][1] - vgeo);
          _st_exp[c_id][2] += -cell_vol[c_id] * cpro_rom[c_id] * nudamp_top[c_id] * (vel[c_id][2]);
        }
      }
    }

    //
    // Compute source term for turbines
    //
    cs_wind_farm_compute_source_term();

    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
    {
      _st_exp[c_id][0] += WF_source_term[c_id][0];
      _st_exp[c_id][1] += WF_source_term[c_id][1];
      _st_exp[c_id][2] += WF_source_term[c_id][2];
      source_term[c_id][0] = _st_exp[c_id][0];
      source_term[c_id][1] = _st_exp[c_id][1];
      source_term[c_id][2] = _st_exp[c_id][2];
    }
  }

  /* For the temperature forcing
   * =========================== */

  if (fld == CS_F_(t))
  {
    if (cs_glob_atmo_option->meteo_dlmo > 0)
    {
      /* u* */
      const cs_real_t ustar = cs_notebook_parameter_value_by_name("ustar");
      /* theta* */
      const cs_real_t tstar = cs_notebook_parameter_value_by_name("tstar");
      const cs_real_t cp0 = cs_glob_fluid_properties->cp0;
      const cs_real_t zi = cs_notebook_parameter_value_by_name("zi");
      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
      {
        cs_real_t z = cell_cen[c_id][2];
        cs_real_t sterm = cp0 * cpro_rom[c_id] * ustar * tstar * 2.0 / zi * pow(1 - z / zi, 2.0 - 1);
        if (z < zi)
        {
          st_exp[c_id] = cell_vol[c_id] * sterm;
        }
      }
    }
  }

  /****************************************/
  /* Dyunkerke model */
  if (    f_id == CS_F_(eps)->id
      &&  cs_notebook_parameter_value_by_name("Dyunkerke") > 0) {

    const cs_lnum_t n_b_faces = m->n_b_faces;
    const cs_lnum_t n_i_faces = m->n_i_faces;
    cs_real_t *cromo = CS_F_(rho)->val;
    cs_real_t *cpro_pcvto = CS_F_(mu_t)->val;
    cs_real_t sigmak = 1.0;
    const cs_equation_param_t *eqp_k = cs_field_get_equation_param_const(CS_F_(k));
    const cs_field_bc_coeffs_t *bc_coeffs;
    cs_real_t *vol_divmugradk = nullptr;
    BFT_MALLOC(vol_divmugradk, n_cells_ext, cs_real_t);
    cs_real_t *w3 = nullptr;
    BFT_MALLOC(w3, n_cells_ext, cs_real_t);
    cs_real_t *viscf, *viscb;
    BFT_MALLOC(viscf, n_i_faces, cs_real_t);
    BFT_MALLOC(viscb, n_b_faces, cs_real_t);
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
      w3[c_id] = cpro_pcvto[c_id] / cromo[c_id] / sigmak;

    cs_face_viscosity(m, mq, eqp_k->imvisf, w3, viscf, viscb);

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
    for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++)
    {
      /* Already contains cell volume */
      st_exp[c_id] = CS_F_(eps)->val_pre[c_id] / CS_F_(k)->val_pre[c_id] * cs_turb_ce1 * cromo[c_id] * fmax(0., -vol_divmugradk[c_id]);
      if (f_eps_transport != nullptr)
      {
        f_eps_transport->val[c_id] = st_exp[c_id] / cell_vol[c_id] / cromo[c_id];
      }
    }
    BFT_FREE(w3);
    BFT_FREE(vol_divmugradk);
    BFT_FREE(viscf);
    BFT_FREE(viscb);
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
