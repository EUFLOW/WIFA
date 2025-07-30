/*============================================================================
 * Definition of the calculation mesh.
 *
 * Mesh-related user functions (called in this order):
 *   1) Manage the exchange of data between code_saturne and the pre-processor
 *   2) Define (conforming or non-conforming) mesh joinings.
 *   3) Define (conforming or non-conforming) periodicity.
 *   4) Define thin walls.
 *   5) Modify the geometry and mesh.
 *   6) Smoothe the mesh.
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
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"
#include "prototypes.h"
/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_mesh.cpp
 *
 * \brief Definition and modification of the calculation mesh.
 *
 * See \ref cs_user_mesh for examples.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Modify geometry and mesh.
 *
 * \param[in,out] mesh  pointer to a cs_mesh_t structure
 */
/*----------------------------------------------------------------------------*/

void yawing(cs_real_t *coords,
            cs_real_t wind_dir,
            cs_real_t *mean_coords,
            cs_real_t *yawed_coords)
{
  const cs_real_t dpi = atan(1.0) * 4.0;
  cs_real_t ra_wind_dir = (270.0 - wind_dir) * dpi / 180.0; // radian absolute wind_dir

  yawed_coords[0] = cos(ra_wind_dir) * (coords[0] - mean_coords[0])
                  - sin(ra_wind_dir) * (coords[1] - mean_coords[1])
                  + mean_coords[0];

  yawed_coords[1] = sin(ra_wind_dir) * (coords[0] - mean_coords[0])
                  + cos(ra_wind_dir) * (coords[1] - mean_coords[1])
                  + mean_coords[1];

  yawed_coords[2] = coords[2];
}

void cs_user_mesh_modify(cs_mesh_t *mesh) {
  /* Refine a selected portion of a mesh */

  /* [mesh_modify_refine_1] */

  /******************************************
   * READ Wind Turbines from coordinates file
   *******************************************/
  cs_real_t wind_dir = cs_glob_atmo_option->meteo_angle;

  // file name and tables
  char WT_file_name[1024];
  sprintf(WT_file_name, "turbines_info.csv");
  cs_real_t *WT_x_coords = nullptr;
  cs_real_t *WT_y_coords = nullptr;
  cs_real_t *WT_z_coords = nullptr;
  cs_real_t *WT_diameters = nullptr;
  cs_lnum_t *WT_types = nullptr;

  // calculate number of lines in the file
  FILE *stream = fopen(WT_file_name, "r");
  char line[1024];
  char *tok0;
  const char sep0[4] = ";";
  size_t n_WT = 0;
  size_t i = -1;
  // read header
  // fgets(line, 1024, stream);
  while (fgets(line, 1024, stream) != nullptr)
  {
    i++;
    if (i > 0)
    {
      // try to parse line
      tok0 = strtok(line, sep0);
      if (tok0 != nullptr)
      {
        n_WT++;
      }
    }
  }
  fclose(stream);

  // allocate memory for all members of emission structure
  BFT_MALLOC(WT_x_coords, n_WT, cs_real_t);
  BFT_MALLOC(WT_y_coords, n_WT, cs_real_t);
  BFT_MALLOC(WT_z_coords, n_WT, cs_real_t);
  BFT_MALLOC(WT_diameters, n_WT, cs_real_t);
  BFT_MALLOC(WT_types, n_WT, cs_lnum_t);

  // read and parse this file
  stream = fopen(WT_file_name, "r");
  char *tok;
  const char sep[4] = ",";
  size_t line_count = -1;
  i = -1;

  while (fgets(line, 1024, stream) != nullptr)
  {
    i++;
    if (++line_count >= n_WT + 1)
      break;

    // skip header
    if (i > 0)
    {
      // parse line
      tok = strtok(line, sep);
      if (tok == nullptr)
        continue;
      WT_x_coords[i - 1] = atof(tok);

      tok = strtok(nullptr, sep);
      if (tok == nullptr)
        continue;
      WT_y_coords[i - 1] = atof(tok);

      tok = strtok(nullptr, sep);
      if (tok == nullptr)
        continue;
      WT_z_coords[i - 1] = atof(tok);

      tok = strtok(nullptr, sep);
      if (tok == nullptr)
        continue;
      WT_diameters[i - 1] = atof(tok);

      tok = strtok(nullptr, sep);
      if (tok == nullptr)
        continue;
      WT_types[i - 1] = atoi(tok);
    }
  }
  fclose(stream);
  //
  /******************************************
   * MARK TURBINE LOCATIONS USING CYLINDERS
   *******************************************/
  cs_real_t mean_coords[3];
  char criteria[100];
  //
  cs_lnum_t n_selected_cells = 0;
  cs_lnum_t *selected_cells = nullptr;
  //
  cs_real_t disk_1_coords[3];
  cs_real_t disk_2_coords[3];
  cs_real_t yawed_disk_1_coords[3];
  cs_real_t yawed_disk_2_coords[3];
  //
  cs_real_t base_corner_coords[3];
  cs_real_t xTranslate_base_corner_coords[3];
  cs_real_t yTranslate_base_corner_coords[3];
  cs_real_t zTranslate_base_corner_coords[3];
  //
  char name[128]; // for Actuator Disk (AD) zone name
  //
  cs_real_t AD_mesh_cell_size = cs_notebook_parameter_value_by_name("AD_mesh_cell_size");
  cs_real_t AD_half_rotor_thickness = 1.2 * AD_mesh_cell_size;
  cs_lnum_t start_WT_count, end_WT_count;
  if (cs_notebook_parameter_value_by_name("isol") > 0)
  {
    start_WT_count = cs_notebook_parameter_value_by_name("isol") - 1;
    end_WT_count = cs_notebook_parameter_value_by_name("isol");
  }
  else
  {
    start_WT_count = 0;
    end_WT_count = n_WT;
  }
  for (cs_lnum_t WT_count = start_WT_count; WT_count < end_WT_count;
       WT_count++)
  {
    //
    cs_real_t WT_d = WT_diameters[WT_count];
    cs_real_t WT_radius = WT_d / 2.;
    /***********************************************************/
    /* Cylinder corresponding to AD zone */
    if (cs_notebook_parameter_value_by_name("control") > 0)
    {
      base_corner_coords[0] = WT_x_coords[WT_count] - 1.1 * WT_radius;
      base_corner_coords[1] = WT_y_coords[WT_count] - 1.1 * WT_radius;
      base_corner_coords[2] = WT_z_coords[WT_count] - 1.1 * WT_radius;
      // dx
      xTranslate_base_corner_coords[0] = 2.2 * WT_radius;
      xTranslate_base_corner_coords[1] = 0.0;
      xTranslate_base_corner_coords[2] = 0.0;
      // dy
      yTranslate_base_corner_coords[0] = 0.0;
      yTranslate_base_corner_coords[1] = 2.2 * WT_radius;
      yTranslate_base_corner_coords[2] = 0.0;
      // dz
      zTranslate_base_corner_coords[0] = 0.0;
      zTranslate_base_corner_coords[1] = 0.0;
      zTranslate_base_corner_coords[2] = 2.2 * WT_radius;

      // select zone
      sprintf(criteria,
              "box[%.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, "
              "%.2f, %.2f]",
              base_corner_coords[0],
              base_corner_coords[1],
              base_corner_coords[2],
              xTranslate_base_corner_coords[0],
              xTranslate_base_corner_coords[1],
              xTranslate_base_corner_coords[2],
              yTranslate_base_corner_coords[0],
              yTranslate_base_corner_coords[1],
              yTranslate_base_corner_coords[2],
              zTranslate_base_corner_coords[0],
              zTranslate_base_corner_coords[1],
              zTranslate_base_corner_coords[2]);
    }
    else
    {
      disk_1_coords[0] = WT_x_coords[WT_count] - 2.0 * AD_half_rotor_thickness;
      disk_1_coords[1] = WT_y_coords[WT_count];
      disk_1_coords[2] = WT_z_coords[WT_count];
      //
      disk_2_coords[0] = WT_x_coords[WT_count] + 2.0 * AD_half_rotor_thickness;
      disk_2_coords[1] = WT_y_coords[WT_count];
      disk_2_coords[2] = WT_z_coords[WT_count];

      // yaw disks centers in the cylinder local referential
      mean_coords[0] = 0.5 * (disk_1_coords[0] + disk_2_coords[0]);
      mean_coords[1] = 0.5 * (disk_1_coords[1] + disk_2_coords[1]);
      mean_coords[2] = 0.5 * (disk_1_coords[2] + disk_2_coords[2]);
      yawing(disk_1_coords, wind_dir, mean_coords, yawed_disk_1_coords);
      yawing(disk_2_coords, wind_dir, mean_coords, yawed_disk_2_coords);

      // select zone
      sprintf(criteria,
              "cylinder[%.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f]",
              yawed_disk_1_coords[0],
              yawed_disk_1_coords[1],
              yawed_disk_1_coords[2],
              yawed_disk_2_coords[0],
              yawed_disk_2_coords[1],
              yawed_disk_2_coords[2],
              1.5 * WT_radius);
    }

    bft_printf("Turbine_%d selection criteria : %s\n", WT_count + 1, criteria);

    BFT_MALLOC(selected_cells, mesh->n_cells, cs_lnum_t);

    cs_selector_get_cell_list(criteria, &n_selected_cells, selected_cells);

    // tag AD zone as volume zone
    sprintf(name, "turbine_%d", WT_count + 1);
    cs_volume_zone_define(
        name,
        criteria,
        CS_VOLUME_ZONE_SOURCE_TERM); /* Cylinder corresponding to AD zone */
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
