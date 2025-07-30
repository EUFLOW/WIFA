#ifndef __PROTOTYPES__H__
#define __PROTOTYPES__H__

/*----------------------------------------------------------------------------*/

/* VERS */

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Global variables definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/

void yawing(cs_real_t *coords,
            cs_real_t wind_dir,
            cs_real_t *mean_coords,
            cs_real_t *yawed_coords);
void yawing2(cs_real_t *coords,
             cs_real_t wind_dir,
             cs_real_t *mean_coords,
             cs_real_t *yawed_coords);

cs_real_t interp_ct_or_cp(cs_real_t *ct_or_cp_values,
                          cs_real_t *ct_or_cp_speeds,
                          cs_lnum_t number_of_values,
                          cs_real_t disk_velocity);
/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __PROTOTYPES__H__ */
