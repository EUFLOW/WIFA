#include "cs_headers.h"
#include "cs_defs.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "cs_wind_farm.h"
#include "bft_mem_usage.h"
#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "fvm_selector.h"

#include "cs_halo.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"

#include "cs_selector.h"


//declaration of resulted variable
static cs_wind_farm _wind_farm = {
  .WT_file_name = NULL,
  .n_WT = 0,
  .n_WT_types = 0,
  .total_power =0.0,
  .WT_coords = NULL,
  .WT_diameters = NULL,
  .WT_types = NULL,
  .WT_volume = NULL,
  .WT_uvect = NULL,
  .WT_u = NULL,
  .WT_rho = NULL,
  .WT_hub_u = NULL,
  .WT_dir = NULL,
  .WT_ra_dir = NULL,
  .WT_ctstar = NULL,
  .WT_cpstar = NULL,
  .WT_thrust = NULL,
  .WT_power = NULL
};

/*============================================================================
 *  Global static variables
 *============================================================================*/

struct _cs_ctstar_table {
  cs_real_t *ctstar_speeds;
  cs_real_t *ctstar_values;
  cs_lnum_t n_ctstar;
};

struct _cs_cpstar_table {
  cs_real_t *cpstar_speeds;
  cs_real_t *cpstar_values;
  cs_lnum_t n_cpstar;
};

static cs_ctstar_table **_cs_ctstar_table_array = NULL;
static cs_cpstar_table **_cs_cpstar_table_array = NULL;

/*============================================================================
 * Private function definitions
 *============================================================================*/
/*----------------------------------------------------------------------------*/

static void
_coeff_tables_create(int n_WT_types)
{
  /******************************************
   * READ Cp and Ct from data files
   *******************************************/
  //for IEC 61400-12 air density correction
  cs_real_t standard_rho = 1.225;
  //
  char cpstar_file_names[n_WT_types][1024];
  char ctstar_file_names[n_WT_types][1024];

  BFT_MALLOC(_cs_ctstar_table_array, n_WT_types, cs_ctstar_table*);
  BFT_MALLOC(_cs_cpstar_table_array, n_WT_types, cs_cpstar_table*);

  for (cs_lnum_t WT_type_count=0; WT_type_count < n_WT_types; WT_type_count ++){
  /*  cs_ctstar_table ctstar_table0 = {
      .ctstar_speeds = NULL,
      .ctstar_values = NULL,
      .n_ctstar = 0
    };
    cs_cpstar_table cpstar_table0 = {
      .cpstar_speeds = NULL,
      .cpstar_values = NULL,
      .n_cpstar = 0
    };

    cs_ctstar_table * ctstar_table = &ctstar_table0;
    cs_cpstar_table * cpstar_table = &cpstar_table0;
*/
    cs_ctstar_table * ctstar_table = NULL;
    cs_cpstar_table * cpstar_table = NULL;
    BFT_MALLOC(ctstar_table, 1, cs_ctstar_table);
    BFT_MALLOC(cpstar_table, 1, cs_cpstar_table);
    sprintf(ctstar_file_names[WT_type_count],"ctstar_table_type%d.csv",WT_type_count+1);
    sprintf(cpstar_file_names[WT_type_count],"cpstar_table_type%d.csv",WT_type_count+1);

    //calculate number of lines in the file
    FILE* stream_ct = fopen(ctstar_file_names[WT_type_count], "r");
    char *tok2;
    const char sep2[4] = ",";
    char line[1024];
    cs_lnum_t local_n_ct = 0;
    //read header
    //fgets(line, 1024, stream);
    while (fgets(line, 1024, stream_ct) != NULL){
      // try to parse line
      tok2 = strtok(line, sep2);
      if (tok2 != NULL)
  	    local_n_ct++;
    }
    ctstar_table->n_ctstar = local_n_ct;
    fclose(stream_ct);

    //allocate memory for all members of structure
    BFT_MALLOC(ctstar_table->ctstar_speeds, local_n_ct, cs_real_t);
    BFT_MALLOC(ctstar_table->ctstar_values, local_n_ct, cs_real_t);

    cs_real_t *ctstar_speeds = ctstar_table->ctstar_speeds;
    cs_real_t *ctstar_values = ctstar_table->ctstar_values;

    //read and parse this file
    stream_ct = fopen(ctstar_file_names[WT_type_count], "r");
    char *tok3;
    const char sep3[4] = ",";
    cs_lnum_t line_count = -1;
    cs_lnum_t i=-1;

    while (fgets(line, 1024, stream_ct) != NULL){
      i++;
      if (++line_count >= local_n_ct)
  	    break;

      // parse line
      tok3 = strtok(line, sep3);
      if (tok3 == NULL)
  	    continue;
      //IEC 61400-12 air density correction
      ctstar_speeds[i]=atof(tok3)*pow(standard_rho/cs_glob_fluid_properties->ro0,1./3.);

      tok3 = strtok(NULL, sep3);
      if (tok3 == NULL)
  	    continue;
      ctstar_values[i]=atof(tok3);
    }
    fclose(stream_ct);

    //cp
    //calculate number of lines in the file
    FILE* stream_cp = fopen(cpstar_file_names[WT_type_count], "r");
    char *tok4;
    const char sep4[4] = ",";
    cs_lnum_t local_n_cp = 0;
    //read header
    //fgets(line, 1024, stream);
    while (fgets(line, 1024, stream_cp) != NULL){
      // try to parse line
      tok4 = strtok(line, sep4);
      if (tok4 != NULL)
  	    local_n_cp++;
    }
    cpstar_table->n_cpstar = local_n_cp;
    fclose(stream_cp);

    //allocate memory for all members of structure
    BFT_MALLOC(cpstar_table->cpstar_speeds, local_n_ct, cs_real_t);
    BFT_MALLOC(cpstar_table->cpstar_values, local_n_ct, cs_real_t);

    cs_real_t *cpstar_speeds = cpstar_table->cpstar_speeds;
    cs_real_t *cpstar_values = cpstar_table->cpstar_values;

    //read and parse this file
    stream_cp = fopen(cpstar_file_names[WT_type_count], "r");
    char *tok5;
    const char sep5[4] = ",";
    line_count = -1;
    i=-1;

    while (fgets(line, 1024, stream_cp) != NULL){
      i++;
      if (++line_count >= local_n_cp)
  	    break;

      // parse line
      tok5 = strtok(line, sep5);
      if (tok5 == NULL)
  	    continue;
      //IEC 61400-12 air density correction
      cpstar_speeds[i]=atof(tok5)*pow(standard_rho/cs_glob_fluid_properties->ro0,1./3.);

      tok5 = strtok(NULL, sep5);
      if (tok5 == NULL)
  	    continue;
      cpstar_values[i]=atof(tok5);
    }
    fclose(stream_cp);
    _cs_ctstar_table_array[WT_type_count] = ctstar_table;
    _cs_cpstar_table_array[WT_type_count] = cpstar_table;
  }
}

cs_wind_farm *cs_glob_wind_farm = &_wind_farm;

/*------------------cs_wind_farm_create_from_file-----------------------------*/
cs_wind_farm *
cs_wind_farm_create_from_file(const char *WT_file_name)
{
  //calculate number of lines in the file
  FILE* stream = fopen(WT_file_name, "r");

  if (stream == NULL)
    return NULL;

  char line[1024];
  char *tok0;
  const char sep0[4] = ";";
  cs_lnum_t n_WT = 0;
  while (fgets(line, 1024, stream) != NULL){
    // try to parse line
    if (line[0] == '#') {
    }
    else {
      tok0 = strtok(line, sep0);
      if (tok0 != NULL)
        n_WT++;
    }
  }
  fclose(stream);

  cs_glob_wind_farm->n_WT = n_WT;

  cs_glob_wind_farm->total_power = 0.0;

  BFT_MALLOC(cs_glob_wind_farm->WT_file_name, strlen(WT_file_name) + 1, char);
  BFT_MALLOC(cs_glob_wind_farm->WT_coords, n_WT, cs_real_3_t);
  BFT_MALLOC(cs_glob_wind_farm->WT_diameters, n_WT, cs_real_t);
  BFT_MALLOC(cs_glob_wind_farm->WT_types, n_WT, cs_lnum_t);
  BFT_MALLOC(cs_glob_wind_farm->WT_volume, n_WT, cs_real_t);
  BFT_MALLOC(cs_glob_wind_farm->WT_uvect, n_WT, cs_real_3_t);
  BFT_MALLOC(cs_glob_wind_farm->WT_u, n_WT, cs_real_t);
  BFT_MALLOC(cs_glob_wind_farm->WT_rho, n_WT, cs_real_t);
  BFT_MALLOC(cs_glob_wind_farm->WT_hub_u, n_WT, cs_real_t);
  BFT_MALLOC(cs_glob_wind_farm->WT_dir, n_WT, cs_real_t);
  BFT_MALLOC(cs_glob_wind_farm->WT_ra_dir, n_WT, cs_real_t);
  BFT_MALLOC(cs_glob_wind_farm->WT_ctstar, n_WT, cs_real_t);
  BFT_MALLOC(cs_glob_wind_farm->WT_cpstar, n_WT, cs_real_t);
  BFT_MALLOC(cs_glob_wind_farm->WT_thrust, n_WT, cs_real_t);
  BFT_MALLOC(cs_glob_wind_farm->WT_power, n_WT, cs_real_t);

  //read and parse the file
  stream = fopen(WT_file_name, "r");
  char *tok;
  const char sep[4] = ",";
  cs_lnum_t line_count = -1;
  cs_lnum_t i=-1;
  while (fgets(line, 1024, stream) != NULL){
    i++;
    if (++line_count >= n_WT+1)
      break;
    //skip header
    if(i>0) {
      // parse line
      tok = strtok(line, sep);
      if (tok == NULL)
	      continue;
      cs_glob_wind_farm->WT_coords[i-1][0]=atof(tok);

      tok = strtok(NULL, sep);
      if (tok == NULL)
	      continue;
      cs_glob_wind_farm->WT_coords[i-1][1]=atof(tok);

      tok = strtok(NULL, sep);
      if (tok == NULL)
	      continue;
      cs_glob_wind_farm->WT_coords[i-1][2]=atof(tok);

      tok = strtok(NULL, sep);
      if (tok == NULL)
      	continue;
      cs_glob_wind_farm->WT_diameters[i-1]=atof(tok);

      tok = strtok(NULL, sep);
      if (tok == NULL)
      	continue;
      cs_glob_wind_farm->WT_types[i-1]=atoi(tok);

      cs_glob_wind_farm->WT_volume[i-1] = 0;
      cs_glob_wind_farm->WT_uvect[i-1][0] = 0;
      cs_glob_wind_farm->WT_uvect[i-1][1] = 0;
      cs_glob_wind_farm->WT_uvect[i-1][2] = 0;
      cs_glob_wind_farm->WT_u[i-1] = 0;
      cs_glob_wind_farm->WT_rho[i-1] = 0;
      cs_glob_wind_farm->WT_hub_u[i-1] = 0;
      cs_glob_wind_farm->WT_dir[i-1] = 0;
      cs_glob_wind_farm->WT_ra_dir[i-1] = 0;
      cs_glob_wind_farm->WT_ctstar[i-1] = 0;
      cs_glob_wind_farm->WT_cpstar[i-1] = 0;
      cs_glob_wind_farm->WT_thrust[i-1] = 0;
      cs_glob_wind_farm->WT_power[i-1] = 0;

    }
  }
  fclose(stream);

  //get number of turbine types
  cs_lnum_t n_WT_types=-1;
  for (cs_lnum_t WT_count=0; WT_count < n_WT; WT_count ++){
    if(cs_glob_wind_farm->WT_types[WT_count] > n_WT_types) {
      n_WT_types = cs_glob_wind_farm->WT_types[WT_count];
    }
  }
  cs_glob_wind_farm->n_WT_types = n_WT_types;

  _coeff_tables_create(n_WT_types);

  return cs_glob_wind_farm;
}
/*-------------------end of function cs_wind_farm_create_from_file------------*/



/*----------------function cs_wind_farm_free----------------------------------*/
cs_wind_farm *
cs_wind_farm_free(cs_wind_farm *wind_farm)
{
  //free memory for the members of the wind farm structure
  BFT_FREE(wind_farm->WT_file_name);
  BFT_FREE(wind_farm->WT_coords);
  BFT_FREE(wind_farm->WT_diameters);
  BFT_FREE(wind_farm->WT_types);
  BFT_FREE(wind_farm->WT_volume);
  BFT_FREE(wind_farm->WT_uvect);
  BFT_FREE(wind_farm->WT_u);
  BFT_FREE(wind_farm->WT_rho);
  BFT_FREE(wind_farm->WT_hub_u);
  BFT_FREE(wind_farm->WT_dir);
  BFT_FREE(wind_farm->WT_ra_dir);
  BFT_FREE(wind_farm->WT_ctstar);
  BFT_FREE(wind_farm->WT_cpstar);
  BFT_FREE(wind_farm->WT_thrust);
  BFT_FREE(wind_farm->WT_power);

  for (cs_lnum_t WT_type_count=0; WT_type_count < wind_farm->n_WT_types; WT_type_count ++){
    cs_ctstar_table *ctstar_table = _cs_ctstar_table_array[WT_type_count];
    cs_cpstar_table *cpstar_table = _cs_cpstar_table_array[WT_type_count];
    BFT_FREE(ctstar_table->ctstar_speeds);
    BFT_FREE(ctstar_table->ctstar_values);
    BFT_FREE(cpstar_table->cpstar_speeds);
    BFT_FREE(cpstar_table->cpstar_values);
    BFT_FREE(ctstar_table);
    BFT_FREE(cpstar_table);
  }
  BFT_FREE(_cs_ctstar_table_array);
  BFT_FREE(_cs_cpstar_table_array);
  return wind_farm;
}
/*-------------------end of function cs_wind_farm_free------------------------*/

/*---------------function cs_wind_farm_compute_source_term--------------------*/
void
cs_wind_farm_compute_source_term()
{
  const cs_mesh_t *m = cs_glob_mesh;
  cs_mesh_quantities_t *mq = cs_glob_mesh_quantities;
  const cs_real_3_t *cell_cen
    = (const cs_real_3_t *)mq->cell_cen;
  const cs_real_t *cell_vol = mq->cell_vol;
  cs_lnum_t n_WT = cs_glob_wind_farm->n_WT;

  const cs_real_t *cpro_rom = CS_F_(rho)->val;
  const cs_real_3_t *vel = (const cs_real_3_t *)CS_F_(vel)->val;

  cs_real_3_t *WF_source_term = (cs_real_3_t *)
                                        cs_field_by_name("WF_source_term")->val;
  cs_real_t *st_coeff = cs_field_by_name("source_term_coeff")->val;

  cs_real_3_t *WT_coords = cs_glob_wind_farm->WT_coords;
  cs_real_t *WT_volume = cs_glob_wind_farm->WT_volume;
  cs_real_t *WT_diameters = cs_glob_wind_farm->WT_diameters;
  cs_lnum_t *WT_types = cs_glob_wind_farm->WT_types;
  cs_real_3_t *WT_uvect= cs_glob_wind_farm->WT_uvect;
  cs_real_t *WT_u = cs_glob_wind_farm->WT_u;
  cs_real_t *WT_rho = cs_glob_wind_farm->WT_rho;
  cs_real_t *WT_hub_u = cs_glob_wind_farm->WT_hub_u;
  cs_real_t *WT_dir = cs_glob_wind_farm->WT_dir;
  cs_real_t *WT_ra_dir = cs_glob_wind_farm->WT_ra_dir;
  cs_real_t *WT_ctstar = cs_glob_wind_farm->WT_ctstar;
  cs_real_t *WT_cpstar = cs_glob_wind_farm->WT_cpstar;
  cs_real_t *WT_thrust = cs_glob_wind_farm->WT_thrust;
  cs_real_t *WT_power = cs_glob_wind_farm->WT_power;
  //
  cs_real_t AD_mesh_cell_size =
                 cs_notebook_parameter_value_by_name("AD_mesh_cell_size");
  cs_real_t AD_half_rotor_thickness = 1.2*AD_mesh_cell_size;
  cs_real_t select_dist = 1.2*AD_mesh_cell_size;
  /******************LOOP ON TURBINES*************************/
  cs_lnum_t start_WT_count, end_WT_count;
  if (cs_notebook_parameter_value_by_name("isol")>0) {
    start_WT_count = cs_notebook_parameter_value_by_name("isol")-1;
    end_WT_count = cs_notebook_parameter_value_by_name("isol");
  }
  else {
    start_WT_count = 0;
    end_WT_count = n_WT;
  }

  const cs_real_t dpi = atan(1.0)*4.0;
  cs_real_t interp_value;
  cs_real_t wind_dir=cs_glob_atmo_option->meteo_angle;
  cs_real_t ra_wind_dir= (270.0-wind_dir)*dpi/180.0; //radian absolute wind_dir
  //
  for (cs_lnum_t WT_count=start_WT_count; WT_count < end_WT_count; WT_count ++){
    //
    cs_real_t WT_d = WT_diameters[WT_count];
    cs_real_t WT_radius = WT_d/2.;
    cs_real_t WT_surf = dpi*pow(WT_radius,2);
    cs_lnum_t local_WT_type = WT_types[WT_count]-1;
    cs_ctstar_table *ctstar_table = _cs_ctstar_table_array[local_WT_type];
    cs_cpstar_table *cpstar_table = _cs_cpstar_table_array[local_WT_type];
    cs_lnum_t local_n_ct = ctstar_table->n_ctstar;
    cs_lnum_t local_n_cp = cpstar_table->n_cpstar;
    //Get AD zone and define source term in it
    /***********************************************************/
    //get Actuator Disk zone by name
    char name[128];
    sprintf(name,"turbine_%d",WT_count+1);
    const cs_zone_t *z = cs_volume_zone_by_name_try(name);
    //
    WT_volume[WT_count]=0.0;
    WT_u[WT_count]=0.0;
    WT_uvect[WT_count][0]=0.0;
    WT_uvect[WT_count][1]=0.0;
    WT_uvect[WT_count][2]=0.0;
    WT_rho[WT_count]=0.0;
    WT_dir[WT_count]=0.0;
    WT_ra_dir[WT_count]=0.0;
    WT_hub_u[WT_count]=0.0;

    //Compute local Wind direction at Turbine hub height
    cs_real_3_t closest_uvect;
    cs_lnum_t closest_id;
    int closest_id_rank;
    cs_real_t xyz_ref[3] = {WT_coords[WT_count][0],
                            WT_coords[WT_count][1],
                            WT_coords[WT_count][2]};
    cs_geom_closest_point(m->n_cells,
  		                    (const cs_real_3_t *)(mq->cell_cen),
  		                    xyz_ref,
  		                    &closest_id,
  		                    &closest_id_rank);

    if (closest_id_rank == cs_glob_rank_id) {
     	closest_uvect[0] = vel[closest_id][0];
     	closest_uvect[1] = vel[closest_id][1];
     	closest_uvect[2] = vel[closest_id][2];
    }
    cs_parall_bcast(closest_id_rank, 3, CS_REAL_TYPE, &closest_uvect);
    WT_hub_u[WT_count]=sqrt(pow(closest_uvect[0],2)+pow(closest_uvect[1],2)
                             +pow(closest_uvect[2],2));
    WT_dir[WT_count]= 270.0 - atan2(closest_uvect[1],closest_uvect[0])*180.0/dpi;

    WT_ra_dir[WT_count]=atan2(closest_uvect[1],closest_uvect[0]);
    cs_real_t mast_center_coords[3];

    mast_center_coords[0] = WT_coords[WT_count][0];
    mast_center_coords[1] = WT_coords[WT_count][1];
    mast_center_coords[2] = 0.0;

    cs_real_t yaw_angle = 0.0;
    cs_real_t tilt_angle = 0.0;
    cs_real_t roll_angle = 0.0;

    if(cs_notebook_parameter_value_by_name("control")>0) {
      yaw_angle = 270 - WT_dir[WT_count];
      tilt_angle = 0;
      roll_angle = 0;
    }
    else {
      yaw_angle = 270 - wind_dir;
      tilt_angle = 0;
      roll_angle = 0;
    }

    //
    for (cs_lnum_t nc = 0; nc < z->n_elts; nc++) {
      cs_lnum_t c_id = z->elt_ids[nc];
      //
      cs_real_t x_dist,y_dist,z_dist;

      cs_real_t rotated_WT_center_coords[3];
      cs_real_t WT_center_coords[3];

      WT_center_coords[0] = WT_coords[WT_count][0];
      WT_center_coords[1] = WT_coords[WT_count][1];
      WT_center_coords[2] = WT_coords[WT_count][2];
      rotate2(WT_center_coords, tilt_angle, yaw_angle, roll_angle, \
        mast_center_coords, rotated_WT_center_coords);
      //
      cs_real_t deyawed_cell_cen[3];
      rotate2(cell_cen[c_id], 360.0 - tilt_angle, 360.0 - yaw_angle, \
      360.0 - roll_angle, mast_center_coords, deyawed_cell_cen);

      x_dist = deyawed_cell_cen[0] - WT_center_coords[0];
      y_dist = deyawed_cell_cen[1] - WT_center_coords[1];
      z_dist = deyawed_cell_cen[2] - WT_center_coords[2];

      cs_real_t cell_radius = sqrt(pow(y_dist,2)+pow(z_dist,2));
      //
      st_coeff[c_id]=0.0;
      if (cs_notebook_parameter_value_by_name("st_method")==0) {
        if (cell_radius<=WT_radius &&
          fabs(x_dist)<=AD_half_rotor_thickness) {
          st_coeff[c_id]=1.0;
        }
        else if (cell_radius>WT_radius+select_dist ||
          fabs(x_dist)>AD_half_rotor_thickness+select_dist) {
          st_coeff[c_id]= 0.0;
        }
        else {
          st_coeff[c_id]= 0.0; //worked with 0.2
        }
        /* else if (fabs(deyawed_x_dist)<=AD_half_rotor_thickness + select_dist)  { */
        /* st_coeff[c_id]=fmax(a_xdist*fabs(deyawed_x_dist) + b_xdist,0.0); */
        /* } */
        /* else if (cell_radius<=WT_radius + select_dist) { */
        /*   st_coeff[c_id]= fmax(a_radius*cell_radius + b_radius,0.0); */
        /* } */
      }
      WT_volume[WT_count] = WT_volume[WT_count]+st_coeff[c_id]*cell_vol[c_id];
      WT_uvect[WT_count][0] += vel[c_id][0]*(st_coeff[c_id]*cell_vol[c_id]);
      WT_uvect[WT_count][1] += vel[c_id][1]*(st_coeff[c_id]*cell_vol[c_id]);
      WT_uvect[WT_count][2] += vel[c_id][2]*(st_coeff[c_id]*cell_vol[c_id]);
      WT_rho[WT_count] += cpro_rom[c_id]*(st_coeff[c_id]*cell_vol[c_id]);
    }
    cs_parall_sum(1, CS_REAL_TYPE, &WT_volume[WT_count]);
    cs_parall_sum(3, CS_REAL_TYPE, &WT_uvect[WT_count]);
    cs_parall_sum(1, CS_REAL_TYPE, &WT_rho[WT_count]);
    //
    WT_uvect[WT_count][0] = WT_uvect[WT_count][0]/WT_volume[WT_count];
    WT_uvect[WT_count][1] = WT_uvect[WT_count][1]/WT_volume[WT_count];
    WT_uvect[WT_count][2] = WT_uvect[WT_count][2]/WT_volume[WT_count];
    WT_rho[WT_count] = WT_rho[WT_count]/WT_volume[WT_count];
    //WT_dir[WT_count]=270.0 - atan2(WT_uvect[WT_count][1],WT_uvect[WT_count][0])*180.0/dpi;
    //

    /* Compute disk normal velocity */
    cs_real_t disk_normal_coords[3];

    cs_real_t wind_normal_coords[3] = {cos(ra_wind_dir), sin(ra_wind_dir), 0};
    cs_real_t origin[3] = {0.0, 0.0, 0.0};
    cs_real_t unitx[3] = {1.0, 0.0, 0.0};
    rotate2(unitx, tilt_angle, yaw_angle, roll_angle, origin, disk_normal_coords);
    for(int n=0; n<3; n++){
       WT_u[WT_count] += disk_normal_coords[n]*WT_uvect[WT_count][n];
     }

     interp_value = interp_ct_or_cp(ctstar_table->ctstar_values, ctstar_table->ctstar_speeds, \
        local_n_ct, WT_u[WT_count]);
     WT_ctstar[WT_count] = interp_value;
     interp_value = interp_ct_or_cp(cpstar_table->cpstar_values, cpstar_table->cpstar_speeds, \
        local_n_cp, WT_u[WT_count]);
     WT_cpstar[WT_count] = interp_value;

     WT_thrust[WT_count] = 0.5*WT_rho[WT_count]*WT_surf*WT_ctstar[WT_count]*cs_math_pow2(WT_u[WT_count]);
     WT_power[WT_count] = 0.5*WT_rho[WT_count]*WT_surf*WT_cpstar[WT_count]*cs_math_pow2(WT_u[WT_count])*WT_u[WT_count];
     for (cs_lnum_t nc = 0; nc < z->n_elts; nc++) {
       cs_lnum_t c_id = z->elt_ids[nc];
       //WT_ctstar as u is taken at wind turbine position, not freestream
       cs_real_t AD_coeff=0.5*WT_rho[WT_count]*WT_surf*WT_ctstar[WT_count];

       //_st_exp[c_id][0]    += -(cell_vol[c_id]*st_coeff[c_id])*AD_coeff*cs_math_pow2(WT_u[WT_count])/WT_volume[WT_count];
       //_st_exp[c_id][1]    += 0.0;
   //
       WF_source_term[c_id][0] = -disk_normal_coords[0]*(cell_vol[c_id]*st_coeff[c_id])*AD_coeff*pow(WT_u[WT_count],2)/WT_volume[WT_count];
       WF_source_term[c_id][1] = -disk_normal_coords[1]*(cell_vol[c_id]*st_coeff[c_id])*AD_coeff*pow(WT_u[WT_count],2)/WT_volume[WT_count];
       WF_source_term[c_id][2] = -disk_normal_coords[2]*(cell_vol[c_id]*st_coeff[c_id])*AD_coeff*pow(WT_u[WT_count],2)/WT_volume[WT_count];
    }
  }

  cs_glob_wind_farm->total_power = 0.0;
  for (cs_lnum_t WT_count=0; WT_count < n_WT; WT_count ++){
    cs_glob_wind_farm->total_power += cs_glob_wind_farm->WT_power[WT_count];
  }

  return;
}
/*------------end of function cs_wind_farm_compute_source_term----------------*/

/*---------------------function cs_wind_farm_write----------------------------*/
void
cs_wind_farm_write()
{
  cs_time_step_t *ts = cs_get_glob_time_step();
  char power_file_name[1024];
  sprintf(power_file_name,"power_iter%d.csv",ts->nt_cur);
  FILE* stream10 = fopen(power_file_name, "w");
  fprintf(stream10,
          "#total wind farm power calculated with cp_star, is %.4f\n",
          cs_glob_wind_farm->total_power);
  fprintf(stream10, "#wind_turbine,xhub,yhub,zhub,turbine_diameter,turbine_type,ux,uy,uz,u,u_hub,dir,rho,ct*,cp*,thrust,power_cpstar\n");
    for (cs_lnum_t WT_count=0; WT_count < cs_glob_wind_farm->n_WT; WT_count ++){
      fprintf(stream10,"%d,%.4f,%.4f,%.4f,%.4f,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n",
              WT_count+1,
              cs_glob_wind_farm->WT_coords[WT_count][0],
              cs_glob_wind_farm->WT_coords[WT_count][1],
              cs_glob_wind_farm->WT_coords[WT_count][2],
              cs_glob_wind_farm->WT_diameters[WT_count],
              cs_glob_wind_farm->WT_types[WT_count],
              cs_glob_wind_farm->WT_uvect[WT_count][0],
              cs_glob_wind_farm->WT_uvect[WT_count][1],
              cs_glob_wind_farm->WT_uvect[WT_count][2],
              cs_glob_wind_farm->WT_u[WT_count],
              cs_glob_wind_farm->WT_hub_u[WT_count],
              cs_glob_wind_farm->WT_dir[WT_count],
              cs_glob_wind_farm->WT_rho[WT_count],
              cs_glob_wind_farm->WT_ctstar[WT_count],
              cs_glob_wind_farm->WT_cpstar[WT_count],
              cs_glob_wind_farm->WT_thrust[WT_count],
              cs_glob_wind_farm->WT_power[WT_count]);
    }
    fclose(stream10);
  return;
}
/*--------------------end of function cs_wind_farm_write----------------------*/

void
rotate2(cs_real_t *coords,
        cs_real_t tilt_angle,
        cs_real_t yaw_angle,
        cs_real_t roll_angle,
        cs_real_t *mean_coords,
        cs_real_t *rotated_coords)
{
  //
  const cs_real_t dpi = atan(1.0)*4.0;
  cs_real_t ra_tilt_angle= tilt_angle*dpi/180.0; //y-rotation
  cs_real_t ra_yaw_angle= yaw_angle*dpi/180.0; //z-rotation
  cs_real_t ra_roll_angle= roll_angle*dpi/180.0; //x-rotation

  //centering
  cs_real_t initial_coords[3];
  initial_coords[0] = coords[0]-mean_coords[0];
  initial_coords[1] = coords[1]-mean_coords[1];
  initial_coords[2] = coords[2]-mean_coords[2];

  //x-rotation
  cs_real_t rolled_coords[3];
  rolled_coords[0] = initial_coords[0];
  rolled_coords[1] = cos(ra_roll_angle)*initial_coords[1]
                            - sin(ra_roll_angle)*initial_coords[2];
  rolled_coords[2] = sin(ra_roll_angle)*initial_coords[1]
                            + cos(ra_roll_angle)*initial_coords[2];


  //y-rotation
  cs_real_t tilted_coords[3];
  tilted_coords[0] = cos(ra_tilt_angle)*rolled_coords[0]
                              + sin(ra_tilt_angle)*rolled_coords[2];
  tilted_coords[1] = rolled_coords[1];
  tilted_coords[2] = -sin(ra_tilt_angle)*rolled_coords[0]
                              + cos(ra_tilt_angle)*rolled_coords[2];


  //z-rotation
  cs_real_t yawed_coords[3];
  yawed_coords[0] = cos(ra_yaw_angle)*tilted_coords[0]
                                 - sin(ra_yaw_angle)*tilted_coords[1];
  yawed_coords[1] = sin(ra_yaw_angle)*tilted_coords[0]
                                 + cos(ra_yaw_angle)*tilted_coords[1];
  yawed_coords[2] = tilted_coords[2];

  //final
  rotated_coords[0] = yawed_coords[0] + mean_coords[0];
  rotated_coords[1] = yawed_coords[1] + mean_coords[1];
  rotated_coords[2] = yawed_coords[2] + mean_coords[2];
}

cs_real_t
interp_ct_or_cp(cs_real_t *ct_or_cp_values, cs_real_t *ct_or_cp_speeds, \
		cs_lnum_t number_of_values, cs_real_t disk_velocity) {
  //
  cs_real_t val_1, val_2, u_1, u_2;
  cs_real_t interpolated_ct_or_cp = 0;
  val_1=ct_or_cp_values[0];
  val_2=ct_or_cp_values[number_of_values-1];
  u_1 = ct_or_cp_speeds[0];
  u_2 = ct_or_cp_speeds[number_of_values-1];
  if(disk_velocity<u_1 || disk_velocity>u_2) {
    interpolated_ct_or_cp = 0.0;
  }
  else {
    val_2 = ct_or_cp_values[0];
    u_2 = ct_or_cp_speeds[0];
    for (cs_lnum_t ct_count=1; ct_count < number_of_values; ct_count ++){
      val_1= val_2;
      u_1=u_2;
      //
      val_2= ct_or_cp_values[ct_count];
      u_2=ct_or_cp_speeds[ct_count];
      //
      if(disk_velocity>=u_1 && disk_velocity<=u_2) {
      	/* //linear interpolation between two values;*/
      	cs_real_t a = (val_2 - val_1)/(u_2-u_1);
      	cs_real_t b = val_1 - a*u_1;
      	interpolated_ct_or_cp = a*disk_velocity+b;
      	break;
      }
    }
  }
  return interpolated_ct_or_cp;
}
