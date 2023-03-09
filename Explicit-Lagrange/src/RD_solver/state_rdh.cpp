/* state.cpp */

#include<iostream>
#include<math.h>

#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"

#define PI 3.14159265

using namespace utils;

void get_state(){
// /*

  for( int cell_gid = 0; cell_gid < mesh.num_cells(); cell_gid++){
    //cell_properties(cell_gid);
    real_t elem_coords_x = 0.0;
    real_t elem_coords_y = 0.0;
    elem_coords_x = mesh.cell_coords(cell_gid,0);
    elem_coords_y = mesh.cell_coords(cell_gid,1);

    cell_state.density(cell_gid) = 0.0; 
   
    cell_state.pressure(cell_gid) = 0.25*( cos(2.0*3.141592653589 * elem_coords_x )  + cos(2.0*3.141592653589 * elem_coords_y ) ) + 1.0;
    
    cell_state.ie(1,cell_gid) = 0.0;    
    for (int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_cell(); gauss_lid++){
      int gauss_gid = mesh.gauss_in_cell(cell_gid, gauss_lid);
      cell_state.ie(1, cell_gid) += 0.125 * mat_pt.sie(num_correction_steps,gauss_gid);
      cell_state.density(cell_gid) += 0.125*mat_pt.density(gauss_gid);
    }

    real_t ie = 0.0;
    real_t ke = 0.0;
    track_rdh(ie,ke);
    
    cell_state.total_energy(1,cell_gid) = ie + ke; 
    
    cell_state.cs(cell_gid) =
        material[cell_state.mat_id(cell_gid)].eos_func(
            sspd,
            cell_state.mat_id(cell_gid),
            cell_state.density(cell_gid),
            cell_state.ie(1, cell_gid)
            );
  }
}
