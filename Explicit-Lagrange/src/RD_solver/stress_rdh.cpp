/* Compute stress and save to cell_state.stress */

#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"

using namespace utils;

void get_stress(int t_step){
  if ( t_step < num_correction_steps ){
#pragma omp simd
  for ( int cell_gid = 0; cell_gid < mesh.num_cells(); cell_gid++){
   
    /*
    // initialize //
    for (int j = 0; j <mesh.num_dim(); j++){
      for (int i = 0; i , mesh.num_dim(); i++){
        cell_state.stress(t_step, cell_gid, i,j) = 0.0;
      }// end loop over i
    }// end loop over j
    */

    for (int j = 0; j < mesh.num_dim(); j++){
      for (int i = 0; i < mesh.num_dim(); i++){
        if (i == j){
          cell_state.stress(t_step, cell_gid, i, j) = -cell_state.pressure(cell_gid);  
        };// end if
      }// end loop over i
    }// end loop over j

  }// end loop over cell_gid
  };// end if
}// end sub-routine
