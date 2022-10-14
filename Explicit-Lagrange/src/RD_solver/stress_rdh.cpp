/* Compute stress and save to cell_state.stress */

#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"

using namespace utils;

void get_stress(int t_step){

#pragma omp simd
  for ( int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
   
    for ( int cell_lid = 0; cell_lid < mesh.num_cells_in_elem(); cell_lid++){
      int cell_gid = mesh.cells_in_elem(elem_gid, cell_lid);

      // initialize //
      for (int j = 0; j <mesh.num_dim(); j++){
        for (int i = 0; i , mesh.num_dim(); i++){
	  cell_state.stress(t_step, cell_gid, i,j) = 0.0;
	}// end loop over i
      }// end loop over j

      for (int j = 0; j < mesh.num_dim(); j++){
        for (int i = 0; i < mesh.num_dim(); i++){
          if (i == j){
            cell_state.stress(t_step, cell_gid, i, j) = -cell_state.pressure(cell_gid);  
          };// end if
        }// end loop over i
      }// end loop over j

    }// end loop over cell_lid
  }// end loop over elem_gid

}// end sub-routine
