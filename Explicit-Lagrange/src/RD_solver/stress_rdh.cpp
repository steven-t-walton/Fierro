/* Compute stress and save to cell_state.stress */

#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"

using namespace utils;

void get_stress_tensor(int t_step){

  for (int gauss_gid = 0; gauss_gid < mesh.num_gauss_pts(); gauss_gid++){
    for (int dim_i = 0; dim_i < mesh.num_dim(); dim_i++){
      for (int dim_j = 0; dim_j < mesh.num_dim(); dim_j++){
        elem_state.stress_tensor(t_step, gauss_gid, dim_i, dim_j) = 0.0;
      }
    }
  }

  for (int gauss_gid = 0; gauss_gid < mesh.num_gauss_pts(); gauss_gid++){
    for (int dim_i = 0; dim_i < mesh.num_dim(); dim_i++){
      for (int dim_j = 0; dim_j < mesh.num_dim(); dim_j++){

     	if ( dim_i == dim_j ){
	  elem_state.stress_tensor(t_step, gauss_gid, dim_i, dim_j) = -1.0*mat_pt.pressure(gauss_gid);
	}// end if
	else if ( dim_i != dim_j){
	  elem_state.stress_tensor(t_step, gauss_gid, dim_i, dim_j) = 0.0;
	}// end else if

      }// end loop over dim_j
    }// end loop over dim_i
  }// end loop over gaus_gid

}// end get_stress
