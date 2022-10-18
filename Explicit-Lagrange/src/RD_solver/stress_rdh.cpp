/* Compute stress and save to cell_state.stress */

#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"

using namespace utils;

void get_stress(){
#pragma omp simd
  auto pressure = ViewCArray <real_t> (&cell_state.pressure(0), mesh.num_cells() );
  for ( int cell_gid = 0; cell_gid < mesh.num_cells(); cell_gid++){
    for (int j = 0; j < mesh.num_dim(); j++){
      for (int i = 0; i < mesh.num_dim(); i++){
        if (i == j){
          cell_state.stress(1, cell_gid, i, j) = -pressure(cell_gid);  
	}// end if
	else if ( i != j ){
	  cell_state.stress(1,cell_gid, i, j) = 0.0;
	};
      }// end loop over i
    }// end loop over j

  }// end loop over cell_gid
}// end sub-routine
