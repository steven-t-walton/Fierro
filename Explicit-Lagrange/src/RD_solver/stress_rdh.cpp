/* Compute stress and save to cell_state.stress */

#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"

using namespace utils;

void get_stress(int t_step){

  for ( int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
   
    for ( int cell_lid = 0; cell_lid < mesh.num_cells_in_elem(); cell_lid++){
      int cell_gid = mesh.cells_in_elem(elem_gid, cell_lid);
      // create view of cell_state.stress() //
      //auto stress = ViewCArray <real_t> (&cell_state.stress(0,cell_gid,0,0), num_correction_steps, mesh.num_cells_in_elem(), mesh.num_dim(), mesh.num_dim());
      for (int j = 0; j < mesh.num_dim(); j++){
        for (int i = 0; i < mesh.num_dim(); i++){
          if (i == j){
            cell_state.stress(t_step, cell_gid, i, j) = -cell_state.pressure(cell_gid);  
          }
	  else if ( i !=j ){
	    cell_state.stress(t_step, cell_gid, i, j) = 0.0;
	  }
        }
      } 
    }
  }
}// end sub-routine
