#include<iostream>

#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"

void update_tn(){
  
  // update nodal velocity and  position //
  for (int node_gid = 0; node_gid < mesh.num_nodes(); node_gid++) {
    for (int dim = 0; dim < mesh.num_dim(); dim++){
      
      node.vel(0, node_gid, dim) = node.vel(num_correction_steps, node_gid, dim);
      node.coords(0, node_gid, dim) = node.coords(num_correction_steps, node_gid, dim);
      
      node.vel(1, node_gid, dim) = node.vel(num_correction_steps, node_gid, dim);
      node.coords(1, node_gid, dim) = node.coords(num_correction_steps, node_gid, dim);
    
    }// end loop over dim
  }// end loop over node_gid 

#pragma omp simd     
  for (int cell_gid = 0; cell_gid < mesh.num_cells(); cell_gid++){

    cell_state.ie( 0, cell_gid ) = cell_state.ie( 1, cell_gid );
    cell_state.total_energy( 0, cell_gid ) = cell_state.total_energy( 1, cell_gid );
    
  }
  
  for (int gauss_gid = 0; gauss_gid < mesh.num_gauss_pts(); gauss_gid++){
    mat_pt.sie(0, gauss_gid) = mat_pt.sie(num_correction_steps, gauss_gid);
    mat_pt.sie(1, gauss_gid) = mat_pt.sie(num_correction_steps, gauss_gid);
  }
  

}// end update


