#include<iostream>

#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"

void update_coeffs(){
  
  // update nodal velocity and  position //
  for (int node_gid = 0; node_gid < mesh.num_nodes(); node_gid++) {
     
    for (int dim = 0; dim < 3; dim++){
      node.vel(0,node_gid, dim) = node.vel(1,node_gid,dim);
    }// end loop over dim      
    
    for (int dim = 0; dim < 3; dim++){
      node.coords(0, node_gid, dim) = node.coords(1, node_gid, dim);
      //mesh.node_coords(node_gid, dim) = node.coords(1, node_gid, dim);      
    }// end loop over dim
    
  }// end loop over node_gid 


#pragma omp simd     
  for (int cell_gid = 0; cell_gid < mesh.num_cells(); cell_gid++){

    cell_state.ie( 0, cell_gid ) = cell_state.ie( 1, cell_gid );
    cell_state.total_energy( 0, cell_gid ) = cell_state.total_energy( 1, cell_gid );
    
  }

}// end update


