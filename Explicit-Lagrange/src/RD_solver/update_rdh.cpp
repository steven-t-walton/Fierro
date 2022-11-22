#include<iostream>

#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"

void update_coeffs(){
  
  //int update = num_correction_steps;

  // Update velocity control coeffs
  for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
    for (int vertex = 0; vertex <  ref_elem.num_basis(); vertex++){
      int node_lid = elem.vert_node_map( vertex );
      int node_gid = mesh.nodes_in_elem( elem_gid, node_lid);
      for (int dim = 0; dim < 3; dim++){
        node.vel(0,node_gid, dim) = node.vel(num_correction_steps,node_gid,dim);
      }// end loop over dim      
    }
  }
  //get_control_coeffs();

  // update nodal position //
  for (int node_gid = 0; node_gid < mesh.num_nodes(); node_gid++) {
    
    for (int dim = 0; dim < 3; dim++){
      node.coords(0, node_gid, dim) = node.coords(1, node_gid, dim);	
    }// end loop over dim
    
  }// end loop over node_gid 
  
  // update ie, stress and total energy //
#pragma omp simd     
  for (int cell_gid = 0; cell_gid < mesh.num_cells(); cell_gid++){

    cell_state.ie( 0, cell_gid ) = cell_state.ie( 1, cell_gid );
  /*     
    for (int dim_j = 0; dim_j < mesh.num_dim(); dim_j++){
      for (int dim_i = 0; dim_i < mesh.num_dim(); dim_i++){
        cell_state.stress(0, cell_gid, dim_i, dim_j) = cell_state.stress(1, cell_gid, dim_i, dim_j);
      }
    }
  */  
    cell_state.total_energy( 0, cell_gid ) = cell_state.total_energy( 1, cell_gid );
  }

  // update control coeffs //
/*
  for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
    for (int basis = 0; basis < ref_elem.num_basis(); basis++){
      for (int dim = 0; dim < mesh.num_dim(); dim++){
        elem_state.BV_vel_coeffs( 0, elem_gid, basis, dim ) = elem_state.BV_vel_coeffs( num_correction_steps, elem_gid, basis, dim );
      }// end loop over dim
    }// end loop over basis
  }// end loop over elem_gid
*/
}// end update


