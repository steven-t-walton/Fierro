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
  // update ie, stress and total energy //
#pragma omp simd     
  for (int cell_gid = 0; cell_gid < mesh.num_cells(); cell_gid++){

    cell_state.ie( 0, cell_gid ) = cell_state.ie( 1, cell_gid );
    cell_state.total_energy( 0, cell_gid ) = cell_state.total_energy( 1, cell_gid );
    
    for (int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_cell(); gauss_lid++){
      int gauss_gid = mesh.gauss_in_cell(cell_gid, gauss_lid);
      mat_pt.sie(0,gauss_gid) = mat_pt.sie(1,gauss_gid);
    }
  }

  get_control_coeffs(); 

/*
 // update control coeffs //
  for (int t_step = 0; t_step < correction_storage; t_step++){
    for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
      for (int basis = 0; basis < ref_elem.num_basis(); basis++){
        for (int dim = 0; dim < mesh.num_dim(); dim++){
          elem_state.vel_coeffs( 0, elem_gid, basis, dim ) = elem_state.vel_coeffs( correction_storage-1, elem_gid, basis, dim );
          elem_state.pos_coeffs( 0, elem_gid, basis, dim ) = elem_state.vel_coeffs( correction_storage-1, elem_gid, basis, dim );
        }// end loop over dim
      }// end loop over basis
      for (int t_basis = 0; t_basis < ref_elem.num_dual_basis(); t_basis++){
        elem_state.sie_coeffs( 0, elem_gid, t_basis ) = elem_state.sie_coeffs( correction_storage-1, elem_gid, t_basis);
      }
    }// end loop over elem_gid
  }
*/
}// end update


