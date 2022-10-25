#include<iostream>

#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"

using namespace utils;

void get_control_coeffs(){
  
  for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
    // Views of coeffs // 
    auto BV_vel_coeffs = ViewCArray <real_t> (&elem_state.BV_vel_coeffs( 0, elem_gid, 0, 0 ), ref_elem.num_basis(), mesh.num_dim() );
    auto BV_pos_coeffs = ViewCArray <real_t> (&elem_state.BV_pos_coeffs( 0, elem_gid, 0, 0 ), ref_elem.num_basis(), mesh.num_dim() );

    // View of BV_inv //
    auto inverse = ViewCArray <real_t> (&elem_state.BV_mat_inv( elem_gid, 0, 0 ), ref_elem.num_basis(), mesh.num_nodes_in_elem() );
    
    // initialize control coeffs to zero //
#pragma omp simd      
    for( int k = 0; k < mesh.num_dim(); k ++){
      for (int j = 0; j < ref_elem.num_basis(); j++){
        BV_vel_coeffs(j,k) = 0.0;
	BV_pos_coeffs(j,k) = 0.0;   
      }// end loop over j
    }// end loop over k

#pragma omp simd      
    for (int dim = 0; dim < mesh.num_dim(); dim++){
      for (int basis_id = 0; basis_id < ref_elem.num_basis(); basis_id++){
      	for (int node_lid = 0; node_lid < mesh.num_nodes_in_elem(); node_lid++){
          int node_gid = mesh.nodes_in_elem( elem_gid, node_lid );
    	  BV_vel_coeffs(basis_id, dim) += inverse( basis_id, node_lid ) * node.vel( 0, node_gid, dim); 
          BV_pos_coeffs(basis_id, dim) += inverse( basis_id, node_lid ) * node.coords( 0, node_gid, dim); 
        }// end loop over node_lid
      }// end loop over basis
    }// end loop over dim 
    
    // push values to all time storage //

#pragma omp simd      
    for (int dim = 0; dim < mesh.num_dim(); dim++){
      for (int basis = 0; basis < ref_elem.num_basis(); basis++){
        for (int t_step = 1; t_step <= num_correction_steps; t_step++){
          elem_state.BV_vel_coeffs(t_step, elem_gid, basis, dim) = BV_vel_coeffs(basis, dim);
	  elem_state.BV_pos_coeffs(t_step, elem_gid, basis, dim) = BV_pos_coeffs(basis, dim);
        }// end loop over t_step
      }// end loop over basis
    }// end loop over dim

  }// end loop over elem_gid

}// end get_control_coeffs()


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    //real_t BV_coeffs_a[ ref_elem.num_basis()*mesh.num_dim() ];

      //for(int vert_id = 0; vert_id < ref_elem.num_basis(); vert_id++){

       // for (int k = 0; k < ref_elem.num_basis(); k++){
       //   int node_basis_id = ref_elem.vert_node_map( k );
       //   int interp_gid = mesh.nodes_in_elem(elem_gid, node_basis_id);

        //}// end loop over k

/*
	int node_id = ref_elem.vert_node_map( vert_id );
        int node_gid_for_control_coeff = mesh.nodes_in_elem( elem_gid, node_id );
	
        node.vel(0, node_gid_for_control_coeff, dim ) = BV_coeffs( vert_id, dim );
	
        for (int t_step = 1; t_step <= num_correction_steps; t_step++){
	  node.vel(t_step, node_gid_for_control_coeff, dim) = node.vel( 0, node_gid_for_control_coeff, dim );
	}
        
      }// end loop over dim
    }// end loop over basis_id
*/
