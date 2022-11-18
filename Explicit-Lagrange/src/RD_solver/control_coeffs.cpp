#include<iostream>

#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"

using namespace utils;

void get_control_coeffs(){

   /////////////// ----------------------------/ ///////////////////
   //------------/////--control only at vertices--/////-----------//
   ////////////////////////////////////////////////////////////////

  for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){

    // initialize control coeffs to zero //
#pragma omp simd     

   for (int t_step = 0; t_step < num_correction_steps+1; t_step++){
      for (int j = 0; j < ref_elem.num_basis(); j++){
        for( int k = 0; k < mesh.num_dim(); k ++){
          elem_state.BV_vel_coeffs(t_step, elem_gid, j, k) = 0.0;
        }// end loop over k
      }// end loop over j
    }// end loop over t_step


#pragma omp simd      
    for (int dim = 0; dim < mesh.num_dim(); dim++){
      for (int basis_id = 0; basis_id < ref_elem.num_basis(); basis_id++){
      	for (int vertex = 0; vertex < ref_elem.num_basis(); vertex++){
	  int node_lid = elem.vert_node_map(vertex);
    	  int node_gid = mesh.nodes_in_elem( elem_gid, node_lid );
	  elem_state.BV_vel_coeffs(0, elem_gid, basis_id, dim) += elem_state.BV_mat_inv( basis_id, vertex ) * node.vel( 0, node_gid, dim );
	}// end loop over vertex
      }// end loop over basis
    }// end loop over dim 

  }// end loop over elem_gid

/*
  //// print statements ///
  for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
  std::cout << " ------- elem id ------- " << std::endl;
  std::cout << elem_gid << std::endl;
  std::cout << " ----------------------- " << std::endl;
    for (int dim = 0; dim < mesh.num_dim(); dim++){
       std::cout << "-------- dim ------" <<std::endl;
       std::cout << dim << std::endl;
      for (int basis_id = 0; basis_id < ref_elem.num_basis(); basis_id++){
        std::cout << elem_state.BV_vel_coeffs( 0, elem_gid, basis_id, dim ) << ", ";  
      }
      std::cout<<std::endl;
    }
    std::cout << " ----------------------- " << std::endl;
    std::cout << " ----------------------- " << std::endl;
  }
*/


}// end get_control_coeffs()


/*
 

/////////////////////////////////////////////////////////////////////////////////////////////
///////////////////// -------- pseudo-inverse --------- /////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

  for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){

    // initialize control coeffs to zero //
#pragma omp simd      
    for( int k = 0; k < mesh.num_dim(); k ++){
      for (int j = 0; j < ref_elem.num_basis(); j++){
        elem_state.BV_vel_coeffs(0, elem_gid, j, k) = 0.0;
	elem_state.BV_pos_coeffs(0, elem_gid, j, k) = 0.0;   
      }// end loop over j
    }// end loop over k

#pragma omp simd      
    for (int dim = 0; dim < mesh.num_dim(); dim++){
      for (int basis_id = 0; basis_id < ref_elem.num_basis(); basis_id++){
      	for (int node_lid = 0; node_lid < mesh.num_nodes_in_elem(); node_lid++){
	  int node_gid = mesh.nodes_in_elem( elem_gid, node_lid );
    	  elem_state.BV_vel_coeffs(0, elem_gid, basis_id, dim) += elem_state.BV_mat_inv( basis_id, node_lid ) * node.vel( 0, node_gid, dim);
          elem_state.BV_pos_coeffs(0, elem_gid, basis_id, dim) += elem_state.BV_mat_inv( basis_id, node_lid ) * node.coords( 0, node_gid, dim); 
        }// end loop over node_lid
      }// end loop over basis
    }// end loop over dim 
    
    // push values to all time storage //

#pragma omp simd      
    for (int dim = 0; dim < mesh.num_dim(); dim++){
      for (int basis = 0; basis < ref_elem.num_basis(); basis++){
        for (int t_step = 1; t_step <= num_correction_steps; t_step++){
          elem_state.BV_vel_coeffs(t_step, elem_gid, basis, dim) = elem_state.BV_vel_coeffs( 0, elem_gid, basis, dim);
	  elem_state.BV_pos_coeffs(t_step, elem_gid, basis, dim) = elem_state.BV_pos_coeffs( 0, elem_gid, basis, dim);
        }// end loop over t_step
      }// end loop over basis
    }// end loop over dim

  }// end loop over elem_gid

*/

/*


*/

/*
// push values to all time storage //
#pragma omp simd      
    for (int dim = 0; dim < mesh.num_dim(); dim++){
      for (int basis = 0; basis < ref_elem.num_basis(); basis++){
        for (int t_step = 1; t_step <= num_correction_steps; t_step++){
          elem_state.BV_vel_coeffs(t_step, elem_gid, basis, dim) = elem_state.BV_vel_coeffs( 0, elem_gid, basis, dim);
	  elem_state.BV_pos_coeffs(t_step, elem_gid, basis, dim) = elem_state.BV_pos_coeffs( 0, elem_gid, basis, dim);
        }// end loop over t_step
      }// end loop over basis
    }// end loop over dim
*/


