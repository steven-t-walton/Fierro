#include<iostream>

#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"

using namespace utils;

void get_control_coeffs(){

  real_t BV_coeffs_a[mesh.num_elems()*elem.num_basis()*mesh.num_dim()];
  auto BV_coeffs = ViewCArray <real_t> (&BV_coeffs_a[0], mesh.num_elems(), elem.num_basis(), mesh.num_dim());
#pragma omp simd
  // initialize control coeffs to zero //
  for( int k = 0; k < mesh.num_dim(); k ++){
    for (int j = 0; j < elem.num_basis(); j++){
      for (int i = 0; i < mesh.num_elems(); i++){
        BV_coeffs(i,j,k) =  0.0;
      }// end loop over i
    }// end loop over j
  }// end loop over k

#pragma omp simd      
  for (int elem_gid = 0 ; elem_gid < mesh.num_elems(); elem_gid++){
    for (int dim = 0; dim < mesh.num_dim(); dim++){
      for(int basis_id = 0; basis_id < elem.num_basis(); basis_id++){

        for (int k = 0; k < elem.num_basis(); k++){
          int node_basis_id = elem.vert_node_map(k);
          int interp_gid = mesh.gauss_in_elem(elem_gid, node_basis_id);
          BV_coeffs(elem_gid,basis_id,dim) += elem_state.BV_mat_inv(basis_id, k)*node.vel(0,interp_gid ,dim);
        }// end loop over k

        //std::cout << "Control coeffs in elem  " << elem_gid << " and dim " << dim << " are " << BV_coeffs(elem_gid, basis_id, dim) << std::endl;
        int node_id = elem.vert_node_map(basis_id);
        int node_gid_for_control_coeff = mesh.nodes_in_elem(elem_gid, node_id);
        
        auto vel = ViewCArray <real_t> (&node.vel(0, node_gid_for_control_coeff, 0), mesh.num_dim() );
       
//        std::cout << " vel before coeff assignment " << vel(dim) << std::endl;
//        std::cout << " " << std::endl;
        vel( dim ) = BV_coeffs(elem_gid, basis_id, dim);
//        std::cout << " vel from BV_coeffs is " << vel(dim) << " at dim " << dim << std::endl;
	for (int t_step = 1; t_step < num_correction_steps; t_step++){
	  node.vel(t_step, node_gid_for_control_coeff, dim) = vel( dim );
	}
        
      }// end loop over basis id
    }// end loop over dim
  }// end loop over elem_gid
}
