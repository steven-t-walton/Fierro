#include<iostream>

#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"

using namespace utils;

void get_control_coeffs(){

  real_t BV_coeffs_a[ ref_elem.num_basis()*mesh.num_dim() ];//mesh.num_elems()*elem.num_basis()*mesh.num_dim()];
  auto BV_coeffs = ViewCArray <real_t> (&BV_coeffs_a[0], ref_elem.num_basis(), mesh.num_dim());//mesh.num_elems(), elem.num_basis(), mesh.num_dim());

  int rid_a[ref_elem.num_basis()];
  auto rid = ViewCArray <int> (&rid_a[0], ref_elem.num_basis());
  for (int m=0; m < ref_elem.num_basis(); m++) rid(m) = 0;


#pragma omp simd      
  for (int elem_gid = 0 ; elem_gid < mesh.num_elems(); elem_gid++){

//#pragma omp simd
    int ind = 0;
    for (int k = 0; k <  cbrt(ref_elem.num_basis()); k++){
      for (int j = 0; j < cbrt(ref_elem.num_basis()); j++){
        for (int i = 0; i < cbrt(ref_elem.num_basis()); i++){
          rid(ind) = ref_elem.node_rid(i,j,k);
          //std::cout << "node_rid at i = " << i << " j = " << j << " k = " << k << " is " << rid(ind) << std::endl;
          ind++;
        };// end loop over i
      };// end loop over j
    };// end loop over k

//#pragma omp simd
    // initialize control coeffs to zero //
    for( int k = 0; k < mesh.num_dim(); k ++){
      for (int j = 0; j < ref_elem.num_basis(); j++){
        //for (int i = 0; i < mesh.num_elems(); i++){
          BV_coeffs(j,k) = 0.0; //BV_coeffs(i,j,k) =  0.0;
        //}// end loop over i
      }// end loop over j
    }// end loop over k

    for(int basis_id = 0; basis_id < ref_elem.num_basis(); basis_id++){
      for (int dim = 0; dim < mesh.num_dim(); dim++){

        for (int k = 0; k < ref_elem.num_basis(); k++){
          int node_basis_id = elem.vert_node_map( k );
          int interp_gid = mesh.nodes_in_elem(elem_gid, node_basis_id);
          BV_coeffs(basis_id, dim) = elem_state.BV_mat_inv( basis_id, k ) * node.vel( 0, interp_gid, dim); //elem_gid, basis_id, dim ) += elem_state.BV_mat_inv( basis_id, k ) * node.vel( 0, interp_gid, dim);
        }// end loop over k

//        std::cout << "Control coeffs in elem  " << elem_gid << " and dim " << dim << "for basis_id " << basis_id << " are " << BV_coeffs(basis_id, dim) << std::endl;
//	std::cout << std::endl;

	int node_id = elem.vert_node_map( basis_id );
        int node_gid_for_control_coeff = mesh.nodes_in_elem( elem_gid, node_id );
        
        //auto vel = ViewCArray <real_t> (&node.vel(0, node_gid_for_control_coeff, 0), mesh.num_dim() );
       /* 
	std::cout << std::endl;
        std::cout << " vel before coeff assignment " << node.vel(0, node_gid_for_control_coeff, dim) << std::endl;
      */
	node.vel(0, node_gid_for_control_coeff, dim ) = BV_coeffs( basis_id, dim );
      /*
	std::cout << " vel from BV_coeffs is " << node.vel(0, node_gid_for_control_coeff, dim ) << " at dim " << dim << std::endl;
        std::cout << std::endl;
      */ 	
	for (int t_step = 1; t_step < num_correction_steps; t_step++){
	  node.vel(t_step, node_gid_for_control_coeff, dim) = node.vel( 0, node_gid_for_control_coeff, dim );
	}
        
      }// end loop over dim
    }// end loop over basis_id
  }// end loop over elem_gid
}
