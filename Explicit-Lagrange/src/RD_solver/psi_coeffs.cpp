#include<iostream>
#include<math.h>
#include<algorithm>

#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"

using namespace utils;

void get_betas(){

    for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
      for (int vertex = 0; vertex < ref_elem.num_basis(); vertex++){
        for (int dim = 0; dim < mesh.num_dim(); dim++){
          elem_state.psi_coeffs( elem_gid, vertex, dim ) = 0.0; 
        }// end loop over dim
      }// end loop over vertex
    }// end loop over elem_gid

    for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
      for (int vertex = 0; vertex < ref_elem.num_basis(); vertex++){
        for (int dim = 0; dim < mesh.num_dim(); dim++){
      	  real_t num = 0.0;
	  num = std::max( 2.5e-30, elem_state.nodal_res(elem_gid, vertex, dim)/elem_state.total_res(elem_gid, dim) );
	  real_t denom = 0.0;
	  for (int k = 0; k < ref_elem.num_basis(); k++){
	    denom += std::max(2.5e-30, elem_state.nodal_res(elem_gid, k, dim)/elem_state.total_res(elem_gid,dim) );
	  }// end loop over k
	 
	  elem_state.psi_coeffs( elem_gid, vertex, dim ) = num/denom;
	  
	  //std::cout << elem_state.psi_coeffs( elem_gid, vertex, dim) << std::endl;
        }// end loop over dim
      }// end loop over vertex
    }// end loop over elem_gid

    /*
    // check that coeffs sum to 1 //

    for (int dim = 0; dim < mesh.num_dim(); dim++){
      for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
	real_t sum = 0.0;
        for (int vertex = 0; vertex < ref_elem.num_basis(); vertex++){
          sum += elem_state.psi_coeffs( elem_gid, vertex, dim);
        }// end loop over vertex
	std::cout << " sum of betas in elem "<< elem_gid << " in dim " << dim<< " is " << sum << std::endl;
      }// end loop over elem_gid
    }// end loop over dim
    */
}// end get betas
