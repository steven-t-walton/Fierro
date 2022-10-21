

#include<iostream>
#include<math.h>

#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "slam.h"
#include "variables.h"

#include "bernstein_polynomials.h"

using namespace utils;


void BV_inv(){

    real_t B_a[ref_elem.num_basis()*ref_elem.num_basis()];
    //real_t check_B_a[ref_elem.num_basis()*ref_elem.num_basis()];
  
    auto B = ViewCArray <real_t> (&B_a[0], ref_elem.num_basis(), ref_elem.num_basis());
    //auto check_B = ViewCArray <real_t> (&check_B_a[0], ref_elem.num_basis(), ref_elem.num_basis());


#pragma omp simd
    for (int j = 0; j < ref_elem.num_basis(); j++){
      for (int i = 0; i < ref_elem.num_basis(); i++){
        B(i, j) = 0.0;
        elem_state.BV_mat_inv( i, j ) = 0.0;
      }
    }  

    for (int basis_id = 0; basis_id < ref_elem.num_basis(); basis_id++){
      for (int index = 0; index < ref_elem.num_basis(); index++){  
        B(index,basis_id) = mesh.gauss_pt_det_j() * ref_elem.ref_nodal_basis(index, basis_id );
      //check_B(index,basis_id) = B(index,basis_id);
      }
    }

    int lu_index_a[ref_elem.num_basis()];
    auto lu_index = ViewCArray <int> (&lu_index_a[0], ref_elem.num_basis());
    for(int i=0; i < ref_elem.num_basis(); i++){
      lu_index(i) = 0;
    };
    int parity = 0;

    real_t col_a[ref_elem.num_basis()];
    auto col = ViewCArray <real_t> (&col_a[0], ref_elem.num_basis());
    for(int i=0; i < ref_elem.num_basis(); i++){
      col(i) = 0;
    };

  
    auto B_lu = ViewCArray(&B(0,0), ref_elem.num_basis(), ref_elem.num_basis());
   
    LU_decompos(B_lu, lu_index, parity, ref_elem.num_basis());
    
    auto B_inv = ViewCArray <real_t> ( &elem_state.BV_mat_inv( 0, 0 ), ref_elem.num_basis(), ref_elem.num_basis());
    LU_invert( B, lu_index, B_inv, col, ref_elem.num_basis() );

}// end BV inverse







/*
*/

  /*  
  real_t val1 = 0.0;
  real_t val2 = 0.0;
  
  for (int i = 0; i < ref_elem.num_basis(); i++){
    for (int j = 0; j < ref_elem.num_basis(); j++){
      for (int k = 0; k < ref_elem.num_basis(); k++){
        val1 += check_B(j,k)*elem_state.BV_mat_inv(k,i);
        val2 += elem_state.BV_mat_inv(j,k)*check_B(k,i);
      }// end loop over k
      std::cout << " B*Binv with j =  "<< j << " and i ="<< i << " is = " << val1 << std::endl;       
      std::cout << " Binv*B with j =  "<< j << " and i ="<< i << " is = " << val2 << std::endl;
      val1 = 0.0;
      val2 = 0.0;
    }// end loop over i
  }// end loop over j
  */

