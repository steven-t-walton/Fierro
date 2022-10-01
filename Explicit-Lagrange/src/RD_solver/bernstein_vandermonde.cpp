#include<iostream>
#include<math.h>

#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "slam.h"
#include "variables.h"

#include "bernstein_polynomials.cpp"

using namespace utils;

void bernstein_vandermonde(ViewCArray <real_t> &B){
  
  real_t rid_a[ref_elem.num_basis()];
  int ind = 0;
  auto rid = ViewCArray <real_t> (rid_a, ref_elem.num_basis());
  for (int k = 0; k <  cbrt(ref_elem.num_basis()); k++){
    for (int j = 0; j < cbrt(ref_elem.num_basis()); j++){
      for (int i = 0; i < cbrt(ref_elem.num_basis()); i++){
        rid(ind) = ref_elem.node_rid(i,j,k);
        //std::cout << "node_rid at i = " << i << " j = " << j << " k = " << k << " is " << ref_elem.node_rid(i,j,k) << std::endl;
        ind++;
      }// end loop over i
    }// end loop over j
  }// end loop over k

  int degree = 0;
  for (int dim = 0; dim < mesh.num_dim(); dim++){      
    for (int index = 0; index < ref_elem.num_basis(); index++){
      for (int basis_id = 0; basis_id < ref_elem.num_basis(); basis_id++){  
        B(index,basis_id,dim) = bernstein::eval(p_order,degree, ref_elem.ref_node_positions(rid(index),dim));
        degree++;
        if ( degree == p_order) degree = 0;
        //std::cout << "B-V mat at node i = "<< index << " and j " << basis_id << " with dim "<< dim <<" is "<< B(index,basis_id,dim) << std::endl;
      }// end loop over degree
    }// end loop over index
  }// end loop over dim
      
}// end B-V matrix


void BV_inv(){//ViewCArray <real_t> &B, ViewCArray <real_t> &B_inv){
 
  real_t B_a[ref_elem.num_basis()*ref_elem.num_basis()*mesh.num_dim()];
  real_t temp_B_a[ref_elem.num_basis()*ref_elem.num_basis()*mesh.num_dim()];
  
  auto B = ViewCArray <real_t> (&B_a[0], ref_elem.num_basis(), ref_elem.num_basis(), mesh.num_dim());
  auto temp_B = ViewCArray <real_t> (&temp_B_a[0], ref_elem.num_basis(), ref_elem.num_basis(), mesh.num_dim());

  real_t rid_a[ref_elem.num_basis()];
  int ind = 0;
  auto rid = ViewCArray <real_t> (&rid_a[0], ref_elem.num_basis());
  for (int k = 0; k <  cbrt(ref_elem.num_basis()); k++){
    for (int j = 0; j < cbrt(ref_elem.num_basis()); j++){
      for (int i = 0; i < cbrt(ref_elem.num_basis()); i++){
        rid(ind) = ref_elem.node_rid(i,j,k);
        //std::cout << "node_rid at i = " << i << " j = " << j << " k = " << k << " is " << ref_elem.node_rid(i,j,k) << std::endl;
        ind++;
      }// end loop over i
    }// end loop over j
  }// end loop over k

  int degree = 0;
  for (int dim = 0; dim < mesh.num_dim(); dim++){
    for (int index = 0; index < ref_elem.num_basis(); index++){
      for (int basis_id = 0; basis_id < ref_elem.num_basis(); basis_id++){
        B(index,basis_id,dim) = bernstein::eval(p_order,degree, ref_elem.ref_node_positions(rid(index),dim));
        temp_B(index,basis_id,dim) = B(index,basis_id, dim);
        degree++;
        if ( degree == p_order+1) degree = 0;
        //std::cout << "B-V mat at node i = "<< index << " and j " << basis_id << " with dim "<< dim <<" is "<< B(index,basis_id,dim) << std::endl;
      }// end loop over degree
    }// end loop over index
  }// end loop over dim
  
  int index_a[ref_elem.num_basis()];
  auto index = ViewCArray <int> (&index_a[0], ref_elem.num_basis());
  for(int i=0; i < ref_elem.num_basis(); i++){
    index(i) = 0;
  }
  int parity = 0;
  int singular = 0;
 
  real_t col_a[ref_elem.num_basis()];
  auto col = ViewCArray <real_t> (&col_a[0], ref_elem.num_basis());
  for(int i=0; i < ref_elem.num_basis(); i++){
    col(i) = 0;
  }

  real_t r_test_a[ref_elem.num_basis()*ref_elem.num_basis()];
  real_t l_test_a[ref_elem.num_basis()*ref_elem.num_basis()];

  auto right_test_mat = ViewCArray <real_t> (&r_test_a[0], ref_elem.num_basis(), ref_elem.num_basis());
  auto left_test_mat = ViewCArray <real_t> (&l_test_a[0], ref_elem.num_basis(), ref_elem.num_basis());

  for (int dim = 0; dim < mesh.num_dim(); dim++){
    auto B_temp1 = ViewCArray(&B(0,0,dim), ref_elem.num_basis(), ref_elem.num_basis());
   
    singular = LU_decompos(B_temp1, index, parity, ref_elem.num_basis());
      
    auto B_inv_temp = ViewCArray <real_t> (&elem_state.BV_mat_inv(0,0,dim), ref_elem.num_basis(), ref_elem.num_basis());
    LU_invert(B, index, B_inv_temp, col, ref_elem.num_basis());
 
    for (int i = 0; i < ref_elem.num_basis(); i++){
      for (int j = 0; j < ref_elem.num_basis(); j++){
        left_test_mat(i,j) = 0.0;
        right_test_mat(i,j) = 0.0;
      }
    }
    
    for (int j = 0; j < ref_elem.num_basis(); j++){
      for (int i = 0; i < ref_elem.num_basis(); i++){
        for (int k = 0; k < ref_elem.num_basis(); k++){
          right_test_mat(i,j) += temp_B(i,k,dim)*elem_state.BV_mat_inv(k,j,dim);
          left_test_mat(i,j) += elem_state.BV_mat_inv(i,k,dim)*temp_B(k,j,dim);
        }// end loop over k
        std::cout << " B*Binv with i =   "<< i << " j ="<< j << " and dim =" << dim << " is " << right_test_mat(i,j) << std::endl;       
        std::cout << " Binv*B with i =   "<< i << " and j ="<< j <<  " and dim =" << dim << " is " << left_test_mat(i,j) << std::endl;
      }// end loop over i
    }// end loop over j
  }
  
}// end BV inverse
