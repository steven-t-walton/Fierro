
#include<iostream>
#include<math.h>

#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "slam.h"
#include "variables.h"

using namespace utils;


void BV_inv(){


    real_t B_a[ref_elem.num_basis() * mesh.num_nodes_in_elem()];
    //real_t check_B_a[ref_elem.num_basis() * mesh.num_nodes_in_elem()];
  
    auto B = ViewCArray <real_t> (&B_a[0], mesh.num_nodes_in_elem(), ref_elem.num_basis());
    //auto check_B = ViewCArray <real_t> (&check_B_a[0], mesh.num_nodes_in_elem(), ref_elem.num_basis());


#pragma omp simd
    
    for (int j = 0; j < ref_elem.num_basis(); j++){
      for (int i = 0; i < mesh.num_nodes_in_elem(); i++){
        B(i, j) = 0.0; // in R^{num_nodes x num_basis}
        elem_state.BV_mat_inv( j, i ) = 0.0; // in R^{num_basis x num_nodes}
       // check_B(i, j) = 0.0;
      }// end loop over i
    }// end loop over j  

#pragma omp simd
    for ( int node_lid = 0; node_lid < mesh.num_nodes_in_elem(); node_lid++){
      for (int k = 0; k < ref_elem.num_basis(); k++){
        B(node_lid, k) = ref_elem.ref_nodal_basis(node_lid, k );
       // check_B(node_lid, k) = B(node_lid, k);
      }
    }
    

     // std::cout << "---- B transpose ----" << std::endl;  
      for (int index = 0; index < ref_elem.num_basis(); index++){
        for (int node_lid = 0; node_lid < mesh.num_nodes_in_elem(); node_lid++){
      //	 std::cout << B(node_lid, index)<< ", ";
	//std::cout << " check at i = " << vertex << " j = " << index << " is " << check_B(vertex, index) << std::endl;
      }
    //  std::cout<<std::endl;
    } 

    // Create B^T //
    real_t BT_a[ref_elem.num_basis() * mesh.num_nodes_in_elem()];
    auto BT = ViewCArray <real_t> (&BT_a[0], ref_elem.num_basis(), mesh.num_nodes_in_elem() );
    for (int j = 0; j < mesh.num_nodes_in_elem(); j++){
      for (int i = 0; i < ref_elem.num_basis(); i++){
        BT(i,j) = 0.0;
      }
    }
    
    for (int j = 0; j < ref_elem.num_basis(); j++){
      for (int i = 0; i < mesh.num_nodes_in_elem(); i++){
        BT(j,i) = B(i,j);
      }
    }
    

    // Create B^T.B //
    real_t BTB_a[ref_elem.num_basis() * ref_elem.num_basis()];
    auto BTB = ViewCArray <real_t> (&BTB_a[0], ref_elem.num_basis(), ref_elem.num_basis() );

#pragma omp simd
    for (int i = 0; i < ref_elem.num_basis(); i++){
      for(int j = 0; j < ref_elem.num_basis(); j++){
        BTB(j,i) = 0.0;// in R^{num_basis x num_basis}
      }// end loop over j
    }//end loop over i

#pragma omp simd
  //  std::cout << "---- B^T.B ---- "<< std::endl;
    for (int j =0; j<ref_elem.num_basis(); j++){
      for (int i = 0; i < ref_elem.num_basis(); i++){
        for (int k = 0; k < mesh.num_nodes_in_elem(); k++){
	   BTB(i,j) += BT(i,k) * B(k, j);
	}// end loop over k
//	std::cout << BTB(i,j) << ", ";
      }// end loop over i
     // std::cout << std::endl;
    }// end loop over j

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

  
    auto BTB_lu = ViewCArray(&BTB(0,0), ref_elem.num_basis(), ref_elem.num_basis());
   
    LU_decompos(BTB_lu, lu_index, parity, ref_elem.num_basis());
    
    // get inverse of BTB //
    real_t BTB_inv_a[ref_elem.num_basis() * ref_elem.num_basis() ];
    auto BTB_inv = ViewCArray <real_t> (&BTB_inv_a[0],  ref_elem.num_basis(), ref_elem.num_basis());
    LU_invert( BTB, lu_index, BTB_inv, col, ref_elem.num_basis() );
    

#pragma omp simd
    // get B^+ = (B^T.B)^(-1).B^T
    //std::cout << "---- pseudo inverse ----" << std::endl;
    for (int dim_j = 0; dim_j < mesh.num_nodes_in_elem(); dim_j++){
      for (int dim_i = 0; dim_i < ref_elem.num_basis(); dim_i++){
        for (int k = 0; k < ref_elem.num_basis(); k++){
	  elem_state.BV_mat_inv( dim_i, dim_j ) +=  BTB_inv( dim_i, k ) * BT( k, dim_j );
        }// end loop over k
      //  std::cout << elem_state.BV_mat_inv( dim_i, dim_j ) <<", ";
      }// end loop over dim_i
    //  std::cout << std::endl;
    }// end loop over dim_j
  
  
/*  
  //real_t val1 = 0.0;
  real_t val2 = 0.0;
    for (int i = 0; i < ref_elem.num_basis(); i++){
      for (int j = 0; j <  ref_elem.num_basis(); j++){
        for (int k = 0; k < mesh.num_nodes_in_elem(); k++){
          //val1 += check_B(j,k)*elem_state.BV_mat_inv( k,i);
          val2 += elem_state.BV_mat_inv(j,k)*check_B(k,i);
        }// end loop over k
        //std::cout << " B*Binv with j =  "<< j << " and i ="<< i << " is = " << val1 << std::endl;       
        std::cout << " Binv*B with j =  "<< j << " and i ="<< i << " is = " << val2 << std::endl;
        //val1 = 0.0;
        val2 = 0.0;
      }// end loop over i
    }// end loop over j
*/  


}// end BV inverse
