
#include<iostream>
#include<math.h>

#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "slam.h"
#include "variables.h"

using namespace utils;


void BV_inv(){

 //-------------------/////////////////////////////////////////////-------------------------
 //----///////////------- control points only at vertices ----------///////////////---------
 //-----------------/////////////////////////////////////////////---------------------------


/// --- get control pooints for kinematic DOFs --- ///	
    real_t B_a[ref_elem.num_basis() * ref_elem.num_basis()];
    real_t check_B_a[ref_elem.num_basis() * ref_elem.num_basis()];
  
    auto B = ViewCArray <real_t> ( &B_a[0], ref_elem.num_basis(), ref_elem.num_basis() );
    auto check_B = ViewCArray <real_t> ( &check_B_a[0], ref_elem.num_basis(), ref_elem.num_basis() );


#pragma omp simd
    
    for (int j = 0; j < ref_elem.num_basis(); j++){
      for (int i = 0; i < ref_elem.num_basis(); i++){
        B(i, j) = 0.0;
        elem_state.BV_mat_inv( i, j ) = 0.0;
        check_B(i, j) = 0.0;
      }// end loop over i
    }// end loop over j  
  
#pragma omp simd
    for ( int index = 0; index < ref_elem.num_basis(); index++){
      for (int vertex = 0; vertex < ref_elem.num_basis(); vertex++){
	int node_lid = ref_elem.vert_node_map(vertex);
        B( vertex, index) = ref_elem.ref_nodal_basis(node_lid, index );
        check_B( vertex, index) = B( vertex, index );
      }
    }
    
/*
      std::cout << "---- B transpose ----" << std::endl;  
      for (int index = 0; index < ref_elem.num_basis(); index++){
        for (int vertex = 0; vertex < ref_elem.num_basis(); vertex++){
      	 std::cout << B( vertex, index)<< ", ";
      }
      std::cout<<std::endl;
    } 
*/
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

  
    auto B_lu = ViewCArray <real_t> (&B(0,0), ref_elem.num_basis(), ref_elem.num_basis());
   
    LU_decompos(B_lu, lu_index, parity, ref_elem.num_basis());
    
    auto B_inv = ViewCArray <real_t> ( &elem_state.BV_mat_inv( 0, 0 ),  ref_elem.num_basis(), ref_elem.num_basis());
    LU_invert( B_lu, lu_index, B_inv, col, ref_elem.num_basis() );
    
/*  
    real_t val1 = 0.0;
    real_t val2 = 0.0;
    for (int j = 0; j < ref_elem.num_basis(); j++ ){
      for (int i = 0; i <  ref_elem.num_basis(); i++ ){
        for (int k = 0; k < ref_elem.num_basis(); k++ ){ 
          val1 += check_B(j,k)*elem_state.BV_mat_inv( k,i);
          val2 += elem_state.BV_mat_inv(j,k)*check_B(k,i);
        }// end loop over k
        std::cout << " B*Binv with i =  "<< i << " and j ="<< j << " is = " << val1 << std::endl;       
        std::cout << " Binv*B with i =  "<< i << " and j ="<< j << " is = " << val2 << std::endl;
        val1 = 0.0;
        val2 = 0.0;
      }// end loop over i
    }// end loop over j  
*/

/// --- get control points for thermodynamic DOFs --- ///
/*	
    real_t dual_B_a[ref_elem.num_dual_basis() * ref_elem.num_dual_basis()];
    real_t check_dual_B_a[ref_elem.num_dual_basis() * ref_elem.num_dual_basis()];
  
    auto dual_B = ViewCArray <real_t> ( &dual_B_a[0], ref_elem.num_dual_basis(), ref_elem.num_dual_basis() );
    auto check_dual_B = ViewCArray <real_t> ( &check_dual_B_a[0], ref_elem.num_dual_basis(), ref_elem.num_dual_basis() );


#pragma omp simd
    
    for (int j = 0; j < ref_elem.num_dual_basis(); j++){
      for (int i = 0; i < ref_elem.num_dual_basis(); i++){
        dual_B(i, j) = 0.0;
        elem_state.dual_BV_mat_inv( i, j ) = 0.0;
        check_dual_B(i, j) = 0.0;
      }// end loop over i
    }// end loop over j  
  
#pragma omp simd
    for ( int index = 0; index < ref_elem.num_dual_basis(); index++){
      for (int vertex = 0; vertex < ref_elem.num_dual_basis(); vertex++){
	int node_lid = ref_elem.vert_node_map(vertex); 
        dual_B( vertex, index) = ref_elem.ref_nodal_dual_basis(node_lid, index );
        check_dual_B( vertex, index) = B( vertex, index );
      }
    }
    

      std::cout << "---- dual B transpose ----" << std::endl;  
      for (int index = 0; index < ref_elem.num_dual_basis(); index++){
        for (int vertex = 0; vertex < ref_elem.num_dual_basis(); vertex++){
      	 std::cout << dual_B( vertex, index)<< ", ";
        }
        std::cout<<std::endl;
      } 


    int dual_lu_index_a[ref_elem.num_dual_basis()];
    auto dual_lu_index = ViewCArray <int> (&dual_lu_index_a[0], ref_elem.num_dual_basis());
    for(int i=0; i < ref_elem.num_dual_basis(); i++){
      dual_lu_index(i) = 0;
    };
    int dual_parity = 0;

    real_t dual_col_a[ref_elem.num_dual_basis()];
    auto dual_col = ViewCArray <real_t> (&dual_col_a[0], ref_elem.num_dual_basis());
    for(int i=0; i < ref_elem.num_dual_basis(); i++){
      dual_col(i) = 0;
    };

  
    auto dual_B_lu = ViewCArray <real_t> (&dual_B(0,0), ref_elem.num_dual_basis(), ref_elem.num_dual_basis());
   
    LU_decompos(dual_B_lu, dual_lu_index, dual_parity, ref_elem.num_dual_basis());
    
    auto dual_B_inv = ViewCArray <real_t> ( &elem_state.dual_BV_mat_inv( 0, 0 ),  ref_elem.num_dual_basis(), ref_elem.num_dual_basis());
    LU_invert( dual_B_lu, dual_lu_index, dual_B_inv, dual_col, ref_elem.num_dual_basis() );
    
*/

           ///////////////////////////////////////////////////////////////
          //------------------------------------------------------------//
          ////// ----- pseudo-inverse for thermodynamic DOFs ----- //////
          //----------------------------------------------------------//
         //////////////////////////////////////////////////////////////
	 
    real_t dual_B_a[ref_elem.num_dual_basis() * mesh.num_nodes_in_elem()];
    real_t check_dual_B_a[ref_elem.num_basis() * mesh.num_nodes_in_elem()];
  
    auto dual_B = ViewCArray <real_t> (&dual_B_a[0], mesh.num_nodes_in_elem(), ref_elem.num_dual_basis());
    auto check_dual_B = ViewCArray <real_t> (&check_dual_B_a[0], mesh.num_nodes_in_elem(), ref_elem.num_dual_basis());


#pragma omp simd
    
    for (int j = 0; j < ref_elem.num_dual_basis(); j++){
      for (int i = 0; i < mesh.num_nodes_in_elem(); i++){
        dual_B(i, j) = 0.0; // in R^{num_nodes x num_basis}
        elem_state.dual_BV_mat_inv( j, i ) = 0.0; // in R^{num_dual_basis x num_nodes}
        check_dual_B(i, j) = 0.0;
      }// end loop over i
    }// end loop over j  
  
#pragma omp simd
    for ( int node_lid = 0; node_lid < mesh.num_nodes_in_elem(); node_lid++){
      for (int k = 0; k < ref_elem.num_dual_basis(); k++){
        dual_B(node_lid, k) = ref_elem.ref_nodal_dual_basis(node_lid, k );
        check_dual_B(node_lid, k) = dual_B(node_lid, k);
      }
    }
    
     /* 
     std::cout << "---- dual B transpose ----" << std::endl;  
      for (int index = 0; index < ref_elem.num_dual_basis(); index++){
        for (int node_lid = 0; node_lid < mesh.num_nodes_in_elem(); node_lid++){
      	std::cout << dual_B(node_lid, index)<< ", ";
	std::cout << " check at i = " << node_lid << " j = " << index << " is " << check_dual_B(node_lid, index) << std::endl;
        }
        std::cout<<std::endl;
      } 
      */

    // Create B^T //
    real_t dual_BT_a[ref_elem.num_dual_basis() * mesh.num_nodes_in_elem()];
    auto dual_BT = ViewCArray <real_t> (&dual_BT_a[0], ref_elem.num_dual_basis(), mesh.num_nodes_in_elem() );
    for (int j = 0; j < mesh.num_nodes_in_elem(); j++){
      for (int i = 0; i < ref_elem.num_dual_basis(); i++){
        dual_BT(i,j) = 0.0;
      }
    }
    
    for (int j = 0; j < ref_elem.num_dual_basis(); j++){
      for (int i = 0; i < mesh.num_nodes_in_elem(); i++){
        dual_BT(j,i) = dual_B(i,j);
      }
    }
    

    // Create B^T.B //
    real_t dual_BTB_a[ref_elem.num_dual_basis() * ref_elem.num_dual_basis()];
    auto dual_BTB = ViewCArray <real_t> (&dual_BTB_a[0], ref_elem.num_dual_basis(), ref_elem.num_dual_basis() );

#pragma omp simd
    for (int i = 0; i < ref_elem.num_dual_basis(); i++){
      for(int j = 0; j < ref_elem.num_dual_basis(); j++){
        dual_BTB(j,i) = 0.0;// in R^{num_dual_basis x num_dual_basis}
      }// end loop over j
    }//end loop over i

#pragma omp simd
  //  std::cout << "---- dual_ B^T.B ---- "<< std::endl;
    for (int j =0; j<ref_elem.num_dual_basis(); j++){
      for (int i = 0; i < ref_elem.num_dual_basis(); i++){
        for (int k = 0; k < mesh.num_nodes_in_elem(); k++){
	   dual_BTB(i,j) += dual_BT(i,k) * dual_B(k, j);
	}// end loop over k
//	std::cout << dual_BTB(i,j) << ", ";
      }// end loop over i
    //  std::cout << std::endl;
    }// end loop over j

    int dual_lu_index_a[ref_elem.num_dual_basis()];
    auto dual_lu_index = ViewCArray <int> (&dual_lu_index_a[0], ref_elem.num_dual_basis());
    for(int i=0; i < ref_elem.num_dual_basis(); i++){
      dual_lu_index(i) = 0;
    };
    int dual_parity = 0;

    real_t dual_col_a[ref_elem.num_dual_basis()];
    auto dual_col = ViewCArray <real_t> (&dual_col_a[0], ref_elem.num_dual_basis());
    for(int i=0; i < ref_elem.num_dual_basis(); i++){
      dual_col(i) = 0;
    };

  
    auto dual_BTB_lu = ViewCArray <real_t> (&dual_BTB(0,0), ref_elem.num_dual_basis(), ref_elem.num_dual_basis());
   
    LU_decompos(dual_BTB_lu, dual_lu_index, dual_parity, ref_elem.num_dual_basis());
    
    // get inverse of BTB //
    real_t dual_BTB_inv_a[ref_elem.num_dual_basis() * ref_elem.num_dual_basis() ];
    auto dual_BTB_inv = ViewCArray <real_t> (&dual_BTB_inv_a[0],  ref_elem.num_dual_basis(), ref_elem.num_dual_basis());
    LU_invert( dual_BTB, dual_lu_index, dual_BTB_inv, dual_col, ref_elem.num_dual_basis() );
    

#pragma omp simd
    // get B^+ = (B^T.B)^(-1).B^T
//    std::cout << "---- pseudo inverse ----" << std::endl;
    for (int dim_i = 0; dim_i < ref_elem.num_dual_basis(); dim_i++){
      for (int dim_j = 0; dim_j < mesh.num_nodes_in_elem(); dim_j++){
        for (int k = 0; k < ref_elem.num_dual_basis(); k++){
	  elem_state.dual_BV_mat_inv( dim_i, dim_j ) +=  dual_BTB_inv( dim_i, k ) * dual_BT( k, dim_j );
        }// end loop over k
      }// end loop over dim_i
    }// end loop over dim_j
/*
    for (int dim_i = 0; dim_i < ref_elem.num_dual_basis(); dim_i++){
      for (int dim_j = 0; dim_j < mesh.num_nodes_in_elem(); dim_j++){
        std::cout << elem_state.BV_mat_inv( dim_i, dim_j ) <<", ";
      }
    } 
    std::cout << std::endl;
*/

    //real_t val1 = 0.0;
    //real_t val2 = 0.0;
/*    
    for (int j = 0; j < mesh.num_nodes_in_elem(); j++ ){
      for (int i = 0; i <  mesh.num_nodes_in_elem(); i++ ){
        for (int k = 0; k < ref_elem.num_dual_basis(); k++ ){ 
          val1 += check_dual_B(j,k)*elem_state.dual_BV_mat_inv( k,i);
        }// end loop over k
        std::cout << " dual_B*dual_Binv with i =  "<< i << " and j ="<< j << " is = " << val1 << std::endl;       
        val1 = 0.0;
      }// end loop over i
    }// end loop over j  
    std::cout << std::endl;
*/
/*
    for (int j = 0; j < ref_elem.num_dual_basis(); j++ ){
      for (int i = 0; i <  ref_elem.num_dual_basis(); i++ ){
        for (int k = 0; k < mesh.num_nodes_in_elem(); k++ ){ 
          val2 += elem_state.dual_BV_mat_inv(j,k)*check_dual_B(k,i);
        }// end loop over k
        std::cout << " dual_Binv*dual_B with i =  "<< i << " and j ="<< j << " is = " << val2 << std::endl;
        val2 = 0.0;
      }// end loop over i
    }// end loop over j  
    std::cout << std::endl;
*/
}// end BV inverse


/*



*/



/*


*/
