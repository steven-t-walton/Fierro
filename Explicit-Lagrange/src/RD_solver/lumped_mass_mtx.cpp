// computes lumped mass and stores at nodes //
// Computes volume of convex hull defined by corners around a node //

#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "slam.h"
#include "variables.h"


using namespace utils;

void lumped_mass(){


  int num_basis = ref_elem.num_basis();

  real_t diag_a[mesh.num_cells_in_elem()*num_basis];
  auto diag_M = ViewCArray <real_t> (diag_a, mesh.num_cells_in_elem(), num_basis);
  
  real_t mass_mat_a[num_basis*num_basis];
  auto mass_mat = ViewCArray <real_t> (mass_mat_a, num_basis, num_basis);

#pragma omp simd 
  for(int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){

    // initialize lumped mass //
    for( int cell_lid = 0; cell_lid < mesh.num_cells_in_elem(); cell_lid++){
      int cell_gid = mesh.cells_in_elem(elem_gid, cell_lid);
      for (int node_lid = 0; node_lid < mesh.num_nodes_in_cell(); node_lid++){
        int node_gid = mesh.nodes_in_cell(cell_gid, node_lid);
        node.lumped_mass(node_gid,cell_gid) = 0.0;
	//std::cout << " # corners around node is " << mesh.num_corners_in_node(node_gid) << std::endl;
      }
    }
     
    // compute lumped mass //
    for(int cell_lid = 0 ; cell_lid < mesh.num_cells_in_elem(); cell_lid++){
    
      int cell_gid = mesh.cells_in_elem(elem_gid,cell_lid);
      for(int j = 0; j < num_basis; j++){
        for(int i = 0; i < num_basis; i++){
          mass_mat(i,j) = 0.0;
        }// end loop over i
      }// end loop over j


      
      for(int basis_n = 0; basis_n < num_basis; basis_n++){
      	 for(int basis_m = 0; basis_m < num_basis; basis_m++){
           for(int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_cell(); gauss_lid++){


	     int gauss_gid = mesh.gauss_in_cell(cell_gid, gauss_lid);
           //for (int vertex = 0; vertex < ref_elem.num_basis(); vertex++){
	     mass_mat(basis_m, basis_n) +=  cell_state.density(cell_gid) 
	    				    * ref_elem.ref_nodal_basis(gauss_lid, basis_m) 
					    * ref_elem.ref_nodal_basis(gauss_lid, basis_n) 
					    * mesh.gauss_pt_det_j(gauss_gid)
					    * ref_elem.ref_node_g_weights(gauss_lid);
 	  //   std::cout<< " mass mat entry = "<< mass_mat(basis_m,basis_n) << std::endl;			
           } // end loop over gauss in element

         } // end loop over basis_m         
      } // end loop over basis_n
      
      for (int basis_m = 0; basis_m < num_basis; basis_m++){ 
        for(int basis_n = 0; basis_n < num_basis; basis_n++){
          int node_lid = elem.vert_node_map(basis_n);
          int node_gid = mesh.nodes_in_cell(cell_gid,node_lid);
          node.lumped_mass(node_gid,cell_gid) += mass_mat(basis_m,basis_n);
      //    std::cout << " lumped mass is " << node.lumped_mass(cell_gid, basis_n) << std::endl;
        }
      }
    }// end loop over cell_lid
		
  } // end loop over elem_gid



}// end lumped mass







/*



  for( int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
    
    for (int node_lid = 0; node_lid < mesh.num_nodes_in_elem(); node_lid++){
    
      int node_gid = mesh.nodes_in_elem( elem_gid, node_lid );

        for (int corner_lid = 0; corner_lid < mesh.num_corners_in_node(node_gid); corner_lid++){
        
          
        
        }// end loop over corner_lid
    
    }// end loop over nodes
    
  }// end loop over elem_gid


*/
