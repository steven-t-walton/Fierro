// computes lumped mass and stores at nodes //

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

  
  for(int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
 
    for(int cell_lid = 0 ; cell_lid < mesh.num_cells_in_elem(); cell_lid++){
    
      int cell_gid = mesh.cells_in_elem(elem_gid,cell_lid);
      for(int j = 0; j < num_basis; j++){
        for(int i = 0; i < num_basis; i++){
          mass_mat(i,j) = 0.0;
        }// end loop over i
      }// end loop over j

      for (int i = 0; i < num_basis; i++){
        for (int j = 0; j < mesh.num_cells_in_elem(); j++){ diag_M(j,i) = 0.0; }
      }

      
   	for(int basis_n = 0; basis_n < num_basis; basis_n++){
      	  for(int basis_m = 0; basis_m < num_basis; basis_m++){
            for(int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_cell(); gauss_lid++){
              int gauss_gid = mesh.gauss_in_cell(cell_gid, gauss_lid);
              mass_mat(basis_m, basis_n) +=  cell_state.density(cell_gid) 
	     				    * ref_elem.ref_nodal_basis(gauss_lid, basis_m) 
					    * ref_elem.ref_nodal_basis(gauss_lid, basis_n) 
					    * mesh.gauss_pt_det_j(gauss_gid) 
					    * ref_elem.ref_node_g_weights(gauss_lid);
 	      //std::cout<< " mass mat entry = "<< mass_mat(basis_m,basis_n) << std::endl;			
            } // end loop over gauss in element

          //diag_M(cell_gid, basis_n) += mass_mat(basis_m,basis_n);
          int node_lid = elem.vert_node_map(basis_n);
          int node_gid = mesh.nodes_in_cell(cell_gid,node_lid);
          node.lumped_mass(node_gid,cell_gid) += mass_mat(basis_m,basis_n);

          //std::cout << " lumped mass is " << diag_M(cell_gid, basis_n) << std::endl;

          } // end loop over basis_m         
        } // end loop over basis_n

      //node.lumped_mass(node_gid,cell_gid) = diag_M(cell_gid, node_lid);
      //std::cout<< "mass at node " << node_gid << " is " << node.lumped_mass(node_gid,cell_gid) << std::endl;
     
      /*
      for (int vertex_id = 0; vertex_id < num_basis; vertex_id++){
        int node_lid = elem.vert_node_map(vertex_id);
        int node_gid = mesh.nodes_in_cell(cell_gid, node_lid);
        node.mass(node_gid) = diag_M(cell_gid, vertex_id);
        std::cout<< "mass at node " << node_gid << " is " << node.mass(node_gid) << std::endl;
      } // end loop over vertex_id
      */
    }// end loop over cell_lid
		
  } // end loop over elem_gid

}

