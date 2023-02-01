// computes lumped mass for a node within a cell. used in low order residual //

#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "slam.h"
#include "variables.h"


using namespace utils;

// Creates Global Mass Matrix

void get_cell_mass(){

  int num_basis = ref_elem.num_basis();

  auto diag_M = CArray <real_t>(mesh.num_elems(), num_basis);

  for(int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
    for (int cell_lid = 0; cell_lid < mesh.num_cells_in_elem(); cell_lid++){
      
      int cell_gid = mesh.cells_in_elem(elem_gid, cell_lid);

      auto mass_mat = CArray <real_t>(num_basis, num_basis);

      // Initialize mass matrix to zero 
      for(int i = 0; i < num_basis; i++){
        for(int j = 0; j < num_basis; j++){
          mass_mat(i,j) = 0.0;
        }
      }

      // Fill mass matrix //
      for(int basis_n = 0; basis_n < num_basis; basis_n++){

        for(int basis_m = 0; basis_m < num_basis; basis_m++){
			
          for(int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_cell(); gauss_lid++){

            int gauss_gid = mesh.gauss_in_cell(elem_gid, gauss_lid);

            mass_mat(basis_m, basis_n) += cell_state.density(cell_gid)
					* ref_elem.ref_nodal_basis(gauss_lid, basis_m) 
					* ref_elem.ref_nodal_basis(gauss_lid, basis_n) 
					* mesh.gauss_pt_det_j(gauss_gid) 
					* ref_elem.ref_node_g_weights(gauss_lid);
          } // end loop over gauss in element
        } // end loop over basis_m
      } // end loop over basis_n
      
     
       
    }// end loop over cells	
  } // end loop over the elements


}

