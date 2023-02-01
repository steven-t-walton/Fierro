// computes lumped mass and stores at nodes //

#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "slam.h"
#include "variables.h"


using namespace utils;

void lumped_mass(){
  
  for (int elem_gid = 0; elem_gid  < mesh.num_elems(); elem_gid++){
    for (int vertex = 0; vertex < ref_elem.num_basis(); vertex++){
      int node_lid = ref_elem.vert_node_map( vertex );
      int node_gid = mesh.nodes_in_elem( elem_gid, node_lid );
      
      real_t temp_sum = 0.0;
      for (int cell_lid = 0; cell_lid < mesh.num_cells_in_node(node_gid); cell_lid++){
        int cell_gid = mesh.cells_in_node(node_gid, cell_lid);
        real_t temp = 0.0;
	for (int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_cell(); gauss_lid++){
          temp += ref_elem.ref_cell_basis(gauss_lid, vertex)
		      * mesh.gauss_cell_pt_det_j(cell_gid)
		      * ref_elem.ref_cell_g_weights(gauss_lid);
	}// end loop over gauss_lid
        temp_sum += temp;
      }// end loop over cells
      node.lumped_mass(node_gid) = temp_sum;
      //std::cout << " lumped mass at node " << node_gid << " is " << node.lumped_mass(node_gid) << std::endl;
    }// end loop over vertices
  }// end elem_gid

}// end lumped mass



/*
 

*/



