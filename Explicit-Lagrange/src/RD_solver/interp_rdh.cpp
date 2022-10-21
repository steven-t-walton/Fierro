#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"

using namespace utils;

void interp_vel(int update){
  for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
      for (int node_lid = 0; node_lid <mesh.num_nodes_in_elem(); node_lid++){

        int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);

        real_t interp_vel[mesh.num_dim()];
        for(int i = 0; i < mesh.num_dim(); i++){
          interp_vel[i] = 0.0;
        }

        for (int dim = 0; dim < mesh.num_dim(); dim++){
          for(int basis_id = 0; basis_id < ref_elem.num_basis(); basis_id++){

            int node_basis_id = elem.vert_node_map( basis_id );
            int interp_gid = mesh.nodes_in_elem( elem_gid, node_basis_id );

            interp_vel[dim] += node.vel( update, interp_gid, dim ) * ref_elem.ref_nodal_basis( node_lid, basis_id );

            }// end loop over basis id
  //        std::cout << " interpolated vel at dim " << dim << " is " << interp_vel[dim] << std::endl;     
        }// end loop over dim

        // Save interpolated velocity back to gauss point
        for (int dim = 0; dim < mesh.num_dim(); dim++){
          node.vel(update, node_gid, dim) = interp_vel[dim];
        }// end loop over dim

      }// end loop over node_lid
    }// end loop over elem_gid
}// end interp_vel()


void interp_pos(int update){
 for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
    for(int node_lid = 0; node_lid < mesh.num_nodes_in_elem(); node_lid++){
      int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);

      real_t interp_pos[mesh.num_dim()];
      for(int i=0; i<mesh.num_dim(); i++) interp_pos[i] = 0.0;

        for (int dim = 0; dim < mesh.num_dim(); dim++){
          for(int basis_id = 0; basis_id < elem.num_basis(); basis_id++){

            int node_basis_id = elem.vert_node_map(basis_id);
            int interp_gid = mesh.nodes_in_elem(elem_gid, node_basis_id);

            interp_pos[dim] += node.coords(update, interp_gid, dim) * ref_elem.ref_nodal_basis(node_lid, basis_id );

          } // end loop over the basis
        } // end loop over dimension

        for (int dim = 0; dim < mesh.num_dim(); dim++){
          node.coords(update, node_gid, dim) = interp_pos[dim];
        }

    }// end loop over node_lid

  }// end loop over elem_gid
}// end interp_pos()


void interp_ie(int update){

}// end interp_ie
