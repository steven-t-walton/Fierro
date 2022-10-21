#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"

using namespace utils;

void interp_vel(int update){

  for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
/* 
     for (int node_lid = 0; node_lid <mesh.num_nodes_in_elem(); node_lid++){

        int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);
*/
        real_t interp_vel[mesh.num_dim()];
        for(int i = 0; i < mesh.num_dim(); i++){
          interp_vel[i] = 0.0;
        }

        for (int dim = 0; dim < mesh.num_dim(); dim++){
          for(int vert_id = 0; vert_id < ref_elem.num_basis(); vert_id++){

            int node_id = ref_elem.vert_node_map( vert_id );
            int interp_gid = mesh.nodes_in_elem( elem_gid, node_id );
            for (int k = 0; k < ref_elem.num_basis(); k++){          
              interp_vel[dim] += node.vel( update, interp_gid, dim ) * ref_elem.ref_nodal_basis( node_id, k );
            }// end loop over k
            node.vel( update, interp_gid, dim ) = interp_vel[dim];
            interp_vel[dim] = 0.0;
          }// end loop over vert_id
          std::cout << " interpolated vel at dim " << dim << " is " << interp_vel[dim] << std::endl;     
        }// end loop over dim
        
/*
        for (int dim = 0; dim < mesh.num_dim(); dim++){
          for(int vertex = 0; vertex < elem.num_basis(); vertex++){
            int node_id = ref_elem.vert_node_map( vertex );
            int interp_gid = mesh.nodes_in_elem( elem_gid, node_id );

            node.vel( update , interp_gid, dim) = interp_vel[dim];
            interp_vel[dim] = 0.0;
          }//end loop over vertex
        }//end loop over num_dim
*/
/*
        // Save interpolated velocity back to gauss point
        for (int dim = 0; dim < mesh.num_dim(); dim++){
          node.vel(update, node_gid, dim) = interp_vel[dim];
        }// end loop over dim

      }// end loop over node_lid
*/
    }// end loop over elem_gid


}// end interp_vel()


void interp_pos(int update){

 for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
 /*   for(int node_lid = 0; node_lid < mesh.num_nodes_in_elem(); node_lid++){
      int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);
*/
      real_t interp_pos[mesh.num_dim()];
      for(int i=0; i<mesh.num_dim(); i++) interp_pos[i] = 0.0;

        for (int dim = 0; dim < mesh.num_dim(); dim++){
          for(int vert_id = 0; vert_id < ref_elem.num_basis(); vert_id++){

            int node_id = ref_elem.vert_node_map(vert_id);
            int interp_gid = mesh.nodes_in_elem(elem_gid, node_id);

            interp_pos[dim] += node.coords(update, interp_gid, dim) * ref_elem.ref_nodal_basis(node_id, vert_id );

          } // end loop over vert_id
        } // end loop over dimension


        for (int dim = 0; dim < mesh.num_dim(); dim++){
          for(int vertex = 0; vertex < elem.num_basis(); vertex++){
            int node_id = elem.vert_node_map( vertex );
            int interp_gid = mesh.nodes_in_elem( elem_gid, node_id );

            node.coords( update , interp_gid, dim) = interp_pos[dim];
            interp_pos[dim] = 0.0;
          }//end loop over vertex
        }//end loop over num_dim
/*
        for (int dim = 0; dim < mesh.num_dim(); dim++){
          node.coords(update, node_gid, dim) = interp_pos[dim];
        }

    }// end loop over node_lid
*/
  }// end loop over elem_gid

}// end interp_pos()


void interp_ie(int update){

}// end interp_ie
