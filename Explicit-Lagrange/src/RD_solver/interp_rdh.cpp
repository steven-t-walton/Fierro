#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"

using namespace utils;

void interp_vel(int update){

  for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){

        real_t interp_vel_a[ ref_elem.num_basis()*mesh.num_dim() ];
        auto interp_vel = ViewCArray <real_t> ( &interp_vel_a[0],  ref_elem.num_basis(), mesh.num_dim() );
        for(int j = 0; j < mesh.num_dim(); j++){
          for (int i = 0; i < ref_elem.num_basis(); i++){
            interp_vel(i,j) = 0.0;
          }
        }

        for (int dim = 0; dim < mesh.num_dim(); dim++){
          for(int vert_id = 0; vert_id < ref_elem.num_basis(); vert_id++){

            int node_id = ref_elem.vert_node_map( vert_id );
            int interp_gid = mesh.nodes_in_elem( elem_gid, node_id );
            for (int k = 0; k < ref_elem.num_basis(); k++){          
              interp_vel(vert_id, dim) += node.vel( update, interp_gid, dim ) * ref_elem.ref_nodal_basis( node_id, k );
            }// end loop over k
            node.vel( update, interp_gid, dim ) = interp_vel(vert_id,dim);
            //std::cout << " interpolated vel at dim " << dim << " is " << node.vel( update, interp_gid, dim ) << std::endl;     
          }// end loop over vert_id
        }// end loop over dim

    }// end loop over elem_gid

}// end interp_vel()


void interp_pos(int update){

  for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
  

        real_t interp_pos_a[ ref_elem.num_basis()*mesh.num_dim() ];
        auto interp_pos = ViewCArray <real_t> ( &interp_pos_a[0],  ref_elem.num_basis(), mesh.num_dim() );
        for(int j = 0; j < mesh.num_dim(); j++){
          for (int i = 0; i < ref_elem.num_basis(); i++){
            interp_pos(i,j) = 0.0;
          }
        }

        for (int dim = 0; dim < mesh.num_dim(); dim++){
          for(int vert_id = 0; vert_id < ref_elem.num_basis(); vert_id++){

            int node_id = ref_elem.vert_node_map( vert_id );
            int interp_gid = mesh.nodes_in_elem( elem_gid, node_id );
            for (int k = 0; k < ref_elem.num_basis(); k++){          
              interp_pos(vert_id, dim) += node.coords( update, interp_gid, dim ) * ref_elem.ref_nodal_basis( node_id, k );
            }// end loop over k
            node.coords( update, interp_gid, dim ) = interp_pos(vert_id,dim);
          }// end loop over vert_id
        }// end loop over dim

  }// end loop over elem_gid

}// end interp_pos()


void interp_ie(int update){

}// end interp_ie

