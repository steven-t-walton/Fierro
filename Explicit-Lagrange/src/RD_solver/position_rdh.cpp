#include<iostream>

#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"

using namespace utils;

void update_position(int t_step){

#pragma omp simd  
/*  
  for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
    for (int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++){
      int gauss_gid = mesh.gauss_in_elem(elem_gid, gauss_lid);
      int node_gid = mesh.nodes_in_elem(elem_gid, gauss_lid);

      real_t interp_a[num_dim];
      for (int i = 0; i < num_dim; i++) interp_a[i] = 0.0;
      auto interp = ViewCArray <real_t> ( &interp_a[0], num_dim);

      for (int dim = 0; dim < mesh.num_dim(); dim++){
        for (int vert = 0; vert < ref_elem.num_basis(); vert++){
          interp(dim) += ref_elem.ref_nodal_basis( gauss_lid, vert ) * elem_state.pos_coeffs(correction_storage-1, elem_gid, vert, dim);
        }// end loop over vertex
      }// end loop over dim
      
      for (int dim = 0; dim < mesh.num_dim(); dim++){
        node.coords(1,node_gid,dim) =  interp(dim);
	mesh.node_coords(node_gid, dim) = node.coords(1, node_gid, dim);
      }

    }// end loop over gauss_lid
 }// end loop over elem_gid
*/

for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
  for (int node_lid = 0; node_lid <mesh.num_nodes_in_elem(); node_lid++){
      int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);
      for (int dim = 0; dim < mesh.num_dim(); dim++){
        node.coords(1, node_gid, dim ) = node.coords(0, node_gid, dim) + 0.5*dt*( node.vel(t_step, node_gid,  dim) + node.vel(0, node_gid, dim) );
        mesh.node_coords(node_gid, dim) = node.coords(1, node_gid, dim);
      }// end loop over dim
  }// end loop over node_gid
}
}// end get_position_rdh()

