#include<iostream>

#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"

using namespace utils;

void interp_position(){
  for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
    for (int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++){
      int gauss_gid = mesh.gauss_in_elem(elem_gid, gauss_lid);
      int node_gid = mesh.nodes_in_elem(elem_gid, gauss_lid);

      CArray <real_t> interp( num_dim );
      for (int i = 0; i < num_dim; i++) interp(i) = 0.0;

      for (int dim = 0; dim < mesh.num_dim(); dim++){
        for (int vert = 0; vert < ref_elem.num_basis(); vert++){
	  int interp_lid = elem.vert_node_map(vert);
	  int interp_gid = mesh.nodes_in_elem(elem_gid, interp_lid);
          interp(dim) += ref_elem.ref_nodal_basis( gauss_lid, vert ) * node.coords(num_correction_steps, interp_gid, dim);
        }// end loop over vertex
      }// end loop over dim
      
      for (int dim = 0; dim < mesh.num_dim(); dim++){
        node.coords(num_correction_steps,node_gid,dim) =  interp(dim);
	mesh.node_coords(node_gid, dim) = interp(dim);//node.coords(num_correction_steps, node_gid, dim);
      }

    }// end loop over gauss_lid
 }// end loop over elem_gid

}// end get_position_rdh()

