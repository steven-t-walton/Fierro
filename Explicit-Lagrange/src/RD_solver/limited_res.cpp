#include<iostream>

#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"

using namespace utils;

void get_limited_res(){
    
    for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
      for (int vertex = 0; vertex < ref_elem.num_basis(); vertex++){
        for (int dim = 0; dim < mesh.num_dim(); dim++){
          elem_state.limited_res( elem_gid, vertex, dim ) = 0.0;
        }// end loop over dim
      }// end loop over vertex
    }// end loop over elem_gid
  
    for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
      for (int vertex = 0; vertex < ref_elem.num_basis(); vertex++){
        for (int dim = 0; dim < mesh.num_dim(); dim++){
          //int node_lid = elem.vert_node_map(vertex);
          //int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);
          //for (int elem_node_lid = 0; elem_node_lid < mesh.num_elems_in_node(node_gid); elem_node_lid++){
          //  int elem_node_gid = mesh.elems_in_node(node_gid, elem_node_lid);
      	    elem_state.limited_res( elem_gid, vertex, dim ) = elem_state.psi_coeffs( elem_gid, vertex, dim ) * elem_state.total_res( elem_gid, dim );
          //}// end loop over elem_node_lid
  	}// end loop over dim
      }// end loop over vertex
    }// end loop over elem_gid
/* 
  // check that limited_res sum to total_res over verts in element 
    for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
      for (int dim = 0; dim < mesh.num_dim(); dim++){
	real_t diff = 0.0;
	real_t sum = 0.0;
      	for (int vertex = 0; vertex < ref_elem.num_basis(); vertex++){
          sum += elem_state.limited_res(elem_gid, vertex, dim );
        }// end loop over vertex
	diff = elem_state.total_res(elem_gid, dim) - sum;
	std::cout <<"The difference between total_res and sum of limited res in elem "<< elem_gid << " in dim "<< dim << " is "<< diff << std::endl;
      }// end loop over dim
    }// end loop over elem_gid
*/
}// end get_limited_res
