#include<iostream>
#include<math.h>

#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"

using namespace utils;

void get_total_res(){

  for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
    for (int dim = 0; dim < mesh.num_dim(); dim++){
	elem_state.total_res( elem_gid, dim ) = 0.0;
    }// end loop over dim
  }// end loop over elem_gid

  for (int dim = 0; dim < mesh.num_dim(); dim++){
    for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
      for (int vertex = 0; vertex < ref_elem.num_basis(); vertex++){
	//int node_lid = elem.vert_node_map(vertex);
	//int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);
	//for (int elem_node_lid = 0; elem_node_lid < mesh.num_elems_in_node(node_gid); elem_node_lid++){
          //int elem_node_gid = mesh.elems_in_node(node_gid, elem_node_lid);
	  //for (int vert = 0; vert < ref_elem.num_basis(); vert++){
	    //std::cout << elem_state.nodal_res( elem_node_gid, vert, dim) << std::endl;
            elem_state.total_res( elem_gid, dim ) += elem_state.nodal_res(elem_gid, vertex, dim);
         // }
       // }
      }// end loop over vertex
     // std::cout << "The sum of nodal_res in elem "<< elem_gid << " in dim " << dim << " is " << elem_state.total_res( elem_gid, dim) << std::endl;
    }// end loop over elem_gid
  }// end loop over dim
}// end get_total_res