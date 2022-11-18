#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"

#define PI 3.14159265

using namespace utils;

void interp_vel(){
  for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
    for (int node_lid = 0; node_lid < mesh.num_nodes_in_elem(); node_lid++){
      int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);
      
      for (int dim = 0; dim < mesh.num_dim(); dim++){
        node.vel(1, node_gid, dim) = 0.0;
      }// end loop over dim


      for (int dim = 0; dim < mesh.num_dim(); dim++){
      	for (int vertex = 0; vertex < ref_elem.num_basis(); vertex++){
	    node.vel( 1, node_gid, dim ) += ref_elem.ref_nodal_basis( node_lid, vertex ) * elem_state.BV_vel_coeffs( num_correction_steps, elem_gid, vertex, dim );
	}// end loop over vertex
      }// end loop over dim



      node.vel(1, node_gid, 0) = sin(PI * mesh.node_coords(node_gid, 0)) * cos(PI * mesh.node_coords(node_gid, 1));
      node.vel(1, node_gid, 1) = -1.0*cos(PI * mesh.node_coords(node_gid, 0)) * sin(PI * mesh.node_coords(node_gid, 1));
      node.vel(1, node_gid, 2) = 0.0;

      
/*
      std::cout << node.vel(1, node_gid, 0) - sin(PI * mesh.node_coords(node_gid, 0)) * cos(PI * mesh.node_coords(node_gid, 1))<< std::endl;
      std::cout << node.vel(1, node_gid, 1) + 1.0*cos(PI * mesh.node_coords(node_gid, 0)) * sin(PI * mesh.node_coords(node_gid, 1))<< std::endl; 
      std::cout << node.vel(1, node_gid, 2) << std::endl; 
*/    
    }// end loop over node_lid
  }// end loop over elem_gid

}// end interp_vel()

/*
void interp_pos(int update){

  for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
    for (int node_lid = 0; node_lid < mesh.num_nodes_in_elem(); node_lid++){
      int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);

      for (int dim = 0; dim < mesh.num_dim(); dim++){
        node.coords(1, node_gid, dim) = 0.0;
      	for (int vertex = 0; vertex < ref_elem.num_basis(); vertex++){
          node.coords( 1, node_gid, dim ) += ref_elem.ref_nodal_basis( node_lid, vertex ) * elem_state.BV_pos_coeffs( update, elem_gid, vertex, dim );
	}// end loop over vertex
      }// end loop over dim
      
    }// end loop over node_lid
  }// end loop over elem_gid
  
}// end interp_pos()
*/

void interp_ie(int update){
  
}// end interp_ie

