#include<iostream>

#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"

using namespace utils;

void get_position_rdh(int correction_step){

  int update = correction_step;

// Update position with v(t^{n+1}) //
#pragma omp simd  
  for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){

    for( int node_lid = 0; node_lid < mesh.num_nodes_in_elem(); node_lid++){ 
      int node_gid = mesh.nodes_in_elem(elem_gid, node_lid );

//    for (int vertex = 0; vertex < ref_elem.num_basis(); vertex++){
//      int node_lid = ref_elem.vert_node_map(vertex);
//     int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);

      auto vel_update = ViewCArray <real_t> ( &node.vel( update, node_gid, 0), mesh.num_dim() );
      //auto vel_n = ViewCArray <real_t> ( &node.vel(0, node_gid, 0 ), mesh.num_dim() );

      for (int dim = 0; dim < mesh.num_dim(); dim++){
        node.coords(update, node_gid, dim) = node.coords(0, node_gid, dim) + dt*vel_update(dim);//0.5*dt * ( vel_update( dim ) + vel_n(dim) );
      }// end loop over dim

    }// end loop over node_lid
  }// end loop over elem_gid
}// end get_position_rdh()
