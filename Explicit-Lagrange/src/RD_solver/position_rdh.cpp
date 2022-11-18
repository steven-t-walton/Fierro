#include<iostream>

#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"

using namespace utils;

void get_position_rdh(){
  
#pragma omp simd  
  for (int node_gid = 0; node_gid <mesh.num_nodes(); node_gid++){

    // View of vel coeffs //
    //auto vel_update = ViewCArray <real_t> ( &node.vel( 1, node_gid, 0 ), mesh.num_dim() );
    //auto vel_n = ViewCArray <real_t> ( &node.vel( 0, node_gid, 0 ), mesh.num_dim() );

    // View of pos coeffs //
    //auto pos_update = ViewCArray <real_t> ( &node.coords( 1, node_gid, 0 ), mesh.num_dim() );
    //auto pos_n = ViewCArray <real_t> ( &node.coords( 0, node_gid,0 ), mesh.num_dim() );

    for (int dim = 0; dim < mesh.num_dim(); dim++){
      node.coords(1, node_gid, dim ) = node.coords(0, node_gid, dim) + 0.5*dt*( node.vel( 1, node_gid, dim) + node.vel(0, node_gid, dim) );//vel_update( dim ) + vel_n( dim ) );
      mesh.node_coords(node_gid, dim) = node.coords(1, node_gid, dim);
    }// end loop over dim
  }// end loop over node_gid

}// end get_position_rdh()




  //int update = correction_step + 1;
  //int current = correction_step;
  //for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
  //

   // for (int vert = 0; vert < ref_elem.num_basis(); vert++){
     // int node_lid = elem.vert_node_map(vert);
     // int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);

      //auto vel_update = ViewCArray <real_t> ( &elem_state.BV_vel_coeffs( current, elem_gid, vert, 0), mesh.num_dim() );
      //auto vel_n = ViewCArray <real_t> ( &elem_state.BV_vel_coeffs( 0, elem_gid, vert, 0), mesh.num_dim() );
      //
      
      //auto pos_update = ViewCArray <real_t> ( &elem_state.BV_pos_coeffs( update, elem_gid, vert, 0), mesh.num_dim() );
      //auto pos_n = ViewCArray <real_t> ( &elem_state.BV_pos_coeffs( current, elem_gid, vert, 0), mesh.num_dim() );
      //

 // }// end loop over vert
 // }// end loop over elem_gid
