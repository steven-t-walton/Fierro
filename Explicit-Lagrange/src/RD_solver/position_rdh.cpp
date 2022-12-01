#include<iostream>

#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"

using namespace utils;

void update_position(){

 
#pragma omp simd  
  for (int node_gid = 0; node_gid <mesh.num_nodes(); node_gid++){
    for (int dim = 0; dim < mesh.num_dim(); dim++){
      node.coords(1, node_gid, dim ) = node.coords(0, node_gid, dim) + 0.5*dt*( node.vel( 1, node_gid, dim) + node.vel(0, node_gid, dim) );
      mesh.node_coords(node_gid, dim) = node.coords(1, node_gid, dim);
    }// end loop over dim
  }// end loop over node_gid

/*
  for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
    for (int vertex = 0; vertex < ref_elem.num_basis(); vertex++){
      for (int dim = 0; dim < mesh.num_dim(); dim++){
        elem_state.pos_coeffs(num_correction_steps, elem_gid, vertex, dim) = elem_state.pos_coeffs(0, elem_gid, vertex, dim)
	                                                                     + 0.5*dt*( elem_state.vel_coeffs(num_correction_steps, elem_gid, vertex, dim) + elem_state.vel_coeffs(0, elem_gid, vertex, dim) );
      }// end loop over dim
    }// end loop over vertex
  }// end loop over elem_gid


  for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
    for (int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++){
      int gauss_gid = mesh.gauss_in_elem(elem_gid, gauss_lid);
      int node_gid = mesh.nodes_in_elem(elem_gid, gauss_lid);

      real_t interp_a[num_dim];
      for (int i =0; i < num_dim; i++) interp_a[i] =0.0;
      auto interp = ViewCArray <real_t> ( &interp_a[0], num_dim);

      for (int dim = 0; dim < mesh.num_dim(); dim++){
        for (int vert = 0; vert < ref_elem.num_basis(); vert++){
          interp( dim) += ref_elem.ref_nodal_basis( gauss_lid, vert ) * elem_state.pos_coeffs(num_correction_steps, elem_gid, vert, dim);
        }// end loop over vertex
      }// end loop over dim

      for (int dim = 0; dim < num_dim; dim++){
        node.coords(1, node_gid, dim) = interp(dim);
	mesh.node_coords(node_gid, dim) = node.coords(1,node_gid, dim);
      }

    }// end loop over gauss_lid
  }// end loop over elem_gid

*/
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
