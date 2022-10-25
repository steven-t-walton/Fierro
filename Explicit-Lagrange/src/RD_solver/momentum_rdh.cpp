/* momentum_rd.cpp*/

#include<iostream>

#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"

using namespace utils;

void get_momentum_rd(int correction_step){

  num_dim = mesh.num_dim();

  int update = correction_step;
  int prev = correction_step - 1;
#pragma omp simd
  
//    for (int node_gid = 0; node_gid < mesh.num_nodes(); node_gid++){

  for (int elem_gid = 0; elem_gid < mesh.num_nodes(); elem_gid++){
      
    for (int vertex = 0; vertex < ref_elem.num_basis(); vertex++){
      
     int node_lid = ref_elem.vert_node_map(vertex);
     int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);
      
      real_t sum_res = 0.0;

      // Create view of vel_{r+1} and vel_r coeffs //
      auto vel_r = ViewCArray <real_t> (&elem_state.BV_vel_coeffs( prev, elem_gid, vertex, 0), num_dim );
      auto vel = ViewCArray <real_t> (&elem_state.BV_vel_coeffs(update, elem_gid, vertex, 0), num_dim );

      // Update each dim of vel_{r+1} //
      for (int dim = 0; dim < num_dim; dim++){
        // Sum res in cells around node //
        for (int cell_lid = 0; cell_lid < mesh.num_cells_in_node(node_gid); cell_lid++){

          // Get cell_gid //
          int cell_gid = mesh.cells_in_node(node_gid, cell_lid);

	  // Perform summation //
          sum_res += node.nodal_res( node_gid, cell_gid, dim );
        }// end loop over cell_lid

        // Update momentum //

        vel( dim ) = vel_r(dim) - dt*sum_res/node.lumped_mass( node_gid );
	
      }// end loop over dim

    }// end loop over vertex 
  }// end loop over elem_gid

}// end get_momentum_rd()
