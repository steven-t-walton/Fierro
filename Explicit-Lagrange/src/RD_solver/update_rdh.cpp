#include<iostream>

#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"

void update_coeffs(){
  
  int update = num_correction_steps - 1;

  for (int node_gid = 0; node_gid < mesh.num_nodes(); node_gid++) {
      
    auto vel = ViewCArray <real_t> (&node.vel( update, node_gid, 0), mesh.num_dim() );
    auto vel_n = ViewCArray <real_t> (&node.vel( 0, node_gid, 0), mesh.num_dim() );

    for (int dim = 0; dim < 3; dim++){
      node.coords(0, node_gid, dim) = node.coords(1, node_gid, dim);
      mesh.node_coords(node_gid, dim) = node.coords(1, node_gid, dim);

      //std::cout<< "pos at time "<< TIME+dt <<" cycle "<< cycle << " node " << node_gid << " and dim " << dim <<" is "<< node.coords(0,node_gid,dim) << std::endl;
      // Update velocity
      
      vel_n( dim ) = vel( dim );

      // update velocity at boundary appropriately //
      boundary_rdh(0);

//      std::cout << " " << std::endl;
//      std::cout<< "vel at time "<< TIME+dt <<" cycle "<< cycle << " node " << node_gid << " and dim " << dim <<" is "<< vel_n(dim) << std::endl;
//      std::cout << " cycle = " << cycle << std::endl;

      for ( int t_step = 1; t_step <= update; t_step++){
        node.vel(t_step, node_gid, dim) = vel_n( dim );
      }
    }// end loop over dim      
  }// end loop over node_gid 

#pragma omp simd     
  for (int cell_gid = 0; cell_gid < mesh.num_cells(); cell_gid++){

    cell_state.ie( 0, cell_gid ) = cell_state.ie( 1, cell_gid );
       
    for (int dim_j = 0; dim_j < mesh.num_dim(); dim_j++){
      for (int dim_i = 0; dim_i < mesh.num_dim(); dim_i++){
        cell_state.stress(0, cell_gid, dim_i, dim_j) = cell_state.stress(1, cell_gid, dim_i, dim_j);
      }
    }
    //cell_state.total_energy( 0, cell_gid ) = cell_state.total_energy( 1, cell_gid );
  }

}


