/* momentum_rd.cpp*/

#include<iostream>

#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"

using namespace utils;

void get_momentum_rd(int correction_step){
  
  int num_basis = ref_elem.num_basis();
  int num_dim = mesh.num_dim();

  int update = correction_step+1;
  int current = correction_step;

  for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
    
    for (int vertex = 0; vertex < num_basis; vertex++){

      auto vel_update = ViewCArray <real_t> ( &elem_state.BV_vel_coeffs( update, elem_gid, vertex, 0 ), num_dim );
      
      auto vel_r = ViewCArray <real_t> ( &elem_state.BV_vel_coeffs( current, elem_gid, vertex, 0 ), num_dim );

      int node_lid = elem.vert_node_map( vertex );
      int node_gid = mesh.nodes_in_elem( elem_gid, node_lid );
     
      real_t global_lumped_mass = 0.0;

      for (int elem_node = 0; elem_node < mesh.num_elems_in_node(node_gid); elem_node++){
        int elem_node_gid = mesh.elems_in_node(node_gid, elem_node);
        //real_t elem_lumped_mass = 0.0;
        for (int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++){
     	  int gauss_gid = mesh.gauss_in_elem(elem_node_gid, gauss_lid);
	  global_lumped_mass += ref_elem.ref_nodal_basis( gauss_lid, vertex )
		              * ref_elem.ref_node_g_weights( gauss_lid )
		              * mesh.gauss_pt_det_j( gauss_gid );
	}// end loop over gauss_lid
	//global_lumped_mass += elem_lumped_mass;
      }// end loop over elem_node

       
      real_t sum_a[num_dim];
      auto sum = ViewCArray <real_t> ( &sum_a[0], num_dim);
      for (int dim  = 0; dim < num_dim; dim++) sum(dim) = 0.0;

      for (int dim = 0; dim < num_dim; dim ++){
        for (int elem_node = 0; elem_node < mesh.num_elems_in_node(node_gid); elem_node++){
          int elem_node_gid = mesh.elems_in_node(node_gid, elem_node);

	  // low order residual //
          sum(dim) += elem_state.nodal_res(elem_node_gid, vertex, dim);
	  //std::cout << elem_state.nodal_res(elem_node_gid, vertex, dim) << std::endl;

	  // limited residual //
	  //sum(dim) += elem_state.limited_res(elem_node_gid, vertex, dim);
	  //std::cout << elem_state.limited_res(elem_node_gid, vertex, dim) << std::endl;

	}// end loop over elem_node
      }// end loop over dim

      for (int dim = 0; dim < num_dim; dim++){
        vel_update(dim) = vel_r(dim) - sum(dim)*(dt/global_lumped_mass);
        //std::cout << vel_update(dim) << std::endl;
      }// end loop over dim

    }// end loop over vertex
  }// end loop over elem_gid

}// end get_momentum_rd()



/*
 

#pragma omp simd
  
//    for (int node_gid = 0; node_gid < mesh.num_nodes(); node_gid++){

  for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
      
    for (int vertex = 0; vertex < ref_elem.num_basis(); vertex++){
      
      int node_lid = elem.vert_node_map(vertex);
      int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);
      
      real_t sum_res = 0.0;

      // Create view of vel_{r+1} and vel_r coeffs //
      auto vel_r = ViewCArray <real_t> (&elem_state.BV_vel_coeffs( prev, elem_gid, vertex, 0), num_dim );
      auto vel = ViewCArray <real_t> (&elem_state.BV_vel_coeffs(update, elem_gid, vertex, 0), num_dim );
      
      // Get lumped mass //
      real_t lumped_mass = 0.0;
      real_t temp_sum = 0.0;

      for (int cell_lid = 0; cell_lid < mesh.num_cells_in_elem(); cell_lid++){
        int cell_gid = mesh.cells_in_elem(node_gid, cell_lid);
        real_t temp = 0.0;
        for (int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_cell(); gauss_lid++){
          temp += ref_elem.ref_cell_basis(gauss_lid, vertex)
                      * mesh.gauss_cell_pt_det_j(cell_gid)
                      * ref_elem.ref_cell_g_weights(gauss_lid);
        }// end loop over gauss_lid
        temp_sum += temp;
      }// end loop over cells
      
      lumped_mass = temp_sum;


      for (int dim = 0; dim < num_dim; dim++){

      	// Sum res in cells around node //
        for (int cell_lid = 0; cell_lid < mesh.num_cells_in_node(node_gid); cell_lid++){

          // Get cell_gid //
          int cell_gid = mesh.cells_in_node(node_gid, cell_lid);

	  // Perform summation //
          sum_res += elem_state.nodal_res( elem_gid, vertex, cell_gid, dim );
        }// end loop over cell_lid

	// Update momentum //

        vel(dim) = vel_r(dim) - dt*sum_res/lumped_mass;
	
      }// end loop over dim

    }// end loop over vertex 
  }// end loop over elem_gid


 */
