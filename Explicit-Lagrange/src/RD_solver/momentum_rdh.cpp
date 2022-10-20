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

  int rid_a[ref_elem.num_basis()];
  auto rid = ViewCArray <int> (&rid_a[0], ref_elem.num_basis());
  for (int m=0; m < ref_elem.num_basis(); m++) rid(m) = 0;

#pragma omp simd
  int ind = 0;
  for (int k = 0; k <  cbrt(ref_elem.num_basis()); k++){
    for (int j = 0; j < cbrt(ref_elem.num_basis()); j++){
      for (int i = 0; i < cbrt(ref_elem.num_basis()); i++){
        rid(ind) = ref_elem.node_rid(i,j,k);
          //std::cout << "node_rid at i = " << i << " j = " << j << " k = " << k << " is " << rid(ind) << std::endl;
        ind++;
      };// end loop over i
    };// end loop over j
  };// end loop over k

#pragma omp simd
  
    for (int node_gid = 0; node_gid < mesh.num_nodes(); node_gid++){

      real_t sum_res = 0.0;
      
      auto vel_r = ViewCArray <real_t> (&node.vel( prev, node_gid, 0), num_dim);

      // Update each dim of vel_{r+1} //
      for (int dim = 0; dim < num_dim; dim++){
        // Sum res in cells around node //
        for (int cell_lid = 0; cell_lid < mesh.num_cells_in_node(node_gid); cell_lid++){

          // Get cell_gid //
          int cell_gid = mesh.cells_in_node(node_gid, cell_lid);

	  // Perform summation //
	  //std::cout << " nodal_res at cell_gid " << cell_gid << " is "  << node.nodal_res(node_gid, cell_gid, dim ) << std::endl;
          sum_res += node.nodal_res( node_gid, cell_gid, dim );
          sum_res = sum_res/node.lumped_mass( node_gid, cell_gid ); 
        }// end loop over cell_lid

       // std::cout << " sum of residuals around node " << node_gid << " in cycle "<< cycle <<" is " << sum_res << std::endl;
        
       // Update momentum //
//	std::cout << " " << std::endl;
//	std::cout << " vel_r at dim " << dim << " is " << vel_r(dim) << std::endl;
//	std::cout << " " << std::endl;

        node.vel(update, node_gid, dim) = vel_r(dim) - dt*sum_res;
	
//	std::cout << " vel before interpolation at " << node_gid << " and dim " << dim << " is " << node.vel(update, node_gid, dim) << std::endl;
//	std::cout << " " << std::endl;

      }// end loop over dim
    }//end loop over nodes

    boundary_rdh(update);

#pragma omp simd  
    for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
    
      for (int node_lid = 0; node_lid <mesh.num_nodes_in_elem(); node_lid++){//for(int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++){

        int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);

        real_t interp_vel[mesh.num_dim()];
        for(int i = 0; i < mesh.num_dim(); i++){
          interp_vel[i] = 0.0;
        }
      
        for (int dim = 0; dim < mesh.num_dim(); dim++){
          for (int basis_n = 0; basis_n < ref_elem.num_basis(); basis_n++){
            for(int basis_id = 0; basis_id < ref_elem.num_basis(); basis_id++){
          
              int node_basis_id = ref_elem.vert_node_map( basis_id );
              int interp_gid = mesh.nodes_in_elem( elem_gid, node_basis_id );

	      auto vel_update = ViewCArray <real_t> (&node.vel(update, interp_gid, 0), num_dim);

              interp_vel[dim] += vel_update( dim ) * ref_elem.ref_nodal_basis( rid(basis_id), basis_n);
            }// end loop over basis id
	  }// end loop over basis_n
  //        std::cout << " interpolated vel at dim " << dim << " is " << interp_vel[dim] << std::endl;     
	}// end loop over dim
     
        // Save interpolated velocity back to gauss point
        for (int dim = 0; dim < mesh.num_dim(); dim++){
          node.vel(update, node_gid, dim) = interp_vel[dim];
          interp_vel[dim] = 0.0;	  
        }// end loop over dim
      } // end loop over gauss_lid
    }// end loop over elem_gid
    
}// end get_momentum_rd()


