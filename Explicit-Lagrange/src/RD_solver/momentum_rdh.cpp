/* momentum_rd.cpp*/

#include<iostream>

#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"

using namespace utils;

void get_momentum_rd(int pred_step, int correction_step){

  num_dim = mesh.num_dim();

#pragma omp simd    
  for (int node_gid = 0; node_gid < mesh.num_nodes(); node_gid++){

    // Create views of vel_r //
    auto vel_r = ViewCArray <real_t> (&node.vel(correction_step, node_gid, 0), num_dim);
    /*
    // Create CArray  to store summed residual //
    auto sum_res = CArray <real_t> (num_dim);

    // Initialize to zero //
    for (int dim = 0; dim < num_dim; dim++){

      sum_res(dim) = 0.0;

    }// end sum_res initialization
    */

    // Update each dim of vel_{r+1} //
    for (int dim = 0; dim < num_dim; dim++){
      real_t sum_res = 0.0;
      // Sum res in cells around node //
      for (int cell_lid = 0; cell_lid < mesh.num_cells_in_node(node_gid); cell_lid++){

        // Get cell_gid //
        int cell_gid = mesh.cells_in_node(node_gid, cell_lid);

        // Perform summation //
        sum_res += node.nodal_res(node_gid, cell_gid, dim);
        sum_res = sum_res/node.lumped_mass(node_gid,cell_gid); 
      }// end loop over cell_lid

      //std::cout << " summed over residuals " << std::endl;
      //std::cout << "dim = " << dim << std::endl;
      //std::cout << " sum_res around node " << node_gid << " is " << sum_res(dim) << std::endl;

      // Update momentum //
      //std::cout << "mass at node "<< node_gid <<" is "<< node.mass(node_gid) << std::endl;
      node.vel(pred_step,node_gid,dim) = vel_r(dim) - dt*sum_res;
      //std::cout << "assignment to node.vel(corr_step, node_gid, dim) " << vel_r(dim) - dt*sum_res << std::endl;
      
    }// end loop over dim
  }//end loop over nodes

  for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
    
    for(int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++){

      int node_gid = mesh.nodes_in_elem(elem_gid, gauss_lid);

      // get the global id of the gauss point
      int gauss_gid = mesh.gauss_in_elem(elem_gid, gauss_lid);

      real_t interp_vel[mesh.num_dim()];
      for(int i=0; i<mesh.num_dim(); i++) interp_vel[i] = 0.0;

      // Sum over the basis times the velocity defined at the basis vertex position
      for (int dim = 0; dim < mesh.num_dim(); dim++){
        for(int basis_id = 0; basis_id < elem.num_basis(); basis_id++){

          int node_basis_id = elem.vert_node_map(basis_id);
          int interp_gid = mesh.gauss_in_elem(elem_gid, node_basis_id);
          interp_vel[dim] += node.vel(pred_step, interp_gid, dim) * ref_elem.ref_nodal_basis(gauss_lid, basis_id);
        }// end loop over basis id
      }// end loop over dim

      // Save interpolated velocity back to gauss point
      for (int dim = 0; dim < mesh.num_dim(); dim++){
        node.vel(pred_step, gauss_gid, dim) =  interp_vel[dim];
      }// end loop over dim
    } // end loop over gauss_lid
  }// end loop over elem_gid

  //std::cout << "finished get momentum" << std::endl;

}// end get_momentum_rd()
