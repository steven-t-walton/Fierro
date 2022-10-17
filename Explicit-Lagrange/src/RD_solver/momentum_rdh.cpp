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
#pragma omp simd
  
    for (int node_gid = 0; node_gid < mesh.num_nodes(); node_gid++){

      real_t sum_res = 0.0;

      // Update each dim of vel_{r+1} //
      for (int dim = 0; dim < num_dim; dim++){
        // Sum res in cells around node //
        for (int cell_lid = 0; cell_lid < mesh.num_cells_in_node(node_gid); cell_lid++){

          // Get cell_gid //
          int cell_gid = mesh.cells_in_node(node_gid, cell_lid);

          // Perform summation //
          sum_res += node.nodal_res(node_gid, cell_gid, dim);
          sum_res = sum_res/node.lumped_mass(node_gid,cell_gid); 
        }// end loop over cell_lid
        //std::cout << " sum of residuals around node " << node_gid << " in cycle "<< cycle <<" is " << sum_res << std::endl;
        // Update momentum //
        node.vel(update, node_gid, dim) = node.vel( correction_step-1, node_gid, dim) - dt*sum_res;
        boundary_rdh(update);
        //sum_res = 0.0;
        //std::cout << " nodal vel at node " << node_gid << " correction step  " << correction_step +1 << " and dim " << dim << " is " << node.vel(correction_step+1, node_gid, dim) <<std::endl; 
      }// end loop over dim
    }//end loop over nodes

#pragma omp simd  
    for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
    
      for(int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++){

        int node_gid = mesh.nodes_in_elem(elem_gid, gauss_lid);

        // get the global id of the gauss point
        int gauss_gid = mesh.gauss_in_elem(elem_gid, gauss_lid);
       
        real_t interp_vel[mesh.num_dim()];
        for(int i = 0; i < mesh.num_dim(); i++){
          interp_vel[i] = 0.0;
        }
      
        for (int dim = 0; dim < mesh.num_dim(); dim++){
          for(int basis_id = 0; basis_id < elem.num_basis(); basis_id++){
          
            int node_basis_id = elem.vert_node_map(basis_id);
            int interp_gid = mesh.gauss_in_elem(elem_gid, node_basis_id);
            interp_vel[dim] += node.vel(update,interp_gid,dim) * ref_elem.ref_nodal_basis(gauss_lid, basis_id);
       
          }// end loop over basis id
        }// end loop over dim
     
        // Save interpolated velocity back to gauss point
        for (int dim = 0; dim < mesh.num_dim(); dim++){
          node.vel(update, node_gid, dim) = interp_vel[dim];
	  interp_vel[dim] = 0.0;
//	  std::cout << " " << std::endl; 
//          std::cout << "assignment to node.vel(t^{r,n}, node_gid, dim) at correction step "<< correction_step+1<< " and node "<< node_gid<< " is " << node.vel(correction_step+1,node_gid, dim) << std::endl;
        }// end loop over dim
       
        if ( update  == num_correction_steps - 1){
             // Update position with v(t^{n+1}) //
#pragma omp simd  
          for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
     
            for(int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++){

              int node_gid = mesh.nodes_in_elem(elem_gid, gauss_lid);

              // get the global id of the gauss point
              int gauss_gid = mesh.gauss_in_elem(elem_gid, gauss_lid);
         
#pragma omp simd
              for (int dim = 0; dim < mesh.num_dim(); dim++){
                node.coords(1, node_gid, dim) = node.coords(0, node_gid, dim) + dt * node.vel(update, node_gid,dim);
              }// end loop over dim

              real_t interp_pos[mesh.num_dim()];
              for(int i=0; i<mesh.num_dim(); i++) interp_pos[i] = 0.0;

              for (int dim = 0; dim < mesh.num_dim(); dim++){
                for(int basis_id = 0; basis_id < elem.num_basis(); basis_id++){

                  int node_basis_id = elem.vert_node_map(basis_id);
                  int interp_gid = mesh.gauss_in_elem(elem_gid, node_basis_id);
                  interp_pos[dim] += node.coords(1,interp_gid,dim) * ref_elem.ref_nodal_basis(gauss_lid, basis_id);

                }// end loop over basis id
              }// end loop over dim

              // Save interpolated position back to gauss point
              for (int dim = 0; dim < mesh.num_dim(); dim++){
                node.coords( 1, node_gid, dim) = interp_pos[dim];
                interp_pos[dim] = 0.0;
                mesh.node_coords(node_gid, dim) = node.coords( 1 , node_gid, dim);
//                 std::cout << " " << std::endl;  
//              std::cout<< "pos at time "<< TIME+dt <<" cycle "<< cycle <<" node " << node_gid << " and dim " << dim <<" is "<< node.coords(update,node_gid,dim) << std::endl;
              }// end loop over dim


     
            } // end loop over gauss_lid
          }// end loop over elem_gid
        };// end if

      } // end loop over gauss_lid
    }// end loop over elem_gid
  
}// end get_momentum_rd()


