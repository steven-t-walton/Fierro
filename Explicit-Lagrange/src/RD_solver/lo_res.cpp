/* nodal residual */
/* R^{r,m}_{p(h)} = \sum_{q}M_{qp}(v^{r,m}_p - v^n_p)
 *                                      + \int^{t^m}_{t^n}(Q^{r,m}_p + \int_{V_h}(\grad\varphi\cdot\sigma)dV)dt  */

#include<iostream>

#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"
#include "jacobi_polynomials.h"

using namespace utils;

void get_lo_res(real_t sub_dt, int t_step, real_t sub_time){
  
  num_dim = mesh.num_dim();

  // Loop over elements //
  for(int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
    // Loop over nodes in element //
    for (int node_lid = 0; node_lid < mesh.num_nodes_in_elem(); node_lid++){
      // Get node global id //
      int node_gid = mesh.nodes_in_elem(elem_gid, node_lid); 
      //std::cout << " node lid = " << node_lid << " in elem " << elem_gid << " with node_gid = " << node_gid << std::endl;
      // Create a view of vel and vel_n //
      auto vel = ViewCArray <real_t> (&node.vel(t_step, node_gid, 0), num_dim);
      auto vel_n = ViewCArray <real_t> (&node.vel(0,node_gid,0), num_dim);
      
      // Create CArray for vel_bar used in artificial viscosity //
      real_t vel_bar_a[num_dim*t_step];
      auto vel_bar = ViewCArray <real_t> (vel_bar_a, num_dim, t_step);

      // Create CArray for Q //
      real_t Q_a[num_dim*t_step];
      auto Q = ViewCArray <real_t> (Q_a,num_dim, t_step);

      // Create CArray for sigma //
      real_t sigma_a[num_dim*num_dim*t_step];
      auto sigma = ViewCArray <real_t> (sigma_a, num_dim, num_dim, t_step);
 
      // Create CArray to store volume integral over cell of force at each sub time //
      auto force_cell_volume = CArray <real_t> (num_dim, t_step);

      real_t alpha_k = 1.0; // set alpha_k to 1 for now. <--- (CHANGE THIS)
      
      // Create CArray to store time integral of force and artificial viscosity
      auto time_integral = CArray <real_t> (num_dim);  
		         
      // Initialize Q, sigma and vel_bar //
#pragma omp simd
      for (int dim = 0; dim < num_dim; dim++){
        for (int prev_times = 0; prev_times <=t_step; prev_times++){
          Q(dim, prev_times) = 0.0;
	  vel_bar(dim, prev_times) = 0.0;
          for (int i = 0; i < num_dim; i++){
            sigma(i,dim,prev_times) = 0.0;
          }
        }
      }
      // Loop over cells in node //
      for (int cell_lid = 0; cell_lid < mesh.num_cells_in_node(node_gid); cell_lid++){
        // get cell_gid //
        int cell_gid = mesh.cells_in_node(node_gid, cell_lid);
        
        // Loop over node_lid in cells to compute vel_bar //
        for (int dim_j = 0; dim_j < num_dim; dim_j++){
          for (int prev_times = 0; prev_times <= t_step; prev_times ++){
            // Loop over node_lid in cells to compute vel_bar //
            for (int node_lid_in_cell = 0; node_lid_in_cell<mesh.num_nodes_in_cell(); node_lid_in_cell++){
              // Get node_gid for each node_lid in cell  //
              int node_gid_from_cell = mesh.nodes_in_cell(cell_gid, node_lid_in_cell);
              vel_bar(dim_j, prev_times) += node.vel(prev_times, node_gid_from_cell, dim_j);
              std::cout << " vel_bar = " << vel_bar(dim_j, prev_times) << std::endl;    
            } // end loop over dim_j for vel_bar
          }// end loop over prev_times for vel_bair
        } // end loop over nodes in cell
         
        //std::cout << "vel_bar filled for node" << node_gid << std::endl;
	
	for (int dim_j=0; dim_j < num_dim; dim_j++){
	  for (int dim_i = 0; dim_i < num_dim; dim_i++){
            for (int prev_times = 0; prev_times <= t_step; prev_times++){
              // Fill sigma //
              sigma(dim_i, dim_j, prev_times) = cell_state.stress(prev_times, cell_gid, dim_i, dim_j);
              //std::cout << " sigma at dim_i = "<< dim_i << " and dim_j = "<< dim_j <<" is equal to " << sigma(dim_i,dim_j,prev_times) << std::endl;

            }// end prev_times sigma
          } // end dim_i sigma
        } // end dim_j sigma
        
        for (int prev_times = 0; prev_times <= t_step; prev_times++){
          for (int dim_j=0; dim_j < num_dim; dim_j++){
            for (int dim_i = 0; dim_i < num_dim; dim_i++){
              for (int basis_id = 0; basis_id < elem.num_basis(); basis_id++){                               
                force_cell_volume(dim_j, prev_times) += 0.0;//mesh.gauss_cell_pt_jacobian_inverse(cell_gid, dim_i, dim_j)
                                                        //*ref_elem.ref_cell_gradient(cell_lid, basis_id, dim_j)
                                                        //*sigma(dim_j, dim_j, prev_times)
                                                        //*ref_elem.ref_cell_g_weights(cell_lid)
                                                        //*mesh.gauss_cell_pt_det_j(cell_gid);
              }// end loop over basis id  
            }// end dim_i for force_cell_volume
          }// end dim_j for force_cell_volume
        }// end loop over prev_times for force_cell_volume
  
        for (int dim_j=0; dim_j < num_dim; dim_j++){ 
          for (int prev_times = 0; prev_times <= t_step; prev_times++){
            // Fill Q //
            Q(dim_j, prev_times) = alpha_k
	                            *(node.vel(prev_times,node_gid,dim_j) - vel_bar(dim_j,prev_times)); 
          }// end loop over prev_times for Q
        }// end loop over dim_j for Q 
        
        for (int dim_j=0; dim_j < num_dim; dim_j++){
          for (int prev_times = 0; prev_times <= t_step; prev_times++){    
            // Begin time integration of force_cell_volume and Q //
            // int^{t^m}_{t^n}(Q^{r,m}_p + \int_{V_h}(\grad\varphi\cdot\sigma)dV)dt 
            time_integral(dim_j) += 0.5*(force_cell_volume(dim_j,prev_times)+Q(dim_j,prev_times));
          }// end loop over prev times for time_integral
        } //end loop over dim_j for time_integral
 
        // Store lo_res with node_gid and cell_gid //
        for (int dim = 0; dim < num_dim; dim++){
          // std::cout << " filling node.lo_res for dim = " << dim << std::endl;
          node.lo_res(node_gid, cell_gid, dim) = vel(dim) - vel_n(dim) + time_integral(dim);//cell_state.lumped_mass(cell_gid,node_gid);
          std::cout << "for node_gid "<< node_gid << " cell_gid "<< cell_gid <<" and dim "<< dim << " lo_res is "<<node.lo_res(node_gid, cell_gid, dim) << std::endl;
          //std::cout << " and lumped_mass is " << cell_state.lumped_mass(cell_gid, node_gid) << std::endl;

        }// end loop over dim  
        //std::cout << "end loop over dim " << std::endl;

      }// end loop over cell_lid      
     // std::cout << "end loop over cell_lid for node "<< node_gid << std::endl;

    }// end loop over node_lid
   // std::cout << " end loop over node_lid for element "<< elem_gid <<std::endl;

  }// end loop over elements
  //std::cout << "end loop over elements " << std::endl;

}// end get_lo_res
