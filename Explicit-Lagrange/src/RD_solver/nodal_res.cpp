/* nodal residual */
/* R^{r,m}_{p(h)} = \sum_{q}M_{qp}(v^{r,m}_p - v^n_p)
 *                                      + \int^{t^m}_{t^n}(Q^{r,m}_p + \int_{V_h}(\grad\varphi\cdot\sigma)dV)dt  */

#include<iostream>
#include<math.h>

#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"


using namespace utils;

void get_nodal_res(real_t sub_dt, int t_step, real_t TIME){
   
  int num_basis = ref_elem.num_basis();
  int num_dim = mesh.num_dim();
 
  // Initialize nodal_res to zero //
#pragma omp simd 
  for (int node_gid = 0; node_gid < mesh.num_nodes(); node_gid++){
    for (int cell_lid = 0; cell_lid < mesh.num_cells_in_node(node_gid);  cell_lid++){
      int cell_gid = mesh.cells_in_node(node_gid, cell_lid);
      for (int dim = 0; dim < num_dim; dim++){
        node.nodal_res(node_gid,cell_gid,dim) = 0.0;
      }
    }
  }



  // Create CArray for vel_bar used in artificial viscosity //
  real_t vel_bar_a[num_dim*(t_step+1)];
  auto vel_bar = ViewCArray <real_t> (&vel_bar_a[0], num_dim, t_step);
  // set to zero //
  // Create CArray for Q //
  real_t Q_a[num_dim*(t_step+1)];
  auto Q = ViewCArray <real_t> (&Q_a[0],num_dim, t_step);

  // Create CArray for sigma //
  real_t sigma_a[num_dim*num_dim*(t_step+1)];
  auto sigma = ViewCArray <real_t> (&sigma_a[0], num_dim, num_dim, t_step);

  // Create CArray to store volume integral over cell of force at each sub time //
  real_t force_a[num_dim*(t_step+1)];
  auto force_cell_volume = ViewCArray <real_t> (&force_a[0], num_dim, t_step);
   
  // Create CArray to store time integral of force and artificial viscosity
  real_t time_int_a[num_dim];
  auto time_integral = ViewCArray <real_t> (&time_int_a[0], num_dim);  

  real_t time_weights[num_correction_steps];
  for (int i = 0; i < num_correction_steps; i++) time_weights[i] = 0.0;
  if (t_step == 0){
    time_weights[0] = 0.5;
    time_weights[1] = 0.5;
  }
  if (t_step == 1){
    time_weights[0] = 0.33;
    time_weights[1] = 1.33;
  }

    int rid_a[ref_elem.num_basis()];
  auto rid = ViewCArray <int> (&rid_a[0], ref_elem.num_basis());
  for (int m=0; m < ref_elem.num_basis(); m++) rid(m) = 0;

#pragma omp simd
  // Loop over elements //
  for(int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
   
    for(int cell_lid = 0; cell_lid < mesh.num_cells_in_elem(); cell_lid++){
      int cell_gid = mesh.cells_in_elem(elem_gid, cell_lid);
      for(int node_lid = 0; node_lid < mesh.num_nodes_in_cell(); node_lid++){
        int node_gid = mesh.nodes_in_cell(cell_gid, node_lid);

        // Create a view of vel and vel_n //
        auto vel = ViewCArray <real_t> (&node.vel(t_step, node_gid, 0), num_dim);
        auto vel_n = ViewCArray <real_t> (&node.vel(0,node_gid,0), num_dim);

        // Initialize Q, sigma, vel_bar volume integral of "force" and time_integral//

        for (int dim = 0; dim < num_dim; dim++){
          for (int prev_times = 0; prev_times <=t_step; prev_times++){
            Q(dim, prev_times) = 0.0;
            vel_bar(dim, prev_times) = 0.0;
            force_cell_volume(dim, prev_times) = 0.0;
            for (int i = 0; i < num_dim; i++){
              sigma(i,dim,prev_times) = 0.0;
            }// end loop over i
          }// end loop over prev_times
          time_integral(dim) = 0.0;
        }// end loop over dim


        // Loop over node_lid (inside the previous loop) in cells to compute vel_bar //
        for (int dim_j = 0; dim_j < num_dim; dim_j++){
          for (int prev_times = 0; prev_times <= t_step; prev_times ++){
            // Loop over node_lid in cells to compute vel_bar //
            for (int node_lid_in_cell = 0; node_lid_in_cell < mesh.num_nodes_in_cell(); node_lid_in_cell++){
              // Get node_gid for each node_lid in cell  //
              int node_gid_from_cell = mesh.nodes_in_cell(cell_gid, node_lid_in_cell);
              vel_bar(dim_j, prev_times) += node.vel(prev_times, node_gid_from_cell, dim_j);
              //std::cout << " vel_bar = " << vel_bar(dim_j, prev_times) << std::endl;   
            }// end loop over dim_j for vel_bar
          }// end loop over prev_times for vel_bar
        }// end loop over nodes in cell
        

	 real_t mass = 0.0;
         real_t mass_a[num_basis];
         auto mass_vec = ViewCArray <real_t> (&mass_a[0], num_basis);
         for (int basis_m = 0; basis_m < num_basis; basis_m++){
           mass_vec(basis_m) = 0.0;
         }// end loop over basis_m
    
         for(int basis_m = 0; basis_m < num_basis; basis_m++){
           for(int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++){
             int gauss_gid = mesh.gauss_in_elem(elem_gid,gauss_lid);//cell(cells_in_node_gid, gauss_lid);
             mass_vec(basis_m) += cell_state.density(cell_gid)//cells_in_node_gid)
                                  *ref_elem.ref_nodal_basis(gauss_lid,basis_m)//cell_basis(gauss_lid, basis_m)
                                  *ref_elem.ref_nodal_basis(gauss_lid, node_lid)//cell_basis(gauss_lid, node_lid)
                                  *mesh.gauss_pt_det_j(gauss_gid)//cell_pt_det_j(gauss_gid)
                                  *ref_elem.ref_node_g_weights(gauss_lid);//cell_g_weights(gauss_lid);
           } // end loop over gauss in element
           mass += mass_vec(basis_m);
         } // end loop over basis_m




        for (int dim_j=0; dim_j < num_dim; dim_j++){ 
          for (int prev_times = 0; prev_times <= t_step; prev_times++){
            real_t alpha_a[3];
            real_t speed = sqrt( node.vel(prev_times, node_gid, 0)*node.vel(prev_times, node_gid, 0) 
			+ node.vel(prev_times, node_gid, 1)*node.vel(prev_times, node_gid, 1)
			+ node.vel(prev_times, node_gid, 2)*node.vel(prev_times, node_gid, 2) );
            alpha_a[0] = speed + cell_state.cs(cell_gid);
            alpha_a[1] = cell_state.cs(cell_gid);
            alpha_a[2] = speed - cell_state.cs(cell_gid);
        
	    real_t alpha = 0.0;


	    for (int i = 0; i < 3; i++){
	      real_t temp1 = alpha_a[0] > alpha_a[1] ? alpha_a[0] : alpha_a[1];
	      alpha = alpha_a[3] > temp1 ? alpha_a[3] : temp1;          	  
            }
        
           // std::cout << "a0 = " << alpha_a[0] << ", a1 = " << alpha_a[1] << ", a2 = " << alpha_a[2] << std::endl;
	   // std::cout << " alpha = " << alpha << std::endl;
		  
	    // Fill Q //
            // Q = alpha_k*(vel - vel_bar) //
            Q(dim_j, prev_times) = alpha
                                   *(node.vel(prev_times,node_gid,dim_j) - vel_bar(dim_j,prev_times));
            //std::cout <<  "Q at dim " << dim_j << " is " << Q(dim_j, prev_times) << std::endl;  
          }// end loop over prev_times for Q
        }// end loop over dim_j for Q 

        
 
       // Loop over cells in node //
       for (int cells_in_node_lid = 0; cells_in_node_lid < mesh.num_cells_in_node(node_gid); cells_in_node_lid++){
         // get cell_gid //
         int cells_in_node_gid = mesh.cells_in_node(node_gid, cells_in_node_lid);
       
       	 for (int dim_j=0; dim_j < num_dim; dim_j++){
	   for (int dim_i = 0; dim_i < num_dim; dim_i++){
             for (int prev_times = 0; prev_times <= t_step; prev_times++){
               // Fill sigma //
               sigma(dim_i, dim_j, prev_times) = cell_state.stress(prev_times, cells_in_node_gid, dim_i, dim_j);
               //std::cout << " sigma at dim_i = "<< dim_i << " and dim_j = "<< dim_j <<" is equal to " << sigma(dim_i,dim_j,prev_times) << std::endl;

             }// end prev_times sigma
           } // end dim_i sigma
         } // end dim_j sigma
        
         
         for (int prev_times = 0; prev_times <= t_step; prev_times++){
           for (int dim_k = 0; dim_k < num_dim; dim_k++){
             for (int dim_j=0; dim_j < num_dim; dim_j++){
               for (int dim_i = 0; dim_i < num_dim; dim_i++){ 
                 for (int gauss_cell_lid = 0; gauss_cell_lid < mesh.num_gauss_in_cell(); gauss_cell_lid++){
               
                   int gauss_gid = mesh.gauss_in_cell(cells_in_node_gid, gauss_cell_lid);
                    
                   force_cell_volume(dim_k, prev_times) += mesh.gauss_cell_pt_jacobian_inverse(gauss_gid, dim_i, dim_j)
                                                           *ref_elem.ref_cell_gradient(gauss_cell_lid, node_lid, dim_i)
                                                           *sigma(dim_j, dim_k, prev_times)
                                                           *ref_elem.ref_cell_g_weights(gauss_cell_lid)
                                                           *mesh.gauss_cell_pt_det_j(gauss_gid);
                   
                 }// end loop over gauss_cell_lid
               }// end dim_i for force_cell_volume
             }// end dim_j for force_cell_volume
             //std::cout << "volume integral of force = "<< force_cell_volume(dim_k,prev_times) << std::endl;
           }// end loop over dim_k
         }// end loop over prev_times for force_cell_volume
         
          
          
         for (int dim_j=0; dim_j < num_dim; dim_j++){
           for (int prev_times = 0; prev_times <= t_step; prev_times++){    
             // Begin time integration of force_cell_volume and Q //
             // int^{t^m}_{t^n}(Q^{r,m}_p + \int_{V_h}(\grad\varphi\cdot\sigma)dV)dt
             // real_t sub_time = prev_times*sub_dt+TIME; 
             time_integral(dim_j) += sub_dt*time_weights[prev_times]*(force_cell_volume(dim_j,prev_times) + Q(dim_j,prev_times));

//	     std::cout<< " force at time_step " << prev_times <<" and dim "<< dim_j  << " is = "<< force_cell_volume(dim_j,prev_times) << std::endl;	
//	     std::cout<< " Q at time_step " << prev_times <<" and dim "<< dim_j  << " is = "<< Q(dim_j,prev_times) << std::endl;                    
//             std::cout << "time integral is = " << time_integral(dim_j) << std::endl;
           }// end loop over prev times for time_integral
         } //end loop over dim_j for time_integral

          
         // Store lo_res with node_gid and cell_gid //
            
         for (int dim = 0; dim < num_dim; dim++){
           node.nodal_res(node_gid, cells_in_node_gid, dim) = mass*(vel(dim) - vel_n(dim)) + time_integral(dim);
//           std::cout << "nodal res value = " << node.nodal_res(node_gid, cells_in_node_gid, dim) << std::endl;
         }// end loop over dim  


       }// end loop over cells_in_node_lid
      

         

      }// end loop over node_lid
    }// end loop over cell_lid
  }// end loop over elements
  //std::cout << "end loop over elements " << std::endl;

}// end get_nodal_res
