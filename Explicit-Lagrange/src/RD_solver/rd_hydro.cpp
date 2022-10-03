/* DeC scheme for RD */


#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"

using namespace utils;

void rd_hydro(){
  
  for (cycle = 1; cycle <= cycle_stop; cycle++){
    
    //std::cout<<"cycle = "<<cycle<<std::endl;

    if (stop_calc == 1) break;
    
        
    if (cycle == 1){
          
      real_t BV_coeffs_a[mesh.num_elems()*elem.num_basis()*mesh.num_dim()];
      auto BV_coeffs = ViewCArray <real_t> (&BV_coeffs_a[0], mesh.num_elems(), elem.num_basis(), mesh.num_dim());

      // initialize control coeffs to zero //
      for( int k = 0; k < mesh.num_dim(); k ++){
        for (int j = 0; j < elem.num_basis(); j++){
          for (int i = 0; i < mesh.num_elems(); i++){
             BV_coeffs(i,j,k) =  0.0;
          }// end loop over i
        }// end loop over j
      }// end loop over k
      
      for (int elem_gid = 0 ; elem_gid < mesh.num_elems(); elem_gid++){
        for (int dim = 0; dim < mesh.num_dim(); dim++){
          for(int basis_id = 0; basis_id < elem.num_basis(); basis_id++){
         
            for (int k = 0; k < elem.num_basis(); k++){
              int node_basis_id = elem.vert_node_map(k);
              int interp_gid = mesh.gauss_in_elem(elem_gid, node_basis_id);
              BV_coeffs(elem_gid,basis_id,dim) += elem_state.BV_mat_inv(basis_id, k)*node.vel(0,interp_gid ,dim); 
            }// end loop over k
                
            //std::cout << "Control coeffs in elem  " << elem_gid << " and dim " << dim << " are " << BV_coeffs(elem_gid, basis_id, dim) << std::endl;
            int node_id = elem.vert_node_map(basis_id);
            int node_gid_for_control_coeff = mesh.nodes_in_elem(elem_gid, node_id);
            node.vel(0,node_gid_for_control_coeff,dim) = BV_coeffs(elem_gid, basis_id, dim);
          }// end loop over basis id
        }// end loop over dim
      }// end loop over elem_gid 
      
      track_rdh(ke, ie, 0);
      te_0 = ie + ke;
    }// end if

    get_timestep();

    dt = fmin(dt, (graphics_time - TIME)+fuzz);
    
    
    { // Time integration scope //
     
      // DeC update //
      for (int correction_step = 0; correction_step < num_correction_steps; correction_step++){
        //std::cout<< " correction step = "<< correction_step << std::endl;          
        real_t sub_dt = (correction_step+1)*dt/num_correction_steps;
        //real_t sub_time = TIME + sub_dt;          
           
	get_stress(correction_step); // assign values to stress
        //std::cout << "called get stress" << std::endl;
    
        //std::cout << "calling get_nodal_res" << std::endl;
        get_nodal_res(sub_dt, correction_step, TIME);

        //std::cout << "calling lumped_mass " << std::endl;
        lumped_mass();// used in momentum
        
        //std::cout << "calling prediction step" << std::endl;
       // prediction_step(sub_dt, pred_step);
          
        //std::cout << " updating momentum " << std::endl;
        get_momentum_rd(num_correction_steps, correction_step);
        
        //std::cout << "Calculating Jacobian at gauss points" << std::endl;
        get_gauss_pt_jacobian(mesh, ref_elem);

        //std::cout << "Before volume from Jacobian"  << std::endl;
        get_vol_jacobi(mesh, ref_elem); 
        
        //std::cout << " calling get state " << std::endl;
        get_state();
      }//end correction steps
      

      //Update position
      for (int node_gid = 0; node_gid < mesh.num_nodes(); node_gid++) {
        // create view of the nodal velocities
        //auto vel = ViewCArray <real_t> (&node.vel(num_correction_steps, node_gid, 0), num_dim);
        for (int dim = 0; dim < 3; dim++){
          node.coords(1, node_gid, dim) = node.coords(0, node_gid, dim) +  dt*node.vel(num_correction_steps, node_gid,dim);//*(vel(dim));
          mesh.node_coords(node_gid, dim) = node.coords(1, node_gid, dim);
        //}// end loop over dim
      //} // end for loop over nodes

      // Update velocity
      //for (int node_gid = 0; node_gid < mesh.num_nodes(); node_gid++) {
        //for (int dim = 0; dim < mesh.num_dim(); dim++){
          node.vel(0,node_gid,dim) = node.vel(num_correction_steps,node_gid,dim);
        }// end loop over dim      
      }// end loop over node_gid 
     
    } // end time integration scope
    
       
    // Track total energy //
    track_rdh(ke,ie, num_correction_steps);
    
        
    
        
    // Increment time //
    TIME += dt;

    if (TIME>=TFINAL) break;
    run_info(cycle);
  }// end loop over time integration cycles

}// end rd_hydro

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
          for (int t_step = 0; t_step < num_correction_steps+1;t_step++){
            for (int dim = 0; dim < mesh.num_dim(); dim++){  
              node.vel(t_step, node_gid, dim) = node.vel(0, node_gid, dim);
            }// end loop over dim
          } // end loop over t_step
          */

