/* DeC scheme for RD */


#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"

using namespace utils;

void rd_hydro(){
  
  for (cycle = 1; cycle <= cycle_stop; cycle++){
    if (stop_calc == 1) break;
    
    track_rdh(ke, ie, 0);
    
    if (cycle == 1) te_0 = ie + ke;
    
    get_timestep();

    dt = fmin(dt, (graphics_time - TIME)+fuzz);
    
    { // Time integration scope //
 
      // Save t^n <--| t^{r+1,m} //   
  
      //for(int gauss_gid = 0; gauss_gid < mesh.num_gauss_pts(); gauss_gid++){
      //  for(int i=0; i<mesh.num_dim(); i++) node.vel(0, gauss_gid, i) = node.vel(num_correction_steps, gauss_gid, i);
     // } // end loop over gauss_gid 
      
      rk_init();

      // DeC update //
      for (int pred_step = 0; pred_step < num_prediction_steps; pred_step++){
        std::cout<<" pred_step = "<< pred_step << std::endl;
        for (int correction_step = 0; correction_step <= num_correction_steps; correction_step++){
          std::cout<< " correction step = "<< correction_step << std::endl;          
	  real_t sub_dt = (correction_step+1)*dt/num_correction_steps;
          real_t sub_time = TIME + sub_dt;          
          
	  get_stress(correction_step); // assign values to stress
          std::cout << "called get stress" << std::endl;
	  
          get_cell_mass(); // lumped mass in a cell. used in lo_res
          std::cout << " called get cell mass" << std::endl;
          
          std::cout << "calling get_lo_res" << std::endl;
          get_lo_res(dt, correction_step, TIME);// low order approx
          

          std::cout << "calling lumped_mass " << std::endl;
	  lumped_mass();// used in momentum
          
          std::cout << "calling prediction step" << std::endl;
          prediction_step(sub_dt, pred_step);// low order approx
          
          std::cout << " updating momentum " << std::endl;
          get_momentum_rd(pred_step, correction_step);// high order correction
          
          std::cout << " calling get state " << std::endl;
          get_state();
        }//end correction steps
      }//end prediction steps
    } // end time integration scope
    
    // Track total energy //
    track_rdh(ke,ie, num_correction_steps);
    
    //Update position
    for (int node_gid = 0; node_gid < mesh.num_nodes(); node_gid++) {
    // create view of the nodal velocities
      auto vel = ViewCArray <real_t> (&node.vel(num_correction_steps, node_gid, 0), num_dim);
      for (int dim = 0; dim < 3; dim++){
        node.coords(1, node_gid, dim) = node.coords(0, node_gid, dim) +  dt*(vel(dim));
        mesh.node_coords(node_gid, dim) = node.coords(1, node_gid, dim);
      }
    } // end for loop over nodes

    // Increment time //
    TIME += dt;

    if (TIME>=TFINAL) break;

  }// end loop over time integration cycles

}// end rd_hydro
