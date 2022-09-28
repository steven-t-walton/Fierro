/* DeC scheme for RD */


#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"

using namespace utils;

void rd_hydro(){
  
  for (cycle = 1; cycle <= cycle_stop; cycle++){
    
    std::cout<<"cycle = "<<cycle<<std::endl;

    if (stop_calc == 1) break;
    if (cycle == 1){
      for (int elem_gid = 0 ; elem_gid < mesh.num_elems(); elem_gid++){
        for (int node_lid = 0; node_lid < mesh.num_nodes_in_elem(); node_lid++){
          int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);
          for (int t_step = 0; t_step < num_correction_steps+1;t_step++){
            for (int dim = 0; dim < mesh.num_dim(); dim++){  
              node.vel(t_step, node_gid, dim) = node.vel(0, node_gid, dim);
            }// end loop over dim
          } // end loop over t_step
        }// end loop over node_lid
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
        
        //std::cout << " calling get state " << std::endl;
        get_state();
      }//end correction steps
      // Update velocity
      for (int node_gid = 0; node_gid < mesh.num_nodes(); node_gid++) {
        for (int dim = 0; dim < mesh.num_dim(); dim++){
          node.vel(0,node_gid,dim) = node.vel(num_correction_steps,node_gid,dim);
        }// end loop over dim       
      }// end loop over node_gid

     
     
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
      }// end loop over dim
    } // end for loop over nodes
    
        
    // Increment time //
    TIME += dt;

    if (TIME>=TFINAL) break;
    run_info(cycle);
  }// end loop over time integration cycles

}// end rd_hydro
