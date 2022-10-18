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
      BV_inv();
      get_control_coeffs();      
      track_rdh(ke, ie, 0);
      te_0 = ie + ke;
      std::cout << "E_tot at t0 = "<< te_0 << std::endl;
      std::cout << " " << std::endl;
    };// end if

    get_timestep();

    dt = fmin(dt, (graphics_time - TIME)+fuzz);
    
    
    { // Time integration scope //
      
      update_coeffs();
#pragma omp simd
       

      // DeC update //
      for (int correction_step = 1; correction_step < num_correction_steps; correction_step++){

	real_t sub_dt = dt/num_correction_steps;
        //real_t sub_time = TIME + sub_dt;          
 
	//std::cout << "calling get_nodal_res" << std::endl;
        get_nodal_res(sub_dt, correction_step);

        //std::cout << "calling lumped_mass " << std::endl;
        lumped_mass();// used in momentum
       
      	// updates velocity at vel(t^{n,m}) //
	get_momentum_rd(correction_step);
        //std::cout << " updated momentum " << std::endl;
       
	// update velocity at boundary appropriately //
       // (inside momentum)	boundary_rdh(correction_step+1);

        //std::cout << "calling prediction step" << std::endl;
       // prediction_step(sub_dt, correction_step);

	// energy update //
//        get_energy_rdh( sub_dt );

      }//end correction steps

      get_position_rdh();   

      //std::cout << " calling get state " << std::endl;
      get_state( cycle );

      get_stress(); // assign values to stress
      //std::cout << "called get stress" << std::endl;
	

    }; // end time integration scope
       
    // Increment time //
    TIME += dt;

    if (TIME>=TFINAL) break;
    run_info(cycle);
  };// end loop over time integration cycles


  // final E_tot //
  track_rdh( ke, ie, num_correction_steps-1);
  std::cout << "E_tot final is "<< ke+ie << std::endl;

}// end rd_hydro


