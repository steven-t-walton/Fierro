/* DeC scheme for RD */


#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"

using namespace utils;

void rd_hydro(){
  
    real_t ke0 = 0.0;
    real_t ie0 = 0.0;

  for (cycle = 1; cycle <= cycle_stop; cycle++){
        
    std::cout<<"cycle = "<<cycle<<std::endl;

    if (stop_calc == 1) break;
    
    if (cycle == 1){
      BV_inv();
      get_control_coeffs();
//      interp_vel(0);      
      track_rdh(ke0, ie0, 0);
      te_0 = ie0 + ke0;
      std::cout << " ke at t0 = " << ke0 << std::endl;
      std::cout << " ie at t0 = " << ie0 << std::endl;
      std::cout << "E_tot at t0 = "<< te_0 << std::endl;
      std::cout << std::endl;
    };// end if

    get_timestep();

    dt = fmin(dt, (graphics_time - TIME)+fuzz);
    
    
    { // Time integration scope //
      
      update_coeffs();
       

#pragma omp simd
      // DeC update //
      for (int correction_step = 1; correction_step <= num_correction_steps; correction_step++){

	real_t sub_dt = dt/num_correction_steps;
        //real_t sub_time = TIME + sub_dt;          
 
        //std::cout << "calling lumped_mass " << std::endl;
        lumped_mass();// used in momentum
       
	//std::cout << "calling get_nodal_res" << std::endl;
        get_nodal_res(sub_dt, correction_step);

      	// updates velocity at vel(t^{n,m}) //
	get_momentum_rd(correction_step);
        //std::cout << " updated momentum " << std::endl;
        
        interp_vel(correction_step);
       
        get_position_rdh(correction_step);   

//        interp_pos(correction_step);

        //std::cout << "Calculating Jacobian at gauss points" << std::endl;
        get_gauss_pt_jacobian(mesh, ref_elem);
    
        //std::cout << "Calculating Jacobian at gauss points in cell" << std::endl;
        get_gauss_cell_pt_jacobian(mesh, ref_elem);

        //std::cout << "Before volume from Jacobian"  << std::endl;
        get_vol_jacobi(mesh, ref_elem);

	// energy update //
//        get_energy_rdh( sub_dt );

      }//end correction steps

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
  track_rdh( ke, ie, num_correction_steps);

  std::cout << " ke at t0 = " << ke0 << std::endl;
  std::cout << " ie at t0 = " << ie0 << std::endl;
  std::cout << " ke at t_final = " << ke << std::endl;
  std::cout << " ie at t_final = " << ie << std::endl;
  std::cout << "E_tot final is "<< ke+ie << std::endl;

}// end rd_hydro


