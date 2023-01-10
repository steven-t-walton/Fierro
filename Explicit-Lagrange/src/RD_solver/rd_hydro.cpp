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
      //interp_vel(0);
      //interp_ie(0);
      track_rdh(ke0, ie0);
      te_0 = ie0 + ke0;
      std::cout << " ke at t0 = " << ke0 << std::endl;
      std::cout << " ie at t0 = " << ie0 << std::endl;
      std::cout << "E_tot at t0 = "<< te_0 << std::endl;
      std::cout << std::endl;
    };// end if
    

    get_timestep();

    dt = fmin(dt, (graphics_time - TIME)+fuzz);

    { // Time integration scope //
      
      //update_coeffs();
#pragma omp simd

      // DeC update //
      for (int correction_step = 0; correction_step < num_correction_steps; correction_step++){

	//real_t sub_dt = correction_step*dt/num_correction_steps;
        //real_t sub_time = TIME + sub_dt;          
 
        update_velocity( correction_step );
        
	// Update internal energy //
	update_energy( correction_step );
      
	// Update position coefficients //
        for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
	  for (int dof = 0; dof < ref_elem.num_basis(); dof++){
	    for (int dim = 0; dim < mesh.num_dim(); dim++){
	      elem_state.pos_coeffs(correction_step+1, elem_gid, dof, dim) = elem_state.pos_coeffs(0, elem_gid,dof, dim)
		             + 0.5*dt*(elem_state.vel_coeffs(correction_step, elem_gid, dof, dim) + elem_state.vel_coeffs(0,elem_gid, dof, dim) );
	    }
	  }
	}
        
      }//end correction steps
      
      // intepolate the velocity with evolved coeffs and save to nodes  //
      interp_vel(num_correction_steps);
      
      // update boundary vel vals //
      boundary_rdh();
     
      // interpolate energy coeffs //
      interp_ie(num_correction_steps);

      // update position //
      update_position();   

   
      get_gauss_pt_jacobian(mesh, ref_elem);
    
      get_gauss_cell_pt_jacobian(mesh, ref_elem);

      get_gauss_patch_pt_jacobian(mesh, ref_elem);

      get_vol_jacobi(mesh, ref_elem);
      
      get_state();

      get_stress(); 

      //for(int gauss_gid=0; gauss_gid<mesh.num_gauss_pts(); gauss_gid++){
      //  gauss_properties(gauss_gid);
      //}// end loop over gauss_gid

      update_coeffs();
       
    }// end time integration scope
       
    // Increment time //
    TIME += dt;

    if (TIME>=TFINAL) break;
    run_info(cycle);
  };// end loop over time integration cycles


  // final E_tot //
  track_rdh( ke, ie );

  std::cout << " ke at t0 = " << ke0 << std::endl;
  std::cout << " ie at t0 = " << ie0 << std::endl;
  std::cout << " ke at t_final = " << ke << std::endl;
  std::cout << " ie at t_final = " << ie << std::endl;
  std::cout << "E_tot final is "<< ke+ie << std::endl;

}// end rd_hydro



	//get_total_res();

	//get_betas();

	//get_limited_res();
        
      	// updates velocity coeffs //
	//get_momentum_rd( correction_step );
	
	// energy update //
        //get_energy_rdh( sub_dt );

