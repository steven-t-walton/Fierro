/* DeC scheme for RD */
#include<iostream>

#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"

using namespace utils;

void rd_hydro(){
  
  real_t ke0 = 0.0;
  real_t ie0 = 0.0;
  
  //test_basis();
  
  for (cycle = 1; cycle <= cycle_stop; cycle++){
        
    std::cout<<" cycle = "<<cycle<<std::endl;

    if (stop_calc == 1) break;

    
    if (cycle == 1){
      track_rdh(ke0, ie0);
      te_0 = ie0 + ke0;
      std::cout << " ke at t0 = " << ke0 << std::endl;
      std::cout << " ie at t0 = " << ie0 << std::endl;
      std::cout << "E_tot at t0 = "<< te_0 << std::endl;
      std::cout << std::endl;
    };// end if
    
    build_corner_normals();
    
    get_timestep();

    dt = fmin(dt, (graphics_time - TIME)+fuzz);

    { // Time integration scope //
      
#pragma omp simd

          
        get_alpha_E();
	get_stress_tensor( 0 );
	get_force_tensor( 0 );

	// Update momentum //
	update_velocity( 0 );
        
	// Update internal energy //
	update_energy( 0 );
      
      // intepolate the velocity with evolved coeffs and save to nodes  //
      interp_vel(1);
      
      // update boundary vel vals //
      boundary_rdh();
     
      // interpolate energy coeffs //
      interp_ie(1);

      // update position //
      update_position(1);   

      get_gauss_pt_jacobian(mesh, ref_elem);
    
      get_gauss_cell_pt_jacobian(mesh, ref_elem);

      get_gauss_patch_pt_jacobian(mesh, ref_elem);

      get_vol_jacobi(mesh, ref_elem);
      
      get_strong_mass();
      
      get_state();
      
      update_coeffs();
      
      for (int gauss_gid = 0; gauss_gid < mesh.num_gauss_pts(); gauss_gid++){
          gauss_properties(gauss_gid);
      }
    }// end time integration scope
       
    // Increment time //
    TIME += dt;

    if (TIME>=TFINAL) break;
    run_info(cycle);
  };// end loop over time integration cycles

}// end rd_hydro

