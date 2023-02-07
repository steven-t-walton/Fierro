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

    /*
    if (cycle == 1){
      track_rdh(ke0, ie0);
      te_0 = ie0 + ke0;
      std::cout << " ke at t0 = " << ke0 << std::endl;
      std::cout << " ie at t0 = " << ie0 << std::endl;
      std::cout << "E_tot at t0 = "<< te_0 << std::endl;
      std::cout << std::endl;
    };// end if
    */
    build_corner_normals();
    
    get_timestep();

    dt = fmin(dt, (graphics_time - TIME)+fuzz);

    { // Time integration scope //
      
#pragma omp simd

          
      get_alpha_E();
      //for (int correction_step = 0; correction_step < num_correction_steps; correction_step++){
	
	//int update = correction_step + 1;
	
	get_stress_tensor(0);// correction_step );
	get_force_tensor(0);// correction_step );

	// Update momentum //
	update_velocity(0);// correction_step );
        
        boundary_rdh();
	
	// Update internal energy //
	update_energy(0);// correction_step );
      
        for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
          for (int dof = 0; dof < ref_elem.num_basis(); dof++){
            int node_lid = elem.vert_node_map(dof);
            int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);
            for (int dim = 0; dim < mesh.num_dim(); dim++){
              node.coords(1, node_gid, dim ) = node.coords( 0, node_gid, dim ) 
		                                 + 0.5*dt*( node.vel(1, node_gid,  dim ) + node.vel( 0, node_gid, dim ) );
            }// end loop over dim
          }// end loop over dof
        }// end loop over elem_gid
      
      //}// end correction step     

      // intepolate the velocity with evolved coeffs and save to nodes  //
      interp_vel(1);//num_correction_steps);
      
      // update boundary vel vals //
      boundary_rdh();
     
      // interpolate energy coeffs //
      interp_ie(1);//num_correction_steps);

      // update position //
      interp_position(1);//num_correction_steps);   

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

