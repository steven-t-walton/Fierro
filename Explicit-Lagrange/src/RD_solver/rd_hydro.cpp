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

    build_corner_normals();
    
    get_timestep();

    dt = fmin(dt, (graphics_time - TIME)+fuzz);

    { // Time integration scope //
      
#pragma omp simd

          
      get_alpha_E();
      for (int k = 0; k < num_correction_steps; k++){
        
	get_stress_tensor( k );
        get_force_tensor( k );

	// Update momentum //
	update_velocity( k );
        
        boundary_rdh();
	
	// Update internal energy //
	update_energy( k );

	real_t vel_cons = 0.0;
	real_t alpha = (real_t)k + 1.0;
	for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
          for (int dof = 0; dof < ref_elem.num_basis(); dof++){
            int node_lid = elem.vert_node_map(dof);
            int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);
            for (int dim = 0; dim < mesh.num_dim(); dim++){
	      vel_cons = 0.5*(node.vel(k, node_gid, dim) + node.vel(0, node_gid, dim));
              node.coords(k+1, node_gid, dim ) = node.coords(0, node_gid,dim) + 0.5*alpha*dt*vel_cons;
	    }// end loop over dim
          }// end loop over dof
        }// end loop over elem_gid
      
      }// end correction step     

      // intepolate the velocity with evolved coeffs and save to nodes  //
      interp_vel();
      
      // interpolate energy coeffs //
      interp_ie();

      // update position //
      interp_position();   

      get_gauss_pt_jacobian(mesh, ref_elem);
    
      get_gauss_cell_pt_jacobian(mesh, ref_elem);

      get_gauss_patch_pt_jacobian(mesh, ref_elem);

      get_vol_jacobi(mesh, ref_elem);
      
      get_strong_mass();
      
      get_state();
      
      update_tn();
      
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

