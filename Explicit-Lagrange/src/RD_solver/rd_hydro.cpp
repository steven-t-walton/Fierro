/* DeC scheme for RD */


#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"

using namespace utils;

void rd_hydro(){
  
  int update = num_correction_steps-1;

  for (cycle = 1; cycle <= cycle_stop; cycle++){
        
    //std::cout<<"cycle = "<<cycle<<std::endl;

    if (stop_calc == 1) break;
    
        
    if (cycle == 1){
      BV_inv();          
      track_rdh(ke, ie, 0);
      te_0 = ie + ke;
      std::cout << "E_tot at t0 = "<< te_0 << std::endl;
      std::cout << " " << std::endl;
    };// end if

    get_timestep();

    dt = fmin(dt, (graphics_time - TIME)+fuzz);
    
    
    { // Time integration scope //

     for (int node_gid = 0; node_gid < mesh.num_nodes(); node_gid++) {
        
       for (int dim = 0; dim < 3; dim++){
         node.coords(0, node_gid, dim) = node.coords(1, node_gid, dim);
         mesh.node_coords(node_gid, dim) = node.coords(1, node_gid, dim);

         //std::cout<< "pos at time "<< TIME+dt <<" cycle "<< cycle << " node " << node_gid << " and dim " << dim <<" is "<< node.coords(0,node_gid,dim) << std::endl;
         // Update velocity
         node.vel( 0, node_gid, dim) = node.vel( update, node_gid, dim);
         if (cycle == 1){         
          real_t BV_coeffs_a[mesh.num_elems()*elem.num_basis()*mesh.num_dim()];
          auto BV_coeffs = ViewCArray <real_t> (&BV_coeffs_a[0], mesh.num_elems(), elem.num_basis(), mesh.num_dim());
#pragma omp simd
          // initialize control coeffs to zero //
          for( int k = 0; k < mesh.num_dim(); k ++){
            for (int j = 0; j < elem.num_basis(); j++){
              for (int i = 0; i < mesh.num_elems(); i++){
                 BV_coeffs(i,j,k) =  0.0;
              }// end loop over i
            }// end loop over j
          }// end loop over k

#pragma omp simd      
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
		
		node.vel(0, node_gid_for_control_coeff, dim) = BV_coeffs(elem_gid, basis_id, dim);
		
              }// end loop over basis id
            }// end loop over dim
          }// end loop over elem_gid 
          }// end if
	 else if ( 1 < cycle ){
           for( int dim = 0; dim < mesh.num_dim(); dim++){
             for (int node_gid = 0; node_gid < mesh.num_nodes(); node_gid ++){
	       for (int t_step = 1 ; t_step < num_correction_steps; t_step++){
                 node.vel(t_step, node_gid, dim) = node.vel(0,node_gid,dim) ;
               }
	     }
	   }
	 };
//	 std::cout << " " << std::endl;
//         std::cout<< "vel at time "<< TIME+dt <<" cycle "<< cycle << " node " << node_gid << " and dim " << dim <<" is "<< node.vel(0,node_gid,dim) << std::endl;
//        std::cout << " cycle = " << cycle << std::endl;
        // update velocity at boundary appropriately //
       	boundary_rdh(0);

       }// end loop over dim      
     }// end loop over node_gid 

#pragma omp simd     
     for (int cell_gid = 0; cell_gid < mesh.num_cells(); cell_gid++){

       cell_state.ie( 0, cell_gid ) = cell_state.ie( update, cell_gid );

       //cell_state.total_energy( 0, cell_gid ) = cell_state.total_energy( 1, cell_gid );
     }
    

#pragma omp simd
       
      // DeC update //
      for (int correction_step = 1; correction_step < num_correction_steps; correction_step++){

	real_t sub_dt = dt/num_correction_steps;
        //real_t sub_time = TIME + sub_dt;          
 
	get_stress(correction_step); // assign values to stress
        //std::cout << "called get stress" << std::endl;

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

        //std::cout << "Calculating Jacobian at gauss points" << std::endl;
        get_gauss_pt_jacobian(mesh, ref_elem);
    
       //std::cout << "Calculating Jacobian at gauss points in cell" << std::endl;
        get_gauss_cell_pt_jacobian(mesh, ref_elem);

       //std::cout << "Before volume from Jacobian"  << std::endl;
        get_vol_jacobi(mesh, ref_elem);      

        //std::cout << " calling get state " << std::endl;
	get_state(cycle, correction_step);

      }//end correction steps

    
    }; // end time integration scope
       
    // Increment time //
    TIME += dt;

    if (TIME>=TFINAL) break;
    run_info(cycle);
  };// end loop over time integration cycles


  // final E_tot //
  track_rdh( ke, ie, update);
  std::cout << "E_tot final is "<< ke+ie << std::endl;

}// end rd_hydro


