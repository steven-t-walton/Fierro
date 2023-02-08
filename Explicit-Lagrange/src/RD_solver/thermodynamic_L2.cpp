#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"
#include<iostream>

using namespace utils;

void get_thermodynamic_L2( int t_step, int dof_gid, real_t& sum_res ){
  
  int num_elems_in_dof = mesh.num_elems_in_node(dof_gid);
  real_t inv_num_dual_basis = 1.0/ref_elem.num_dual_basis();
  
  for (int elem_lid = 0; elem_lid < num_elems_in_dof; elem_lid++){
    int elem_gid = mesh.elems_in_node(dof_gid, elem_lid);

    auto energy_res = CArray <real_t> ( ref_elem.num_dual_basis() );
  
    for (int j = 0; j < ref_elem.num_dual_basis(); j++){
      energy_res(j) = 0.0;
    }// end loop over j
  
    for (int t_dof = 0; t_dof < ref_elem.num_dual_basis(); t_dof++){
      //int node_lid = ref_elem.dual_vert_node_map(t_dof);
      //int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);

      auto res_mass = CArray <real_t> (ref_elem.num_dual_basis());
      for (int i = 0; i < ref_elem.num_dual_basis(); i++) res_mass(i) = 0.0;

      for (int basis_id = 0; basis_id < ref_elem.num_dual_basis(); basis_id++){
        for (int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++){
	  int gauss_gid = mesh.gauss_in_elem(elem_gid, gauss_lid);
	  res_mass(basis_id) += ref_elem.ref_nodal_dual_basis(gauss_lid, t_dof)
		                * ref_elem.ref_nodal_dual_basis(gauss_lid, basis_id)
				* mat_pt.density(gauss_gid)
				* ref_elem.ref_node_g_weights(gauss_lid)
				* mesh.gauss_pt_det_j(gauss_gid);
	}// end loop over gauss_lid
      }// end loop over basis_id
      
      real_t M_dot_e = 0.0;
      for (int basis_id = 0; basis_id < ref_elem.num_dual_basis(); basis_id++){
        M_dot_e += res_mass(basis_id)*(elem_state.sie_coeffs(t_step, elem_gid, basis_id) - elem_state.sie_coeffs(0, elem_gid, basis_id));
      }// end loop over basis id
      
      real_t force = 0.0;
      
      for (int dim = 0; dim < mesh.num_dim(); dim++){
        for (int k_dof = 0; k_dof < ref_elem.num_basis(); k_dof++){
          force += 0.5*(0.5*elem_state.force_tensor(t_step, elem_gid, k_dof, t_dof, dim)
		      * (elem_state.vel_coeffs(0, elem_gid, k_dof, dim) + elem_state.vel_coeffs(t_step,elem_gid,k_dof,dim)) 
		       + elem_state.force_tensor(0, elem_gid, k_dof, t_dof, dim)*elem_state.vel_coeffs(0, elem_gid, k_dof, dim));
        }// end loop over k_dof
      }// end loop over dim
      
      //--- Artificial Viscosity ---//
      real_t sie_bar = 0.0;
      real_t sie_bar0 = 0.0;

      real_t Q = 0.0;
      
       
      // Compute sie_bar //
      for (int dof = 0; dof < ref_elem.num_dual_basis(); dof++){
        sie_bar += elem_state.sie_coeffs(t_step, elem_gid, dof)*inv_num_dual_basis; 
        sie_bar0 += elem_state.sie_coeffs(0, elem_gid, dof)*inv_num_dual_basis; 
      }
  
      // Fill Q //
      Q = 0.5*elem_state.alpha_E(elem_gid)*(elem_state.sie_coeffs(t_step, elem_gid, t_dof) - sie_bar)
	  + 0.5*elem_state.alpha_E(elem_gid)*(elem_state.sie_coeffs(0, elem_gid, t_dof) - sie_bar0);
  
      //--- end Artificial Viscosity ---//

      real_t source_int = 0.0;
      
      for (int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++){
        int gauss_gid = mesh.gauss_in_elem(elem_gid, gauss_lid);
	int node_gid = mesh.nodes_in_elem(elem_gid, gauss_lid);
        
	real_t source = 0.0;	
        
	// the usual source term //
        source = 0.5*(3.141592653589/(4.0*(0.66666667))* ( cos(3.0*3.141592653589 * node.coords(0,node_gid,0))  * cos( 3.141592653589 * node.coords(0,node_gid,1)) - cos( 3.141592653589 * node.coords(0,node_gid,0) ) * cos( 3.0*3.141592653589 * node.coords(0,node_gid, 1) ) )
		       	+ 3.141592653589/(4.0*(0.66666667))* ( cos(3.0*3.141592653589 * node.coords(1,node_gid,0))  * cos( 3.141592653589 * node.coords(1,node_gid,1)) - cos( 3.141592653589 * node.coords(1,node_gid,0) ) * cos( 3.0*3.141592653589 * node.coords(1,node_gid, 1) ) ) ); 
	    
        source_int += (source)*ref_elem.ref_nodal_dual_basis(gauss_lid,t_dof) * mesh.gauss_pt_det_j(gauss_gid)* ref_elem.ref_node_g_weights(gauss_lid);
      }// end loop over gauss_lid 
      //std::cout << source_int << std::endl;
      
      real_t sub_factor = 1.0/((real_t)num_correction_steps - (real_t)t_step);
      energy_res(t_dof) =  M_dot_e + sub_factor*dt*(Q - force - source_int);

    }// end loop over t_dof

 
 /* 
    for (int t_dof = 0; t_dof< ref_elem.num_dual_basis(); t_dof++){
      std::cout << energy_res(t_dof) << std::endl;
    }
  } 
 */ 
    real_t total_energy_res = 0.0;

    for (int dof = 0; dof < ref_elem.num_dual_basis(); dof++){
      total_energy_res += energy_res( dof);
    }

    real_t inv_total_res = 1.0/total_energy_res;

    for (int dof = 0; dof < ref_elem.num_dual_basis(); dof++){
      
      real_t numerator = 0.0;
      numerator = std::max(0.0,( energy_res(dof)*inv_total_res ) );
      real_t denom = 0.0;
      for (int dof_id = 0; dof_id < ref_elem.num_dual_basis(); dof_id++){
        denom += std::max( 0.0, ( energy_res(dof_id)*inv_total_res ) ); 
      }

      int node_lid = ref_elem.dual_vert_node_map(dof);
      int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);

      if (dof_gid == node_gid){
        sum_res += energy_res(dof);
	//sum_res += numerator/denom*total_res(dof);
      }

    }
  
  }// end loop over elem_lid
  
} // end subroutine


  /*
    std::cout << " total thermo residual in elem " << elem_gid << " is " << total_energy_res << std::endl;
  */
/*
    auto betas = CArray <real_t> ( ref_elem.num_dual_basis());
    for(int dof = 0; dof < ref_elem.num_dual_basis(); dof++){
       betas(dof) = 0.0; 
    }
 
    for(int dof = 0; dof < ref_elem.num_dual_basis(); dof++){
      real_t numerator = 0.0;
      numerator = std::max(0.0,( energy_res(dof)/total_energy_res ) );
      real_t denom = 0.0;
      for (int dof_id = 0; dof_id < ref_elem.num_dual_basis(); dof_id++){
        denom += std::max( 0.0, ( energy_res(dof_id)/total_energy_res ) ); 
      }
      betas( dof) = numerator/denom; 
    }
 */
  /* 
  // check that betas sum to 1 //
    real_t sum = 0.0;
    for (int vertex = 0; vertex < ref_elem.num_dual_basis(); vertex++){
      sum += betas(vertex);
    }// end loop over vertex
    std::cout << " sum of betas in elem "<< elem_gid << " is " << sum << std::endl;
  */
/*
    auto limited_energy_res = CArray <real_t> (ref_elem.num_dual_basis());
  
    for(int dof=0; dof < ref_elem.num_dual_basis(); dof++){
       limited_energy_res(dof) = 0.0; 
    }
    
    for(int dof=0; dof < ref_elem.num_dual_basis(); dof++){

      limited_energy_res(dof) = betas(dof)*total_energy_res;

    }
*/
/*
    for(int dof=0; dof < ref_elem.num_dual_basis(); dof++){
      std::cout << " limited thermo residual in elem " << elem_gid << " at node " << dof << " is " << limited_energy_res(dof) << std::endl;
    }
*/
