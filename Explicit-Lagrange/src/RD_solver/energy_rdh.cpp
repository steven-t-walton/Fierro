// -----------------------------------------------------------------------------
//  \rho de/dt + \tau : div( u ) = 0 
//
//  currently for
//
//  de/dt + p div( u ) = 0, 
//------------------------------------------------------------------------------
#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"

using namespace utils;

void update_energy( int t_step ){
  
  int update = t_step + 1;

  auto energy_res = CArray <real_t> (mesh.num_elems(), mesh.num_nodes() );
  
  for (int i = 0; i < mesh.num_elems(); i++){
    for (int j = 0; j < mesh.num_nodes(); j++){
      energy_res(i,j) = 0.0;
    }// end loop over j
  }// end loop over i


  for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
  
    for (int t_dof = 0; t_dof < ref_elem.num_dual_basis(); t_dof++){
      int node_lid = ref_elem.dual_vert_node_map(t_dof);
      int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);

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
      
      real_t Me = 0.0;
      for (int basis_id = 0; basis_id < ref_elem.num_dual_basis(); basis_id++){
	//std::cout <<  elem_state.sie_coeffs(t_step, elem_gid, basis_id) - elem_state.sie_coeffs(0, elem_gid, basis_id) << std::endl;
        Me += res_mass(basis_id)*elem_state.sie_coeffs(t_step, elem_gid, basis_id) - res_mass(basis_id)*elem_state.sie_coeffs(0, elem_gid, basis_id);
      }// end loop over basis id
      
      auto F = CArray <real_t> ( num_correction_steps, ref_elem.num_basis(), mesh.num_dim() );
    
        for (int i = 0; i < ref_elem.num_basis(); i++){
          for (int dim = 0; dim < mesh.num_dim(); dim++){
            F(0,i,dim) = 0.0;
            F(1,i,dim) = 0.0;
	  }// end loop over dim
        }// end loop over i
    
      for (int k_dof = 0; k_dof < ref_elem.num_basis(); k_dof++){
        for (int dim = 0; dim < mesh.num_dim(); dim++){
      
          for (int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++){
            int gauss_gid = mesh.gauss_in_elem(elem_gid, gauss_lid);
	    int node_gid = mesh.nodes_in_elem(elem_gid, gauss_lid);

	    real_t J_inv_dot_grad_phi = 0.0;
	    
	    for (int k = 0; k < mesh.num_dim(); k++){
	      J_inv_dot_grad_phi += mesh.gauss_pt_jacobian_inverse(gauss_gid, k, dim)*ref_elem.ref_nodal_gradient(gauss_lid, k_dof, k);
	    }// end loop over k
	    
	    F( 0, k_dof,dim) += J_inv_dot_grad_phi * mat_pt.pressure(gauss_gid) * ref_elem.ref_nodal_dual_basis(gauss_lid, t_dof) 
			            * mesh.gauss_pt_det_j(gauss_gid)
			            * ref_elem.ref_node_g_weights(gauss_lid);
	    F(t_step,k_dof,dim) += J_inv_dot_grad_phi * mat_pt.pressure(gauss_gid) * ref_elem.ref_nodal_dual_basis(gauss_lid, t_dof)
			            * mesh.gauss_pt_det_j(gauss_gid)
			            * ref_elem.ref_node_g_weights(gauss_lid);
          }// end loop over gauss_lid 

      	}// end loop over dim
      }// end loop over k_dof
      
      real_t force_0 = 0.0;
      real_t force_k = 0.0;

      for (int k_dof = 0; k_dof < ref_elem.num_basis(); k_dof++){
        for (int dim = 0; dim < mesh.num_dim(); dim++){
          force_0 += 0.5*F(0,k_dof,dim)*elem_state.vel_coeffs(0 , elem_gid, k_dof, dim) + 0.5*F(0,k_dof,dim)*elem_state.vel_coeffs(t_step , elem_gid, k_dof, dim) ;
          force_k += 0.5*F(t_step,k_dof,dim)*elem_state.vel_coeffs(0, elem_gid, k_dof, dim)
		     + 0.5*F(t_step,k_dof, dim)*elem_state.vel_coeffs(t_step+1,elem_gid,k_dof,dim);	  
 	}// end loop over dim
      }// end loop over k_dof
      
      
      real_t source_int = 0.0;

      for (int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++){
        int gauss_gid = mesh.gauss_in_elem(elem_gid, gauss_lid);
	int node_gid = mesh.nodes_in_elem(elem_gid, gauss_lid);
        
	real_t source_0 = 0.0;	
	real_t source_k = 0.0;
        
	// the usual source term //
        source_0 = 3.141592653589/(4.0*(0.66666667))* ( cos(3.0*3.141592653589 * node.coords(0,node_gid,0))  * cos( 3.141592653589 * node.coords(0,node_gid,1)) - cos( 3.141592653589 * node.coords(0,node_gid,0) ) * cos( 3.0*3.141592653589 * node.coords(0,node_gid, 1) ) ); 
	    
        source_k = 3.141592653589/(4.0*(0.66666667))* ( cos(3.0*3.141592653589 * node.coords(1,node_gid,0))  * cos( 3.141592653589 * node.coords(1,node_gid,1)) - cos( 3.141592653589 * node.coords(1,node_gid,0) ) * cos( 3.0*3.141592653589 * node.coords(1,node_gid, 1) ) ); 
        
	// source from Dobrev et al. //
       // source_0 = ( 0.25*( cos(2.0*3.141592653589*node.coords(0,node_gid,0) ) + cos(2*3.141592653589*node.coords(0, node_gid, 1) ) ) +1 )/(mat_pt.density(gauss_gid)*0.666666667);
 
       // source_k = ( 0.25*( cos(2.0*3.141592653589*node.coords(1,node_gid,0) ) + cos(2*3.141592653589*node.coords(1, node_gid, 1) ) ) +1 )/(mat_pt.density(gauss_gid)*0.666666667);

        source_int += (source_0+source_k)*ref_elem.ref_nodal_dual_basis(gauss_lid,t_dof) * mesh.gauss_pt_det_j(gauss_gid)* ref_elem.ref_node_g_weights(gauss_lid);
      }// end loop over gauss_lid 
     
      //std::cout << source_int << std::endl;
      
      energy_res(elem_gid, node_gid) =  Me/dt + 0.5*(force_0 + force_k) - 0.5*source_int;

      //std::cout << energy_res(elem_gid, node_gid) << std::endl;
      
    }// end loop over t_dof

  }// end loop over elem_gid

  for (int elem_gid = 0; elem_gid <  mesh.num_elems(); elem_gid++){
    for (int t_dof = 0; t_dof < ref_elem.num_dual_basis(); t_dof++){
      
      int node_lid = ref_elem.dual_vert_node_map(t_dof);
      int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);

      real_t lumped_mass = 0.0;
      
      //std::cout << mesh.num_elems_in_node(node_gid) << std::endl;

      for (int elem_node_lid = 0; elem_node_lid < mesh.num_elems_in_node(node_gid); elem_node_lid++){
        int elem_node_gid = mesh.elems_in_node(node_gid, elem_node_lid);
	for (int g_lid = 0; g_lid < mesh.num_gauss_in_elem(); g_lid++){
	  int g_gid = mesh.gauss_in_elem(elem_node_gid, g_lid);
	  lumped_mass += mat_pt.density(g_gid)
		         * ref_elem.ref_nodal_dual_basis(g_lid, t_dof)
			 * ref_elem.ref_node_g_weights(g_lid)
			 * mesh.gauss_pt_det_j(g_gid);
	}// end loop over g_lid
      }// end loop over elem_node_lid
      
      real_t res_sum = 0.0;
      for (int elem_node_lid = 0; elem_node_lid < mesh.num_elems_in_node(node_gid); elem_node_lid++){
        int elem_node_gid = mesh.elems_in_node(node_gid, elem_node_lid);
	res_sum += energy_res(elem_node_gid, node_gid);
      }// end loop over elem_node_lid

      elem_state.sie_coeffs(update, elem_gid, t_dof) = elem_state.sie_coeffs(t_step,elem_gid,t_dof) - (dt/lumped_mass)*res_sum;

    }// end loop over t_dof
  }// end loop over elem_gid

} // end subroutine
