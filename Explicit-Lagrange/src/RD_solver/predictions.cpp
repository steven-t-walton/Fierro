// Computes forward eurler prediction //


#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"

using namespace utils;


void prediction(int prediction_step, real_t xi_dt){
  
  for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
  
    for (int k_dof = 0; k_dof < ref_elem.num_basis(); k_dof++){
      int node_lid = ref_elem.vert_node_map(k_dof);
      int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);
      int num_elems_in_vert = mesh.num_elems_in_node(node_gid);	    
      
      // Compute lumped mass
      real_t lumped_mass = 0.0;
      
      for (int elems_in_vert = 0; elems_in_vert < num_elems_in_vert; elems_in_vert++){
        int elems_in_vert_gid = mesh.elems_in_node(node_gid, elems_in_vert);
        for (int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++){
          int gauss_gid = mesh.gauss_in_elem(elems_in_vert_gid, gauss_lid);
          lumped_mass +=  ref_elem.ref_node_g_weights(gauss_lid)
 		          * mat_pt.density(gauss_gid)
		          * mesh.gauss_pt_det_j(gauss_gid)
	                  * ref_elem.ref_nodal_basis(gauss_lid, k_dof);
        }// end loop over gauss_lid
      }// end loop over elems_in_node_lid
      for (int dim = 0; dim < mesh.num_dim(); dim++){ 
        elem_state.vel_coeffs(prediction_step, elem_gid, k_dof, dim) = elem_state.vel_coeffs(0, elem_gid, k_dof, dim) 
		                                                       + xi_dt*elem_state.K_res_sum(0,node_gid,dim)/lumped_mass;
        //std::cout << "vel_coeffs at pred step " << elem_state.vel_coeffs(prediction_step, elem_gid, k_dof,dim) << std::endl;	
	
      }// end loop over dim
    }// end loop over k_dof
  }//end loop over elem_gid
   
  for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
  
    for (int t_dof = 0; t_dof < ref_elem.num_dual_basis(); t_dof++){
      int dual_node_lid = ref_elem.dual_vert_node_map(t_dof);
      int dual_node_gid = mesh.nodes_in_elem(elem_gid, dual_node_lid);
      int num_elems_in_dual_vert = mesh.num_elems_in_node(dual_node_gid);	    
      real_t dual_lumped_mass = 0.0;

      for (int elem_node_lid = 0; elem_node_lid < num_elems_in_dual_vert; elem_node_lid++){
        int elem_node_gid = mesh.elems_in_node(dual_node_gid, elem_node_lid);
        for (int g_lid = 0; g_lid < mesh.num_gauss_in_elem(); g_lid++){
          int g_gid = mesh.gauss_in_elem(elem_node_gid, g_lid);
          dual_lumped_mass += mat_pt.density(g_gid)
    	  	             * ref_elem.ref_nodal_dual_basis(g_lid, t_dof)
			     * ref_elem.ref_node_g_weights(g_lid)
			     * mesh.gauss_pt_det_j(g_gid);
        }// end loop over g_lid
      }// end loop over elem_node_lid
      

      elem_state.sie_coeffs(prediction_step, elem_gid, t_dof) = elem_state.sie_coeffs(0,elem_gid,t_dof) + (xi_dt/dual_lumped_mass)*elem_state.T_res_sum(0,dual_node_gid);
    
      //std::cout << "sie_coeffs at pred step " << elem_state.sie_coeffs(prediction_step, elem_gid, t_dof) << std::endl;	
    
    }//end loop over t_dof
  }// end loop over elem_gid
}