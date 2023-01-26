
#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"
#include<iostream>

using namespace utils;

void get_force_tensor(int t_step){

  //for (int t_step = 0; t_step < correction_storage-1; t_step++){
    for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
      for (int k_dof = 0; k_dof < ref_elem.num_basis(); k_dof++){
	for (int t_dof = 0; t_dof < ref_elem.num_dual_basis(); t_dof++){
          for (int dim = 0; dim < mesh.num_dim(); dim++){
            for (int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++){
	       elem_state.force_tensor(t_step, elem_gid, k_dof, t_dof, dim) = 0.0;  
	    }
	  }
	}
      }
    }
  //}


  //for (int t_step = 0; t_step < correction_storage-1; t_step++){
    for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
      for (int k_dof = 0; k_dof < ref_elem.num_basis(); k_dof++){
	for (int t_dof = 0; t_dof < ref_elem.num_dual_basis(); t_dof++){
          for (int dim = 0; dim < mesh.num_dim(); dim++){
            for (int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++){
              
	      //auto J_inv_dot_grad_phi = CArray <real_t> (mesh.num_dim());
              //real_t sigma_dot_grad_phi = 0.0;
              
	      int gauss_gid = mesh.gauss_in_elem(elem_gid, gauss_lid);
	      /*
	      for (int l =0; l < mesh.num_dim(); l++){
	        for (int k = 0; k < mesh.num_dim(); k++){
	          J_inv_dot_grad_phi(l) += mesh.gauss_pt_jacobian_inverse(gauss_gid, k, l)*ref_elem.ref_nodal_gradient(gauss_lid, k_dof, k);
	        }
	      }
               
	      for (int l = 0; l < mesh.num_dim(); l++){
	        sigma_dot_grad_phi += elem_state.stress_tensor(t_step,gauss_gid, dim, l)*J_inv_dot_grad_phi(l); 
	      }
	      */
              //std::cout << " sigma.grad varphi = " << sigma_dot_grad_phi << std::endl;
	      for (int k = 0; k < mesh.num_dim(); k++){
	        for (int l = 0; l < mesh.num_dim(); l++){
	          elem_state.force_tensor(t_step, elem_gid, k_dof, t_dof, dim) += ref_elem.ref_nodal_gradient(gauss_lid, k_dof,k)
		                                                              * mesh.gauss_pt_jacobian_inverse(gauss_gid, k, l)
									      * elem_state.stress_tensor(t_step, gauss_gid, dim, l)
								              * ref_elem.ref_nodal_dual_basis(gauss_lid, t_dof) 
								              * mesh.gauss_pt_det_j(gauss_gid) 
								              * ref_elem.ref_node_g_weights(gauss_lid);
		}// end loop over l
              }// end loop over k
            }// end loop over gauss_lid

      	  }// end loop over dim
	}// end loop over t_dof
      }// end loop over k_dof
    }// end loop over elem_gid
 // }// end loop over t_step

}// end get_force_tensor

