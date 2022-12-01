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
  
  for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
  
    // --- Assemble Force Matrix, a la Dobrev et al.  --- //
    // F = \int_{\Omega_t} (\sigma:\nabla \varphi)\psi dx //
    // \varphi are kinematic basis function //
    // \psi are thermodynamic basis function //
    // The Galerkin residual is then R(e) = M.de/dt + F^T.v // 
    
    auto F = CArray <real_t> ( num_correction_steps, ref_elem.num_basis(), ref_elem.num_dual_basis() );
    
    for (int t_step = 0; t_step < num_correction_steps; t_step++){
      for (int i = 0; i < ref_elem.num_basis(); i++){
        for (int j = 0; j < ref_elem.num_dual_basis(); j++){
          for (int dim = 0; dim < mesh.num_dim(); dim++){
            F(t_step,i,j,dim) = 0.0;
	  }// end loop over dim
        }// end loop over j
      }// end loop over i
    }// end loop over t_step
    
    for (int k_dof = 0; k_dof < ref_elem.num_basis(); k_dof++){
      for (int t_dof = 0; t_dof < ref_elem.num_dual_basis(); t_dof++){
        for (int dim = 0; dim < mesh.num_dim(); dim++){
        
	auto F0 = ViewCArray <real_t> (&F(0,k_dof, t_dof, dim), ref_elem.num_basis(), ref_elem.num_dual_basis(), dim);
        auto Fk = ViewCArray <real_t> (&F(1,k_dof, t_dof, dim), ref_elem.num_basis(), ref_elem.num_dual_basis(), dim );
      
        for (int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++){
          int gauss_gid = mesh.gauss_in_elem(elem_gid, gauss_lid);
	  real_t J_inv_dot_grad_phi = 0.0;
	  for (int k = 0; k < mesh.num_dim(); k++){
	    J_inv_dot_grad_phi += mesh.gauss_pt_jacobian_inverse(gauss_gid, k, dim)*ref_elem.ref_nodal_gradient(gauss_lid, k_dof, k);
	  }// end loop over k
	  F0(k_dof, t_dof,dim) -= J_inv_dot_grad_phi * mat_pt.pressure(gauss_gid)// <-- stress(0, gauss_gid, dim, dim)
	            	       * ref_elem.ref_nodal_dual_basis(gauss_lid, t_dof)
			       * mesh.gauss_pt_det_j(gauss_gid)
			       * ref_elem.ref_node_g_weights(gauss_lid);
	  Fk(k_dof, t_dof,dim) -= J_inv_dot_grad_phi * mat_pt.pressure(gauss_gid)// <-- stress(t_step, gauss_gid,dim,dim)
	            	       * ref_elem.ref_nodal_dual_basis(gauss_lid, t_dof)
			       * mesh.gauss_pt_det_j(gauss_gid)
			       * ref_elem.ref_node_g_weights(gauss_lid);
        }// end loop over gauss_lid 

        }// end loop over dim
      }// end loop over t_dof
    }// end loop over k_dof


  }// end loop over elem_gid


} // end subroutine
