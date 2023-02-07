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

 
  FOR_ALL( elem_gid , 0, mesh.num_elems(),{
    for (int t_dof = 0; t_dof < ref_elem.num_dual_basis(); t_dof++){
      
      int node_lid = ref_elem.dual_vert_node_map(t_dof);
      int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);
      int gauss_gid = mesh.gauss_in_elem(elem_gid, node_lid);
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
      if (lumped_mass <= 0.0){
        std::cout << " thermodynamic lumped mass is negative " << lumped_mass << std::endl;
      }      

      real_t sum_res = 0.0;
      get_thermodynamic_L2( t_step, node_gid, sum_res); 
      // L^1(e^{k+1}) = L^1(e^k) - L^2(e^k) //
      mat_pt.sie(update, gauss_gid) = mat_pt.sie(t_step, gauss_gid) - (dt/lumped_mass)*sum_res;
/*
      if (elem_state.sie_coeffs(update, elem_gid, t_dof) <= 0.0){
        std::cout << " sie is negative in elem "<< elem_gid << " and dof "<< t_dof << "with value "<< elem_state.sie_coeffs(update, elem_gid, t_dof) << std::endl; 
      }
*/
    }// end loop over t_dof
  });// end loop over elem_gid
  
} // end subroutine

