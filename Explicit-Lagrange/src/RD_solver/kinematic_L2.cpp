#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"

using namespace utils;

void get_kinematic_L2(int t_step, int dof_gid, int dim, real_t& sum_res){

  int num_basis = ref_elem.num_basis();
  int num_dim = mesh.num_dim();
  
  int num_elems_in_dof = mesh.num_elems_in_node(dof_gid);

  for(int elem_lid = 0; elem_lid < num_elems_in_dof; elem_lid++){
    
    int elem_gid = mesh.elems_in_node(dof_gid, elem_lid);
    
    auto galerkin_res = CArray <real_t> (num_basis);
    for (int i = 0; i < num_basis; i++) galerkin_res(i) = 0.0;

    // Compute residual //
    for (int vertex = 0; vertex < num_basis; vertex++){
      
      int node_lid = elem.vert_node_map(vertex);
      int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);
      int g_gid = mesh.gauss_in_elem(elem_gid, node_lid);
	
      // View for mass_vector //
      CArray <real_t> res_mass(  num_basis );
      for (int i = 0; i < num_basis; i++) res_mass(i) = 0.0;

      // Create view for \sum_{q} M_{pq} \delta u^r_q //      
      real_t Mv = 0.0;

      //--- Mass Vector ---//
        
      for (int index = 0; index < num_basis; index++){
        for (int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++){
	  int gauss_gid = mesh.gauss_in_elem(elem_gid, gauss_lid);
	  res_mass(index) +=   ref_elem.ref_nodal_basis(gauss_lid, vertex)
		               * mat_pt.density(gauss_gid) 
                               * ref_elem.ref_node_g_weights(gauss_lid)
		               * mesh.gauss_pt_det_j(gauss_gid)
		               * ref_elem.ref_nodal_basis(gauss_lid, index);
	}// end loop over gauss_lid
      }// end loop over index
      
      //--- end Mass Vector ---//


      //--- M_l.(v^{r}-v^{n}) ---/
      for (int index = 0; index < num_basis; index++){
	int dof_lid = ref_elem.vert_node_map(index);
	int dof_gid = mesh.nodes_in_elem(elem_gid, dof_lid);
	Mv += res_mass( index )*( node.vel( t_step, dof_gid, dim) - node.vel( 0, dof_gid, dim) );
      }// end loop over index
      //std::cout << " M.v at dim " << dim << " is " << Mv(dim) << std::endl;
      
      //-- end M_l.(v^{r}-v^{n}) ---//

      CArray <real_t> ones(ref_elem.num_dual_basis());
      for (int i = 0; i < ref_elem.num_dual_basis(); i++) ones(i) = 1.0;
      
      real_t force = 0.0;      

      for (int dof = 0; dof < ref_elem.num_dual_basis(); dof++){
	force += ones(dof)*0.5*( elem_state.force_tensor(t_step,elem_gid, vertex, dof, dim) + elem_state.force_tensor(0, elem_gid, vertex, dof, dim));
      }// end loop over dof
      //std::cout << " force at dim " << dim << " is " << force(dim) << std::endl;
      
      //--- Assign Values to Galerkin/Rusanov Residual ---//
      real_t sub_factor = 1.0/((real_t)num_correction_steps - (real_t)t_step);
      galerkin_res(vertex) = Mv + sub_factor*dt*force;
      
      //std::cout << "kinematic nodal res = "<< galerkin_res(vertex)<< " at node "<< node_gid << " elem "<< elem_gid << " and dim " << dim << std::endl;

    }// end loop over vertex
      
  /// end nodal res computation in element///


  // compute total residual //
    real_t total_res = 0.0;
      
    for (int dof = 0; dof < ref_elem.num_basis(); dof++){
        total_res += galerkin_res(dof);
    }// end loop over dof 

    //std::cout << total_res << std::endl;
    real_t inv_total_res = 1.0/total_res;
    for (int dof = 0; dof < ref_elem.num_basis(); dof++){
      
      real_t num = 0.0;
      num = std::max( 0.0, galerkin_res( dof )*inv_total_res );
      real_t denom = 0.0;
      for (int k = 0; k < ref_elem.num_basis(); k++){
        denom += std::max(0.0, galerkin_res(k)*inv_total_res );
      }// end loop over k
      
      int node_lid = ref_elem.vert_node_map(dof);
      int node_gid = mesh.nodes_in_elem(elem_gid,node_lid);
      if (node_gid == dof_gid){
        sum_res += galerkin_res(dof);
        //sum_res += (num/denom)*total_res;
      }
    }// end loop over dof

  }// end loop over elem_lid
}


/*
    // compute psi coeffs //
    auto psi_coeffs = CArray <real_t> (ref_elem.num_basis());
    for (int i = 0; i < num_basis; i++) psi_coeffs(i) = 0.0;

    for (int dof = 0; dof < ref_elem.num_basis(); dof++){
      real_t num = 0.0;
      num = std::max( 0.0, galerkin_res( dof )/total_res );
      real_t denom = 0.0;
      for (int k = 0; k < ref_elem.num_basis(); k++){
        denom += std::max(0.0, galerkin_res(k)/total_res );
      }// end loop over k

      psi_coeffs( dof ) = num/denom;
 
    }// end loop over dof
*/ 
/* 
  // check that betas sum to 1 //
      real_t sum = 0.0;
      for (int dof = 0; dof < ref_elem.num_basis(); dof++){
        sum += psi_coeffs(dof);
      }// end loop over dof
      std::cout << " sum of betas in elem "<< elem_gid << " in dim " << dim<< " is " << sum << std::endl;
*/ 

/*
    // compute limited residual //
    auto limited_res = CArray <real_t> ( ref_elem.num_basis() );
    for (int i = 0; i < num_basis; i++) limited_res(i) = 0.0;

    for (int dof = 0; dof < ref_elem.num_basis(); dof++){
      limited_res(dof) = psi_coeffs(dof)*total_res;
      //std::cout << limited_res(dof) << std::endl;
    }// end loop over dof
*/    
