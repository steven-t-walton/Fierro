#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"

using namespace utils;

void get_kinematic_L2(int t_step, int dof_gid, int dim, real_t& sum_res){

  int num_basis = ref_elem.num_basis();
  int num_dim = mesh.num_dim();
  
  int num_elems_in_dof = mesh.num_elems_in_node(dof_gid);

  real_t inv_num_basis = 1.0/ref_elem.num_basis(); 
  
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
      
      // Compute artificial viscosity //
      real_t Q = 0.0;
      real_t u_bar = 0.0;
      real_t u_bar0 = 0.0;

      // Compute u_bar //
      for (int dof = 0; dof < ref_elem.num_basis(); dof++){
	  int n_lid = ref_elem.vert_node_map(dof);
	  int n_index = mesh.nodes_in_elem(elem_gid, n_lid);
          u_bar += node.vel(t_step, n_index, dim)*inv_num_basis; 
          u_bar0 += node.vel(0, n_index, dim)*inv_num_basis; 
      }
  
      // Compute Q //
      Q = 0.5*elem_state.alpha_E(elem_gid)*(node.vel(t_step, dof_gid, dim) - u_bar)
	  + 0.5*elem_state.alpha_E(elem_gid)*(node.vel(0, dof_gid, dim) - u_bar0);
      //if (dim == 2){ 
      //  std::cout << " artificial viscosity for dim "<< dim << " is " << Q << std::endl;
      //}

      //--- Assign Values to Galerkin/Rusanov Residual ---//
      galerkin_res(vertex) = Mv + dt*(force + Q);
      
      //std::cout << "kinematic nodal res = "<< galerkin_res(vertex)<< " at node "<< node_gid << " elem "<< elem_gid << " and dim " << dim << std::endl;

    }// end loop over vertex
      
  /// end nodal res computation in element///


  // compute total residual //
    real_t total_res = 0.0;
    real_t abs_total_res = 0.0;
    for (int dof = 0; dof < ref_elem.num_basis(); dof++){
        total_res += galerkin_res(dof);
	abs_total_res += abs(galerkin_res(dof));
    }// end loop over dof 
    //abs_total_res += 1.0e-10;

    //std::cout << "kinematic total res in elem "<< elem_gid << " and dof " << dof_gid << " is " << total_res << std::endl;
    //std::cout << "abs kinematic total res in elem "<< elem_gid << " and dof " << dof_gid << " is " << abs(total_res) << std::endl;
    
    real_t inv_total_res = 0.0;
    real_t inv_abs_total_res = 0.0;

    inv_total_res = 1.0/(total_res);
    inv_abs_total_res = 1.0/abs_total_res;

    // compute psi coeffs //
    auto psi_coeffs = CArray <real_t> (ref_elem.num_basis());
    for (int i = 0; i < num_basis; i++) psi_coeffs(i) = 0.0;
    
    for (int dof = 0; dof < ref_elem.num_basis(); dof++){
      if ( total_res > 0.0){
        
        real_t num = 0.0;
        num = std::max( 0.0, galerkin_res( dof )*inv_total_res );
      
        real_t denom = 0.0;
        for (int k = 0; k < ref_elem.num_basis(); k++){
          denom += std::max(0.0, galerkin_res(k)*inv_total_res );
        }// end loop over k
        psi_coeffs(dof) = num/denom;
      }
      else{
        psi_coeffs(dof) = inv_num_basis;
      };      
    }  

/* 
  // check that betas sum to 1 //
      real_t sum = 0.0;
      for (int dof = 0; dof < ref_elem.num_basis(); dof++){
        sum += psi_coeffs(dof);
      }// end loop over dof
      std::cout << " sum of betas in elem "<< elem_gid << " in dim " << dim<< " is " << sum << std::endl;
*/

    for (int dof = 0; dof < ref_elem.num_basis(); dof++){
      int node_lid = ref_elem.vert_node_map(dof);
      int node_gid = mesh.nodes_in_elem(elem_gid,node_lid);
      
      if (node_gid == dof_gid){
        if ( 0.0 < abs(total_res) ){
	  //std::cout<< " executing blended " << std::endl;
	  
	  //sum_res += galerkin_res(dof);
	  
          real_t limited_res = 0.0;
          real_t factor = 0.0;
	  
	  limited_res = psi_coeffs(dof)*total_res;
	  
	  //std::cout << " limited res is "<< limited_res << std::endl;
	  
	  factor = abs(total_res)*inv_abs_total_res;
	  
	  //std::cout << " factor is " << factor << std::endl;
	  //std::cout << " blended res is " <<  factor*galerkin_res(dof) + (1.0-factor)*limited_res << std::endl;
           
	  sum_res += (factor)*galerkin_res(dof) + (1.0-factor)*limited_res;
	}
	else{
	  
          //std::cout << " else condition met " << std::endl;
	  sum_res += 0.0;
	  
	  //sum_res += galerkin_res(dof);
	
	};
      }
    }// end loop over dof

  }// end loop over elem_lid

}


/*
    // compute limited residual //
    auto limited_res = CArray <real_t> ( ref_elem.num_basis() );
    for (int i = 0; i < num_basis; i++) limited_res(i) = 0.0;

    for (int dof = 0; dof < ref_elem.num_basis(); dof++){
      limited_res(dof) = psi_coeffs(dof)*total_res;
      //std::cout << limited_res(dof) << std::endl;
    }// end loop over dof
*/    
