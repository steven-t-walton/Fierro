#include<iostream>
#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"

using namespace utils;

void get_alpha_E(){

  int num_basis = ref_elem.num_basis();
  int num_dim = mesh.num_dim();
  
  for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
    elem_state.alpha_E(elem_gid) = 0.0;
  }  

  // Compute alpha_E //	
  for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){    
    for (int node_lid = 0; node_lid < mesh.num_nodes_in_elem(); node_lid++){
      int n_gid = mesh.nodes_in_elem(elem_gid, node_lid);
      int gauss_gid = mesh.gauss_in_elem(elem_gid, node_lid);

      real_t alpha_a = 0.0;

      real_t speed = sqrt( node.vel(0, n_gid, 0)*node.vel(0, n_gid, 0) 
			+ node.vel(0, n_gid,  1)*node.vel(0, n_gid, 1)
			+ node.vel(0, n_gid,  2)*node.vel(0, n_gid, 2) );
      
      alpha_a = speed + mat_pt.sspd(gauss_gid);
      elem_state.alpha_E(elem_gid) = elem_state.alpha_E(elem_gid) > alpha_a ? elem_state.alpha_E(elem_gid) : alpha_a;

    }// end loop over node_lid	  

    real_t max_density = 0.0;

    for(int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++){
 
      int gauss_gid = mesh.gauss_in_elem(elem_gid, gauss_lid);

      //real_t vol_gauss = ref_elem.ref_node_g_weights(gauss_lid)*mesh.gauss_pt_det_j(gauss_gid);
      //mat_pt.density(gauss_gid) = mat_pt.mass(gauss_gid)/vol_gauss;
      //std::cout << mat_pt.density(gauss_gid) << std::endl;
      max_density = max_density > mat_pt.density(gauss_gid) ? max_density : mat_pt.density(gauss_gid);
    } // end loop gauss
    
    elem_state.alpha_E(elem_gid) = elem_state.alpha_E(elem_gid)*max_density;
    
    real_t max_length = 0.0;
    for (int vertex = 0; vertex < num_basis; vertex++){
      real_t length[num_dim];
      for (int i = 0; i < num_dim; i++) length[i] = 0.0;
      real_t mag_length = 0.0;
      for (int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++){
        int gauss_gid = mesh.gauss_in_elem(elem_gid, gauss_lid); 
        for (int dim = 0; dim < num_dim; dim++){
          length[dim] += ref_elem.ref_nodal_gradient(gauss_lid, vertex, dim)
	  	       * mesh.gauss_pt_det_j(gauss_gid)
		       * ref_elem.ref_node_g_weights(gauss_lid);
        }// end loop over dim
      }// end loop over gauss_lid

      mag_length = sqrt( length[0]*length[0] + length[1]*length[1] +  length[2]*length[2]  );
      
      max_length = max_length > mag_length ? max_length : mag_length;


    }// end loop over vertex
    
    elem_state.alpha_E(elem_gid) = max_length*elem_state.alpha_E(elem_gid); 
    
    // FOR TG multiply by small coefficient //
    elem_state.alpha_E(elem_gid) = 0.0*elem_state.alpha_E(elem_gid);
    //std::cout << elem_state.alpha_E(elem_gid) << std::endl;

  }// end loop over elem_gid
}// end get_alpha_E
