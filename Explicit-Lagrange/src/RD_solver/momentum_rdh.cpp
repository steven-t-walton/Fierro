/* nodal residual */

#include<iostream>
#include<math.h>
#include<algorithm>
#include<vector>

#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"

#define PI 3.14159265

using namespace utils;

void update_velocity(int t_step){
   
  int num_basis = ref_elem.num_basis();
  int num_dim = mesh.num_dim();

  int current = t_step;
  int update = t_step+1;

#pragma omp simd

  /// update velocity coefficients //
  for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
    for (int vertex = 0; vertex < ref_elem.num_basis(); vertex++){
      int node_lid = elem.vert_node_map( vertex );
      int g_gid = mesh.gauss_in_elem(elem_gid, node_lid);
      int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);
      int num_elems_in_vert = mesh.num_elems_in_node(node_gid);
      
      // Compute lumped mass
      real_t lumped_mass = 0.0;
      
      for (int elems_in_vert = 0; elems_in_vert < num_elems_in_vert; elems_in_vert++){
        int elems_in_vert_gid = mesh.elems_in_node(node_gid, elems_in_vert);
        for (int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++){
	  int gauss_gid = mesh.gauss_in_elem(elems_in_vert_gid, gauss_lid);
	  lumped_mass += ref_elem.ref_node_g_weights(gauss_lid)
		         * mat_pt.density(gauss_gid)
			 * mesh.gauss_pt_det_j(gauss_gid)
	                 * ref_elem.ref_nodal_basis(gauss_lid, vertex);
        }// end loop over gauss_lid
      }// end loop over elems_in_node_lid
      if (lumped_mass <= 0.0){
        std::cout << " kinematic lumped mass is negative " << lumped_mass << std::endl;
      }
      // L^1(u^{k+1}) = L^1(u^k) - L^2(u^k) //
      for (int dim = 0; dim < num_dim; dim++){
        mat_pt.velocity(update, g_gid, dim) = mat_pt.velocity(current,g_gid, dim) 
		                                               - (dt/lumped_mass)*elem_state.kinematic_L2(current, node_gid, dim);
      }
    }// end loop over vertex
  }// end loop over elem_gid

}// end get_nodal_res

