#include<iostream>

#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"

using namespace utils;

void get_position_rdh(int correction_step){

  int update = correction_step;//num_correction_steps;// - 1;
/*
  int rid_a[ref_elem.num_basis()];
  auto rid = ViewCArray <int> (&rid_a[0], ref_elem.num_basis());
  for (int m=0; m < ref_elem.num_basis(); m++) rid(m) = 0;

#pragma omp simd
  int ind = 0;
  for (int k = 0; k <  cbrt(ref_elem.num_basis()); k++){
    for (int j = 0; j < cbrt(ref_elem.num_basis()); j++){
      for (int i = 0; i < cbrt(ref_elem.num_basis()); i++){
        rid(ind) = ref_elem.node_rid(i,j,k);
        //std::cout << "node_rid at i = " << i << " j = " << j << " k = " << k << " is " << rid(ind) << std::endl;
        ind++;
      };// end loop over i
    };// end loop over j
  };// end loop over k
*/

// Update position with v(t^{n+1}) //
#pragma omp simd  
          for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){

            for( int node_lid = 0; node_lid < mesh.num_nodes_in_elem(); node_lid++){ //for(int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++){

              int node_gid = mesh.nodes_in_elem(elem_gid, node_lid );// gauss_lid);

              auto vel_update = ViewCArray <real_t> ( &node.vel( update, node_gid, 0), mesh.num_dim() );
              auto vel_n = ViewCArray <real_t> ( &node.vel(0, node_gid, 0 ), mesh.num_dim() );

              // get the global id of the gauss point
              //int gauss_gid = mesh.gauss_in_elem(elem_gid, gauss_lid);

#pragma omp simd
              for (int dim = 0; dim < mesh.num_dim(); dim++){
                node.coords(update, node_gid, dim) = node.coords(0, node_gid, dim) + 0.5*dt * ( vel_update( dim ) + vel_n(dim) );
              }// end loop over dim

  	 /*     
              real_t interp_pos[mesh.num_dim()];
              for(int i=0; i<mesh.num_dim(); i++) interp_pos[i] = 0.0;

              for (int dim = 0; dim < mesh.num_dim(); dim++){
	        for( int basis_n = 0; basis_n < ref_elem.num_basis(); basis_n++){
                  for(int basis_id = 0; basis_id < ref_elem.num_basis(); basis_id++){

                    int node_basis_id = ref_elem.vert_node_map(basis_id);
                    int interp_gid = mesh.nodes_in_elem(elem_gid, node_basis_id);
                    interp_pos[dim] += node.coords( update, interp_gid, dim) * ref_elem.ref_nodal_basis( rid(basis_id), basis_n );

                  }// end loop over basis id
		}// end loop over basis_n
              }// end loop over dim

              // Save interpolated position back to gauss point
              for (int dim = 0; dim < mesh.num_dim(); dim++){
                node.coords( update, node_gid, dim) = interp_pos[dim];
                interp_pos[dim] = 0.0;
                mesh.node_coords(node_gid, dim) = node.coords( update , node_gid, dim);
//                 std::cout << " " << std::endl;  
//              std::cout<< "pos at time "<< TIME+dt <<" cycle "<< cycle <<" node " << node_gid << " and dim " << dim <<" is "<< node.coords(update,node_gid,dim) << std::endl;

	      }// end loop over dim
	   */   

	    }// end loop over gauss_lid
	  }// end loop over elem_gid

}
