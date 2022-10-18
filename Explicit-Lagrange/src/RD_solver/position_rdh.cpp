#include<iostream>

#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"

using namespace utils;

void get_position_rdh(){

  int update = num_correction_steps - 1;

// Update position with v(t^{n+1}) //
#pragma omp simd  
          for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){

            for(int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++){

              int node_gid = mesh.nodes_in_elem(elem_gid, gauss_lid);

              auto vel_update = ViewCArray <real_t> ( &node.vel( update, node_gid, 0), mesh.num_dim() );

              // get the global id of the gauss point
              int gauss_gid = mesh.gauss_in_elem(elem_gid, gauss_lid);

#pragma omp simd
              for (int dim = 0; dim < mesh.num_dim(); dim++){
                node.coords(1, node_gid, dim) = node.coords(0, node_gid, dim) + dt * vel_update( dim );
              }// end loop over dim

              real_t interp_pos[mesh.num_dim()];
              for(int i=0; i<mesh.num_dim(); i++) interp_pos[i] = 0.0;

              for (int dim = 0; dim < mesh.num_dim(); dim++){
                for(int basis_id = 0; basis_id < elem.num_basis(); basis_id++){

                  int node_basis_id = elem.vert_node_map(basis_id);
                  int interp_gid = mesh.gauss_in_elem(elem_gid, node_basis_id);
                  interp_pos[dim] += node.coords(1,interp_gid,dim) * ref_elem.ref_nodal_basis(gauss_lid, basis_id);

                }// end loop over basis id
              }// end loop over dim

              // Save interpolated position back to gauss point
              for (int dim = 0; dim < mesh.num_dim(); dim++){
                node.coords( 1, node_gid, dim) = interp_pos[dim];
                interp_pos[dim] = 0.0;
                mesh.node_coords(node_gid, dim) = node.coords( 1 , node_gid, dim);
//                 std::cout << " " << std::endl;  
//              std::cout<< "pos at time "<< TIME+dt <<" cycle "<< cycle <<" node " << node_gid << " and dim " << dim <<" is "<< node.coords(update,node_gid,dim) << std::endl;

	      }// end loop over dim
            }// end loop over gauss_lid
	  }// end loop over elem_gid

        //std::cout << "Calculating Jacobian at gauss points" << std::endl;
        get_gauss_pt_jacobian(mesh, ref_elem);
    
       //std::cout << "Calculating Jacobian at gauss points in cell" << std::endl;
        get_gauss_cell_pt_jacobian(mesh, ref_elem);

       //std::cout << "Before volume from Jacobian"  << std::endl;
        get_vol_jacobi(mesh, ref_elem);
}
