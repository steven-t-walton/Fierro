#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"

#define PI 3.14159265

using namespace utils;

void interp_vel(){

   for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
      for (int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++){
        int gauss_gid = mesh.gauss_in_elem(elem_gid, gauss_lid);
        int node_gid = mesh.nodes_in_elem(elem_gid, gauss_lid);
        
        real_t interp_a[num_dim];
        for (int i =0; i < num_dim; i++) interp_a[i] =0.0;
        auto interp = ViewCArray <real_t> ( &interp_a[0], num_dim);
        
        for (int dim = 0; dim < mesh.num_dim(); dim++){
          for (int vert = 0; vert < ref_elem.num_basis(); vert++){
            interp(dim) += ref_elem.ref_nodal_basis( gauss_lid, vert ) * elem_state.vel_coeffs(num_correction_steps, elem_gid, vert, dim);
          }// end loop over vertex
        }// end loop over dim
      
          
        for (int dim = 0; dim < num_dim; dim++){
          node.vel(1, node_gid, dim) = interp(dim);
        }

/*      
        node.vel(1, node_gid, 0) = sin(PI * mesh.node_coords(node_gid, 0)) * cos(PI * mesh.node_coords(node_gid, 1));
        node.vel(1, node_gid, 1) = -1.0*cos(PI * mesh.node_coords(node_gid, 0)) * sin(PI * mesh.node_coords(node_gid, 1));
        node.vel(1, node_gid, 2) = 0.0;
*/

/*
        std::cout << node.vel(1, node_gid, 0) - sin(PI * mesh.node_coords(node_gid, 0)) * cos(PI * mesh.node_coords(node_gid, 1))<< std::endl;
        std::cout << node.vel(1, node_gid, 1) + 1.0*cos(PI * mesh.node_coords(node_gid, 0)) * sin(PI * mesh.node_coords(node_gid, 1))<< std::endl; 
        std::cout << node.vel(1, node_gid, 2) << std::endl; 
*/
      }// end loop over node_lid
    }// end loop over elements

}// end interp_vel()

/*
void interp_pos(int update){

  for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
    for (int node_lid = 0; node_lid < mesh.num_nodes_in_elem(); node_lid++){
      int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);

      for (int dim = 0; dim < mesh.num_dim(); dim++){
        node.coords(1, node_gid, dim) = 0.0;
      	for (int vertex = 0; vertex < ref_elem.num_basis(); vertex++){
          node.coords( 1, node_gid, dim ) += ref_elem.ref_nodal_basis( node_lid, vertex ) * elem_state.BV_pos_coeffs( update, elem_gid, vertex, dim );
	}// end loop over vertex
      }// end loop over dim
      
    }// end loop over node_lid
  }// end loop over elem_gid
  
}// end interp_pos()
*/

void interp_ie(){
  
  for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
    for (int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++){
      
      int gauss_gid = mesh.gauss_in_elem(elem_gid, gauss_lid);
      int node_gid = mesh.nodes_in_elem(elem_gid, gauss_lid);
      real_t interp = 0.0;

      for (int vert = 0; vert < ref_elem.num_dual_basis(); vert++){
        interp += ref_elem.ref_nodal_dual_basis( gauss_lid, vert ) * elem_state.sie_coeffs(num_correction_steps, elem_gid, vert);
      }// end loop over vertex
      
      real_t source = 0.0;

      source = 3.141592653589/(4.0*(0.66666667))* ( cos(3.0*3.141592653589 * node.coords(1,node_gid,0))  * cos( 3.141592653589 * node.coords(1,node_gid,1)) - cos( 3.141592653589 * node.coords(1,node_gid,0) ) * cos( 3.0*3.141592653589 * node.coords(1,node_gid, 1) ) );
      
      mat_pt.sie(1, gauss_gid) = interp;//+source;

      //std::cout<< mat_pt.sie(1,gauss_gid) << std::endl;

      mat_pt.ie(gauss_gid) = mat_pt.sie(1,gauss_gid);   
    
    }// end loop over gauss_lid
  }// end loop over elem_gid


}// end interp_ie

