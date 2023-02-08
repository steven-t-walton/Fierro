#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"

#define PI 3.14159265

using namespace utils;

void interp_vel(){

   for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
        for(int node_lid = 0; node_lid < mesh.num_gauss_in_elem(); node_lid++){

            int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);

	    CArray <real_t> interp_vel(mesh.num_dim());
            for(int i=0; i<mesh.num_dim(); i++) interp_vel(i) = 0.0;
         
            for (int dim = 0; dim < mesh.num_dim(); dim++){
    	        for(int basis_id = 0; basis_id < elem.num_basis(); basis_id++){
  
                    int node_basis_id = ref_elem.vert_node_map(basis_id);
                    int interp_gid = mesh.nodes_in_elem(elem_gid, node_basis_id);
                    interp_vel(dim) += node.vel(num_correction_steps, interp_gid, dim) * ref_elem.ref_nodal_basis(node_lid, basis_id);
                }// end loop over basis_id
            }// end loop over dim
            
            for (int dim = 0; dim < mesh.num_dim(); dim++){   
                node.vel(num_correction_steps, node_gid, dim) =  interp_vel(dim);
            }// end loop over dim

   	}//end loop over gauss_lid
   }// end loop over elements

}// end interp_vel()

void interp_ie(){
  
  for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
        for(int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++){

            int node_gid = mesh.nodes_in_elem(elem_gid, gauss_lid);

            int gauss_gid = mesh.gauss_in_elem(elem_gid, gauss_lid);

            real_t interp_sie = 0.0; 
         
    	    for(int basis_id = 0; basis_id < elem.num_dual_basis(); basis_id++){
  
                    int node_basis_id = ref_elem.dual_vert_node_map(basis_id);
                    int interp_gid = mesh.gauss_in_elem(elem_gid, node_basis_id);
                    interp_sie += mat_pt.sie(num_correction_steps, interp_gid) * ref_elem.ref_nodal_dual_basis(gauss_lid, basis_id);
            }// end looop over basis_id
            
            mat_pt.sie(num_correction_steps,gauss_gid) =  interp_sie;
            
	}// end loop over gauss_lid
  }// end loop over elems

}// end sie


      /*
      for (int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++){
        int gauss_gid = mesh.gauss_in_elem(elem_gid, gauss_lid);
        int node_gid = mesh.nodes_in_elem(elem_gid, gauss_lid);
        
        real_t interp_a[num_dim];
        for (int i =0; i < num_dim; i++) interp_a[i] =0.0;
        auto interp = ViewCArray <real_t> ( &interp_a[0], num_dim);
        
        for (int dim = 0; dim < mesh.num_dim(); dim++){
          for (int vert = 0; vert < ref_elem.num_basis(); vert++){
            interp(dim) += ref_elem.ref_nodal_basis( gauss_lid, vert ) * elem_state.vel_coeffs(t_step, elem_gid, vert, dim);
          }// end loop over vertex
        }// end loop over dim
      
          
        for (int dim = 0; dim < num_dim; dim++){
          node.vel(1, node_gid, dim) = interp(dim);
        }
	*/
/*
        for (int basis = 0; basis < ref_elem.num_basis(); basis++){
          for (int dim = 0; dim < mesh.num_dim(); dim++){
            elem_state.vel_coeffs( 0, elem_gid, basis, dim ) = elem_state.vel_coeffs( t_step, elem_gid, basis, dim );
          }// end loop over dim
        }// end loop over basis
*/      
/*      
        node.vel(1, node_gid, 0) = sin(PI * mesh.node_coords(node_gid, 0)) * cos(PI * mesh.node_coords(node_gid, 1));
        node.vel(1, node_gid, 1) = -1.0*cos(PI * mesh.node_coords(node_gid, 0)) * sin(PI * mesh.node_coords(node_gid, 1));
        node.vel(1, node_gid, 2) = 0.0;
*/

/*
        std::cout << node.vel(1, node_gid, 0) - sin(PI * mesh.node_coords(node_gid, 0)) * cos(PI * mesh.node_coords(node_gid, 1))<< std::endl;
        std::cout << node.vel(1, node_gid, 1) + 1.0*cos(PI * mesh.node_coords(node_gid, 0)) * sin(PI * mesh.node_coords(node_gid, 1))<< std::endl; 
        std::cout << node.vel(1, node_gid, 2) << std::endl; 

      }// end loop over node_lid
    */



/*
    for (int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++){
      
      int gauss_gid = mesh.gauss_in_elem(elem_gid, gauss_lid);
      int node_gid = mesh.nodes_in_elem(elem_gid, gauss_lid);
      real_t interp = 0.0;

      for (int vert = 0; vert < ref_elem.num_dual_basis(); vert++){
        interp += ref_elem.ref_nodal_dual_basis( gauss_lid, vert ) * elem_state.sie_coeffs(t_step, elem_gid, vert);
      }// end loop over vertex
      
      mat_pt.sie(1, gauss_gid) = interp;

      //std::cout<< mat_pt.sie(1,gauss_gid) << std::endl;

      mat_pt.ie(gauss_gid) = mat_pt.sie(1,gauss_gid);   
     // /
      for (int t_basis = 0; t_basis < ref_elem.num_dual_basis(); t_basis++){
        elem_state.sie_coeffs( 0, elem_gid, t_basis ) = elem_state.sie_coeffs( t_step, elem_gid, t_basis);
      }
      
    }// end loop over gauss_lid
*/	

