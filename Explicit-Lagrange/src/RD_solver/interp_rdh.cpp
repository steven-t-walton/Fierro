#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"

using namespace utils;

void interp_vel(int update){

  for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
    for (int node_lid = 0; node_lid < mesh.num_nodes_in_elem(); node_lid++){
      int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);

      for (int dim = 0; dim < mesh.num_dim(); dim++){
        node.vel(1, node_gid, dim) = 0.0;
      	for (int vertex = 0; vertex < ref_elem.num_basis(); vertex++){
          node.vel( 1, node_gid, dim ) += ref_elem.ref_nodal_basis( node_lid, vertex ) * elem_state.BV_vel_coeffs( update, elem_gid, vertex, dim );
  	}// end loop over vertex
      }// end loop over dim
      
    }// end loop over node_lid
  }// end loop over elem_gid

}// end interp_vel()

/*
void interp_pos(int update){


  for (int dim = 0; dim < mesh.num_dim(); dim++){
    for (int node_gid = 0; node_gid < mesh.num_nodes(); node_gid++){
      node.coords(1, node_gid,dim) = 0.0;
    }
  }

  for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
    for (int node_lid = 0; node_lid < mesh.num_nodes_in_elem(); node_lid++){
      int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);
      // interpolate //
      for (int dim = 0; dim < mesh.num_dim(); dim++){
        for (int vertex = 0; vertex < ref_elem.num_basis(); vertex++){
          node.coords( 1, node_gid, dim) += ref_elem.ref_nodal_basis( node_lid, vertex ) * elem_state.BV_pos_coeffs( update, elem_gid, vertex, dim ) ;
        }// end loop over vertex
      }// end loop over dim

    }// end loop over node_lid
  }// end loop over elem_gid
}// end interp_pos()
*/

void interp_ie(int update){

}// end interp_ie

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
      real_t coeff_a[ref_elem.num_basis() * mesh.num_dim()];
      auto coeff = ViewCArray <real_t> ( &coeff_a[0], ref_elem.num_basis(), mesh.num_dim() );
    
      // get coeffs for basis expansion //
      for (int dim = 0; dim < mesh.num_dim(); dim++){
        for (int vertex = 0; vertex < ref_elem.num_basis(); vertex++){
	  for (int k = 0; k < ref_elem.num_basis(); k++){
    	    for (int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++){
	      int gauss_gid = mesh.gauss_in_elem(elem_gid, gauss_lid);
	      coeff(vertex, dim)  += ref_elem.ref_nodal_basis(gauss_lid,k) * node.vel(update, node_gid, dim)
		                     * ref_elem.ref_node_g_weights(gauss_lid) * mesh.gauss_pt_det_j(gauss_gid)
		                     /( ref_elem.ref_nodal_basis(gauss_lid,k)*ref_elem.ref_nodal_basis(gauss_lid, vertex)
				     * ref_elem.ref_node_g_weights(gauss_lid) *mesh.gauss_pt_det_j(gauss_gid) );
	    }// end loop over gauss_lid
	  }// end loop over k
        }// end loop over vertex
      }// end loop over dim
*/
   /*
      real_t coeff_a[ref_elem.num_basis() * mesh.num_dim()];
      auto coeff = ViewCArray <real_t> ( &coeff_a[0], ref_elem.num_basis(), mesh.num_dim() );
   
      // get coeffs for basis expansion //
      for (int dim = 0; dim < mesh.num_dim(); dim++){
        for (int vertex = 0; vertex < ref_elem.num_basis(); vertex++){
	  for (int k = 0; k < ref_elem.num_basis(); k++){
            for (int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++){
	    
	      int gauss_gid = mesh.gauss_in_elem(elem_gid, gauss_lid);
	      coeff(vertex, dim)  += ref_elem.ref_nodal_basis(gauss_lid,k) * node.coords(update, node_gid, dim)
		                     * ref_elem.ref_node_g_weights(gauss_lid) * mesh.gauss_pt_det_j(gauss_gid)
		                     /( ref_elem.ref_nodal_basis(gauss_lid,k)*ref_elem.ref_nodal_basis(gauss_lid, vertex)
				     * ref_elem.ref_node_g_weights(gauss_lid) *mesh.gauss_pt_det_j(gauss_gid) );
	    }// end loop over gauss_lid
	  }// end loop over k
        }// end loop over vertex
      }// end loop over dim
    */


  // interpolation only over vertices //
  /*
        real_t interp_vel_a[ ref_elem.num_basis()*mesh.num_dim() ];
        auto interp_vel = ViewCArray <real_t> ( &interp_vel_a[0],  ref_elem.num_basis(), mesh.num_dim() );
        for(int j = 0; j < mesh.num_dim(); j++){
          for (int i = 0; i < ref_elem.num_basis(); i++){
            interp_vel(i,j) = 0.0;
          }
        }

        for (int dim = 0; dim < mesh.num_dim(); dim++){
          for(int vert_id = 0; vert_id < ref_elem.num_basis(); vert_id++){

            int node_id = ref_elem.vert_node_map( vert_id );
            int interp_gid = mesh.nodes_in_elem( elem_gid, node_id );
            for (int k = 0; k < ref_elem.num_basis(); k++){          
              interp_vel(vert_id, dim) += node.vel( update, interp_gid, dim ) * ref_elem.ref_nodal_basis( node_id, k );
            }// end loop over k
            node.vel( update, interp_gid, dim ) = interp_vel(vert_id,dim);
            //std::cout << " interpolated vel at dim " << dim << " is " << node.vel( update, interp_gid, dim ) << std::endl;     
          }// end loop over vert_id
        }// end loop over dim
  */
  
/*
  // interpolation over all nodes //
    
    real_t interp_vel_a[ mesh.num_nodes()*mesh.num_dim() ];
    auto interp_vel = ViewCArray <real_t> ( &interp_vel_a[0],  mesh.num_nodes(), mesh.num_dim() );
    for(int j = 0; j < mesh.num_dim(); j++){
      for (int i = 0; i < mesh.num_nodes(); i++){
        interp_vel(i,j) = 0.0;
      }
    }

    for (int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++){
      int node_gid = mesh.nodes_in_elem( elem_gid, gauss_lid );
      for (int dim = 0; dim < mesh.num_dim(); dim++){
        for (int k = 0; k < ref_elem.num_basis(); k++){
        
          interp_vel(node_gid, dim) = node.vel( update, node_gid, dim ) * ref_elem.ref_nodal_basis( gauss_lid, k );
	
	}// end loop over k        
	node.vel( update, node_gid, dim) = interp_vel( node_gid, dim); 
      }// end loop over dim

    }// end loop over node_lid
    
  }// end loop over elem_gid
*/


  /*
   for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
  // interpolation only at vertices //
/
        real_t interp_pos_a[ ref_elem.num_basis()*mesh.num_dim() ];
        auto interp_pos = ViewCArray <real_t> ( &interp_pos_a[0],  ref_elem.num_basis(), mesh.num_dim() );
        for(int j = 0; j < mesh.num_dim(); j++){
          for (int i = 0; i < ref_elem.num_basis(); i++){
            interp_pos(i,j) = 0.0;
          }
        }

        for (int dim = 0; dim < mesh.num_dim(); dim++){
          for(int vert_id = 0; vert_id < ref_elem.num_basis(); vert_id++){

            int node_id = ref_elem.vert_node_map( vert_id );
            int interp_gid = mesh.nodes_in_elem( elem_gid, node_id );
            for (int k = 0; k < ref_elem.num_basis(); k++){          
              interp_pos(vert_id, dim) += node.coords( update, interp_gid, dim ) * ref_elem.ref_nodal_basis( node_id, k );
            }// end loop over k
            node.coords( update, interp_gid, dim ) = interp_pos(vert_id,dim);
          }// end loop over vert_id
        }// end loop over dim


  // interpolation over all nodes //

    real_t interp_pos_a[ mesh.num_nodes()*mesh.num_dim() ];
    auto interp_pos= ViewCArray <real_t> ( &interp_pos_a[0],  mesh.num_nodes(), mesh.num_dim() );
    for(int j = 0; j < mesh.num_dim(); j++){
      for (int i = 0; i < mesh.num_nodes(); i++){
        interp_pos(i,j) = 0.0;
      }
    }

    for (int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++){
      int node_gid = mesh.nodes_in_elem( elem_gid, gauss_lid );
      for (int dim = 0; dim < mesh.num_dim(); dim++){
        for (int k = 0; k < ref_elem.num_basis(); k++){
          
          interp_pos(node_gid, dim) = node.coords( update, node_gid, dim ) * ref_elem.ref_nodal_basis( gauss_lid, k );
	
	}// end loop over k        
	node.coords( update, node_gid, dim) = interp_pos( node_gid, dim); 
      }// end loop over dim

    }// end loop over node_lid
    
  }// end loop over elem_gid
*/
