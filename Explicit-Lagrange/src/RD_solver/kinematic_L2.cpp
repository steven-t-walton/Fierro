


#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"

using namespace utils;

void get_kinematic_L2(int t_step){


  int num_basis = ref_elem.num_basis();
  int num_dim = mesh.num_dim();

#pragma omp simd

  for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
    for (int vertex = 0; vertex < ref_elem.num_basis(); vertex++){
      int node_lid = elem.vert_node_map( vertex );
      int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);
      for (int dim = 0; dim < num_dim; dim++){
          elem_state.kinematic_L2(t_step, node_gid, dim) = 0.0;
      }// end loop over dim
    }// end loop over vertex
  }// end loop over elem_gid

  // Create Nodal Res Tensor//
  int nodal_res_size = mesh.num_elems()*mesh.num_nodes()*mesh.num_dim();
  real_t nodal_res_a[nodal_res_size];
  for (int i = 0; i < nodal_res_size; i++) nodal_res_a[i] = 0.0;
  auto nodal_res = ViewCArray <real_t> ( &nodal_res_a[0], mesh.num_nodes(), mesh.num_elems(), mesh.num_dim() );

  for(int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
  
      // Compute alpha_E //	
      real_t alpha = 0.0;
      
      for (int node_lid = 0; node_lid < mesh.num_nodes_in_elem(); node_lid++){
        int n_gid = mesh.nodes_in_elem(elem_gid, node_lid);
	int gauss_gid = mesh.gauss_in_elem(elem_gid, node_lid);

        real_t alpha_a = 0.0;//[3];
        //for (int i = 0; i < 3; i++) alpha_a[i] = 0.0;

        real_t speed = sqrt( node.vel(0, n_gid, 0)*node.vel(0, n_gid, 0) 
			+ node.vel(0, n_gid,  1)*node.vel(0, n_gid, 1)
			+ node.vel(0, n_gid,  2)*node.vel(0, n_gid, 2) );
	
	alpha_a = speed + mat_pt.sspd(gauss_gid);
	alpha = alpha > alpha_a ? alpha : alpha_a;

    }// end loop over node_lid	  

    real_t max_density = 0.0;

    for(int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++){
 
      int gauss_gid = mesh.gauss_in_elem(elem_gid, gauss_lid);

      real_t vol_gauss = ref_elem.ref_node_g_weights(gauss_lid)*mesh.gauss_pt_det_j(gauss_gid);
      mat_pt.density(gauss_gid) = mat_pt.mass(gauss_gid)/vol_gauss;
      //std::cout << mat_pt.density(gauss_gid) << std::endl;
      max_density = max_density > mat_pt.density(gauss_gid) ? max_density : mat_pt.density(gauss_gid);
    } // end loop gauss
    
    alpha = alpha*max_density;
    
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
    
    alpha = max_length*alpha; 
    
    // FOR TG multiply by small coefficient //
    alpha = 0.001*alpha;

    // Compute residual //
    for (int vertex = 0; vertex < num_basis; vertex++){
      
      int node_lid = elem.vert_node_map(vertex);
      int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);
      int g_gid = mesh.gauss_in_elem(elem_gid, node_lid);
	
      // View for mass_vector //
      real_t res_mass_a[num_basis];
      for (int i = 0; i < num_basis; i++) res_mass_a[i] = 0.0;
      auto res_mass = ViewCArray <real_t> ( &res_mass_a[0], num_basis );

      // Create view for \sum_{q} M_{pq} \delta u^r_q //      
      real_t Mv_a[num_dim];
      for (int i = 0; i < num_dim; i++) Mv_a[i] = 0.0;
      auto Mv = ViewCArray <real_t> ( &Mv_a[0], num_dim );

      // Create ViewCArray to store volume integral over cells in elem of force //
      real_t force_a[num_dim];
      for (int i = 0; i < num_dim; i++) force_a[i] = 0.0;
      auto force = ViewCArray <real_t> (&force_a[0], num_dim);

      //--- Mass Vector ---//
        
      for (int index = 0; index < num_basis; index++){
        for (int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++){
	  int gauss_gid = mesh.gauss_in_elem(elem_gid, gauss_lid);
	  res_mass(index) += ref_elem.ref_nodal_basis(gauss_lid, vertex)
		               * mat_pt.density(gauss_gid) 
                               * ref_elem.ref_node_g_weights(gauss_lid)
		               * mesh.gauss_pt_det_j(gauss_gid)
		               * ref_elem.ref_nodal_basis(gauss_lid, index);
	}// end loop over gauss_lid
      }// end loop over index
        
      //--- end Mass Vector ---//


      //--- Force ---//

      real_t volume_int_a[ num_dim ];
      auto volume_int = ViewCArray <real_t> (&volume_int_a[0], num_dim);
      for (int i = 0; i < num_dim; i++) volume_int(i) = 0.0;
      
      // Volume integral //
      for (int dim = 0; dim < num_dim; dim++){
	for (int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++){
	  int gauss_vert_gid = mesh.gauss_in_elem(elem_gid, gauss_lid);
	  real_t J_inv_dot_grad_phi = 0.0;
	  for (int k = 0; k < num_dim; k++){
	    J_inv_dot_grad_phi += mesh.gauss_pt_jacobian_inverse(gauss_vert_gid, dim, k) * ref_elem.ref_nodal_gradient(gauss_lid, vertex, k);
	  }// end loop over k
	  volume_int(dim) += mat_pt.pressure(gauss_vert_gid)
	        	       * J_inv_dot_grad_phi
			       * ref_elem.ref_node_g_weights(gauss_lid)
			       * mesh.gauss_pt_det_j(gauss_vert_gid);
	}// end loop over gauss_lid
      }// end loop over dim

      real_t surface_int_a[ num_dim ];
      auto surface_int = ViewCArray <real_t> (&surface_int_a[0], num_dim);
      for (int i = 0; i < num_dim; i++) surface_int(i) = 0.0;        
 
      build_corner_normals(); 	
      // Surface integral
      for(int cell_lid = 0; cell_lid < mesh.num_cells_in_elem(); cell_lid++){

        int cell_gid = mesh.cells_in_elem(elem_gid, cell_lid);
        for(int node_lid = 0; node_lid < mesh.num_nodes_in_cell(); node_lid++){
          int node_gid = mesh.nodes_in_cell(cell_gid, node_lid);
	  int gauss_gid = mesh.gauss_in_cell(cell_gid, node_lid);
          int corner_lid = node_lid;
          int corner_gid = mesh.corners_in_cell(cell_gid, corner_lid);  // node_lid = corner_lid

          // Get the id for the node in the reference element
          int node_rid = ref_elem.cell_nodes_in_elem(cell_lid, node_lid);

          for (int dim = 0; dim < num_dim; dim ++){
            surface_int(dim) += ref_elem.ref_nodal_basis(node_rid, vertex)
		                * corner.normal(corner_gid, dim) 
				* mat_pt.pressure(gauss_gid);
	  }// end loop over dim
        }// end loop over nodes/corners in a cell
      } // end loop over cells in an element

      for (int dim = 0; dim < num_dim; dim++){
        force(dim) =  surface_int(dim) - volume_int(dim);
	//std::cout << force(dim) << std::endl;
      }// end loop over dim
      //--- end Force ---//
      
      //--- Artificial Viscosity ---//
      real_t vel_bar_a[num_dim];
      for (int i = 0; i < num_dim; i++) vel_bar_a[i] = 0.0;
      auto vel_bar = ViewCArray <real_t> ( &vel_bar_a[0], num_dim );

      real_t Q_a[num_dim];
      for (int i = 0; i < num_dim; i++) Q_a[i] = 0.0;
      auto Q = ViewCArray <real_t> ( &Q_a[0], num_dim );
       
      // Compute vel_bar //
        for (int dim = 0; dim < num_dim; dim++){
          for (int vertex = 0; vertex < ref_elem.num_basis(); vertex++){
	      //int vert_lid = elem.vert_node_map(vertex);
	      //int vert_gid = mesh.nodes_in_elem(elem_gid, vert_lid);
	      //vel_bar(dim, current) += node.vel(0, vert_gid, dim)/ref_elem.num_basis();
            vel_bar(dim) += elem_state.vel_coeffs(t_step, elem_gid, vertex, dim)/ref_elem.num_basis(); 
	  }
	}
      // Fill Q //
          
        for (int dim = 0; dim < num_dim; dim++){ 
	  // Q = alpha_E*(vel - vel_bar) //
          //Q(dim, current) = alpha*(node.vel(0, node_gid, dim) - vel_bar(dim, current)); 
	  Q(dim) = alpha*(elem_state.vel_coeffs(t_step,elem_gid, vertex,dim) - vel_bar(dim));
          //std::cout <<  "Q at dim " << dim << " is " << Q(dim, current) << std::endl;  
        }// end loop over dim for Q 
      //--- end Artificial Viscosity ---//

      //--- M_l.(v^{r}-v^{n}) ---/
      for (int dim = 0; dim < num_dim; dim++){
        for (int index = 0; index < num_basis; index++){
	  Mv(dim) += res_mass( index )*(elem_state.vel_coeffs( t_step, elem_gid, index, dim) - elem_state.vel_coeffs( 0, elem_gid, index, dim));
        }// end loop over index
      }// end loop over dim
      //-- end M_l.(v^{r}-v^{n}) ---//

        /*
        for (int dim = 0; dim < num_dim; dim++){
          std::cout<< "Mv/dt at dim " << dim << " and correction_step " << t_step << " is " << Mv(dim) << std::endl;
	  std::cout << "force at dim " << dim << " is  " << dt*force(dim,current) << std::endl;
        }
        */
      //--- Assign Values to Galerkin/Rusanov Residual ---//
      for (int dim = 0; dim < num_dim; dim++){
        // Assign values to nodal res
        nodal_res(node_gid, elem_gid, dim ) = Mv(dim)/dt +  force( dim) + Q(dim);
      }// end loop over dim


    }// end loop over vertex
  }// end loop over elem_gid
      
  /// end nodal res computation ///


  // compute total residual //
  int total_res_size = num_dim*mesh.num_elems();
  real_t total_res_a[total_res_size];
  for (int i = 0; i < total_res_size; i++) total_res_a[i] = 0.0;
  auto total_res = ViewCArray <real_t> ( &total_res_a[0], mesh.num_elems(), num_dim);
      
  for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
    for (int dim = 0; dim < num_dim; dim++){
      for (int vertex = 0; vertex < ref_elem.num_basis(); vertex++){
        int node_lid = elem.vert_node_map( vertex);
        int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);
        total_res(elem_gid, dim) += nodal_res(node_gid, elem_gid, dim);
      }// end loop over node_lid
    }// end loop over dim
  }// end loop over elem_gid

/*
  for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
    for (int dim = 0; dim < num_dim; dim++){
      std::cout << total_res(elem_gid, dim) << std::endl;
    }
  }
*/

  // compute psi coeffs //
  int psi_coeffs_size = num_dim*ref_elem.num_basis()*mesh.num_elems();
  real_t psi_coeffs_a[psi_coeffs_size];
  for (int i = 0; i < psi_coeffs_size; i++) psi_coeffs_a[i] = 0.0;
  auto psi_coeffs = ViewCArray <real_t> ( &psi_coeffs_a[0], mesh.num_elems(), ref_elem.num_basis(), num_dim);

  for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
    for (int vertex = 0; vertex < ref_elem.num_basis(); vertex++){
      int node_lid = elem.vert_node_map(vertex);
      int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);
      for (int dim = 0; dim < mesh.num_dim(); dim++){
        real_t num = 0.0;
        num = std::max( 0.0, nodal_res(node_gid, elem_gid, dim)/total_res(elem_gid, dim) );
        real_t denom = 0.0;
        for (int k = 0; k < ref_elem.num_basis(); k++){
          int node_k_lid = elem.vert_node_map( k );
          int node_k_gid = mesh.nodes_in_elem(elem_gid, node_k_lid);
          denom += std::max(0.0, nodal_res( node_k_gid, elem_gid, dim)/total_res(elem_gid,dim) );
        }// end loop over k

        psi_coeffs( elem_gid, vertex, dim ) = num/denom;

      }// end loop over dim
    }// end loop over vertex
  }// end loop over elem_gid
 
/* 
  // check that betas sum to 1 //
  for (int dim = 0; dim < mesh.num_dim(); dim++){
    for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
      real_t sum = 0.0;
      for (int vertex = 0; vertex < ref_elem.num_basis(); vertex++){
        sum += psi_coeffs( elem_gid, vertex, dim);
      }// end loop over vertex
      std::cout << " sum of betas in elem "<< elem_gid << " in dim " << dim<< " is " << sum << std::endl;
    }// end loop over elem_gid
  }// end loop over dim
*/ 


  // compute limited residual //
  int limited_res_size = num_dim*ref_elem.num_basis()*mesh.num_elems();
  real_t limited_res_a[limited_res_size];
  for (int i = 0; i < limited_res_size; i++) limited_res_a[i] = 0.0;
  auto limited_res = ViewCArray <real_t> ( &limited_res_a[0], mesh.num_elems(), ref_elem.num_basis(), num_dim);

  for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
    for (int vertex = 0; vertex < ref_elem.num_basis(); vertex++){
      for (int dim = 0; dim < num_dim; dim++){
        limited_res(elem_gid, vertex, dim) = psi_coeffs( elem_gid, vertex, dim)*total_res(elem_gid, dim);
        //std::cout << limited_res(elem_gid, vertex, dim) << std::endl;
      }// end loop over dim
    }// end loop over vertex
  }// end loop over elem_gid


  for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
    for (int vertex = 0; vertex < ref_elem.num_basis(); vertex++){
      int node_lid = elem.vert_node_map( vertex );
      int g_gid = mesh.gauss_in_elem(elem_gid, node_lid);
      int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);
      int num_elems_in_vert = mesh.num_elems_in_node(node_gid);
      
      for (int dim = 0; dim < num_dim; dim++){
        for (int elems_in_vert = 0; elems_in_vert < num_elems_in_vert; elems_in_vert++){
          int elems_in_vert_gid = mesh.elems_in_node(node_gid, elems_in_vert);
          elem_state.kinematic_L2(t_step, node_gid, dim) += nodal_res(node_gid, elems_in_vert_gid, dim);
          //elem_state.kinematic_L2(t_step, node_gid, dim) += limited_res(elems_in_vert_gid, vertex, dim);       
        }// end loop over elems_in_vert
      }// end loop over dim
    
    }// end loop over vertex
  }// end loop over elem_gid

}
