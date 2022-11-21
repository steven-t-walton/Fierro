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

void get_nodal_res(int t_step){
   
  int num_basis = ref_elem.num_basis();
  int num_dim = mesh.num_dim();

  int current = t_step;
  int update = t_step+1;

#pragma omp simd
  
  std::vector<int> verts_temp; 
  for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
    for (int verts = 0; verts < num_basis; verts++){
      int node_lid = elem.vert_node_map(verts);
      int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);
      verts_temp.push_back(node_gid);
    }
  }

  std::sort(verts_temp.begin(), verts_temp.end());
  auto temp = std::unique(verts_temp.begin(), verts_temp.end());
  verts_temp.erase(temp, verts_temp.end());

  int num_verts = verts_temp.size();
  auto vert_gid_array = ViewCArray <int> ( &verts_temp[0], num_verts );
  
  

  //for(int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
  for (int index = 0; index < num_verts; index++){
    int vert_gid = vert_gid_array(index); 

    //for (int vertex = 0; vertex < num_basis; vertex++){
    //  int node_lid = elem.vert_node_map(vertex);
    //  int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);
    int num_elems_in_vert = mesh.num_elems_in_node(vert_gid);

    
      real_t lumped_mass = 0.0;
      
      for (int elems_in_vert = 0; elems_in_vert < num_elems_in_vert; elems_in_vert++){
        int elems_in_vert_gid = mesh.elems_in_node(vert_gid, elems_in_vert);
        for (int basis_id = 0; basis_id < num_basis; basis_id++){        
          for (int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++){
	    int gauss_gid = mesh.gauss_in_elem(elems_in_vert_gid, gauss_lid);
	    lumped_mass += ref_elem.ref_nodal_basis(gauss_lid, vertex)
		           * ref_elem.ref_node_g_weights(gauss_lid)
			   * mesh.gauss_pt_det_j(gauss_gid);
	  }// end loop over gauss_lid
	}// end loop over basis_id
      }// end loop over elems_in_node_lid
      
      int nodal_res_size = num_elems_in_vert*num_dim;
      real_t nodal_res_a[nodal_res_size];
      for (int i = 0; i < nodal_res_size; i++) nodal_res_a[i] = 0.0;
      auto nodal_res = ViewCArray <real_t> ( &nodal_res_a[0], num_elems_in_vert, num_dim);

      real_t sum_nodal_res_a[num_dim];
      for (int i = 0; i < num_dim; i++) sum_nodal_res_a[i] = 0.0;
      auto sum_nodal_res = ViewCArray <real_t> ( &sum_nodal_res_a[0], num_dim);

      for (int elems_in_vert = 0; elems_in_vert < num_elems_in_vert; elems_in_vert++){

        int elems_in_vert_gid = mesh.elems_in_node(vert_gid, elems_in_vert);
        
	real_t res_mass_a[num_basis];
        for (int i = 0; i < num_basis; i++) res_mass_a[i] = 0.0;
	auto res_mass = ViewCArray <real_t> ( &res_mass_a[0], num_basis );

        // Create variable for \sum_{q} M_{pq} \delta u^r_q //      
        real_t Mv_a[num_dim];
        for (int i = 0; i < num_dim; i++) Mv_a[i] = 0.0;
	auto Mv = ViewCArray <real_t> ( &Mv_a[0], num_dim );

        // Create CArray to store volume integral over cells in elem of force at each sub time //
	int const force_size = num_dim*num_correction_steps;
        real_t force_a[force_size];
	for (int i = 0; i < force_size; i++) force_a[i] = 0.0;
        auto force = ViewCArray <real_t> (&force_a[0], num_dim, num_correction_steps);

        for (int index = 0; index < num_basis; index++){
          for (int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++){
	    int gauss_gid = mesh.gauss_in_elem(elems_in_vert_gid, gauss_lid);
	    res_mass(index) += ref_elem.ref_nodal_basis(gauss_lid, vertex) 
                                 * ref_elem.ref_nodal_basis(gauss_lid, index)
                                 * ref_elem.ref_node_g_weights(gauss_lid)
		                 * mesh.gauss_pt_det_j(gauss_gid);
	  }// end loop over gauss_lid
        }// end loop over index
     
        /* 
        for (int index = 0; index < num_basis; index++){
          std::cout << " mass vec in res at j = " << index << " is " << res_mass(index) << std::endl;
	}
        */

        real_t volume_int_a[ num_dim ];
        auto volume_int = ViewCArray <real_t> (&volume_int_a[0], num_dim);
        for (int i = 0; i < num_dim; i++) volume_int(i) = 0.0;
      
        // Volume integral //
        for (int dim = 0; dim < num_dim; dim++){
	  for (int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++){
	    int gauss_gid = mesh.gauss_in_elem(elems_in_vert_gid, gauss_lid);
	    real_t J_inv_dot_grad_phi = 0.0;
	    for (int k = 0; k < num_dim; k++){
	      J_inv_dot_grad_phi += mesh.gauss_pt_jacobian_inverse(gauss_gid, k, dim)*ref_elem.ref_nodal_gradient(gauss_lid,vertex,k);
	    }// end loop over k
	    volume_int(dim) += mat_pt.pressure(gauss_lid)
	        	       * J_inv_dot_grad_phi
			       * ref_elem.ref_node_g_weights(gauss_lid)
			       * mesh.gauss_pt_det_j(gauss_gid);
	  }// end loop over gauss_lid
        }// end loop over dim

        real_t surface_int_a[ num_dim ];
        auto surface_int = ViewCArray <real_t> (&surface_int_a[0], num_dim);
        for (int i = 0; i < num_dim; i++) surface_int(i) = 0.0;
      
        // Surface Integral //
	// compute patch normal
	int const normal_size = num_dim*num_dim*mesh.num_patches_in_elem();
	real_t patch_normal_a[normal_size];
	for (int i = 0; i < normal_size; i++) patch_normal_a[i] = 0.0;
	auto patch_normal = ViewCArray <real_t> ( &patch_normal_a[0], num_dim, num_dim, mesh.num_patches_in_elem());

        for (int gauss_patch_lid = 0; gauss_patch_lid < mesh.num_patches_in_elem(); gauss_patch_lid++){
          int gauss_patch_gid = mesh.gauss_patch_pt_in_elem( elems_in_vert_gid, gauss_patch_lid );
            
	  patch_normal(0,0, gauss_patch_lid) = mesh.gauss_patch_pt_jacobian(gauss_patch_gid, 1, 1)*mesh.gauss_patch_pt_jacobian(gauss_patch_gid, 2,2) - mesh.gauss_patch_pt_jacobian(gauss_patch_gid, 2, 1)*mesh.gauss_patch_pt_jacobian(gauss_patch_gid,1,2);
          patch_normal(0,1, gauss_patch_lid) = mesh.gauss_patch_pt_jacobian(gauss_patch_gid, 2, 1)*mesh.gauss_patch_pt_jacobian(gauss_patch_gid, 0, 2) - mesh.gauss_patch_pt_jacobian(gauss_patch_gid, 0,1)*mesh.gauss_patch_pt_jacobian(gauss_patch_gid, 2, 2);
          patch_normal(0,2, gauss_patch_lid) = mesh.gauss_patch_pt_jacobian(gauss_patch_gid,0,1)*mesh.gauss_patch_pt_jacobian(gauss_patch_gid,1,2) - mesh.gauss_patch_pt_jacobian( gauss_patch_gid,1,1)*mesh.gauss_patch_pt_jacobian(gauss_patch_gid, 0,2);

	  patch_normal(1,0, gauss_patch_lid) = mesh.gauss_patch_pt_jacobian(gauss_patch_gid, 1, 0)*mesh.gauss_patch_pt_jacobian(gauss_patch_gid, 2,2) - mesh.gauss_patch_pt_jacobian(gauss_patch_gid, 2, 0)*mesh.gauss_patch_pt_jacobian(gauss_patch_gid,1,2);
          patch_normal(1,1, gauss_patch_lid) = mesh.gauss_patch_pt_jacobian(gauss_patch_gid, 2, 0)*mesh.gauss_patch_pt_jacobian(gauss_patch_gid, 0, 2) - mesh.gauss_patch_pt_jacobian(gauss_patch_gid, 0,0)*mesh.gauss_patch_pt_jacobian(gauss_patch_gid, 2, 2);
          patch_normal(1,2, gauss_patch_lid) = mesh.gauss_patch_pt_jacobian(gauss_patch_gid,0,0)*mesh.gauss_patch_pt_jacobian(gauss_patch_gid,1,2) - mesh.gauss_patch_pt_jacobian(gauss_patch_gid,1,0)*mesh.gauss_patch_pt_jacobian(gauss_patch_gid, 0,2);

	  patch_normal(2,0, gauss_patch_lid) = mesh.gauss_patch_pt_jacobian(gauss_patch_gid, 1, 1)*mesh.gauss_patch_pt_jacobian(gauss_patch_gid, 2,0) - mesh.gauss_patch_pt_jacobian(gauss_patch_gid, 2, 1)*mesh.gauss_patch_pt_jacobian(gauss_patch_gid,1,0);
          patch_normal(2,1, gauss_patch_lid) = mesh.gauss_patch_pt_jacobian(gauss_patch_gid, 2, 1)*mesh.gauss_patch_pt_jacobian(gauss_patch_gid, 0, 0) - mesh.gauss_patch_pt_jacobian(gauss_patch_gid, 0,1)*mesh.gauss_patch_pt_jacobian(gauss_patch_gid, 2, 0);
          patch_normal(2,2, gauss_patch_lid) = mesh.gauss_patch_pt_jacobian(gauss_patch_gid,0,1)*mesh.gauss_patch_pt_jacobian(gauss_patch_gid,1,0) - mesh.gauss_patch_pt_jacobian(gauss_patch_gid,1,1)*mesh.gauss_patch_pt_jacobian(gauss_patch_gid, 0,0);

	}// end loop over patch_lid

	// normalize
	for (int dim = 0; dim < num_dim; dim++){
	  for( int patch_lid = 0; patch_lid < mesh.num_patches_in_elem(); patch_lid++){
	    patch_normal(dim, 0, patch_lid) = patch_normal(dim, 0, patch_lid)/sqrt(patch_normal(dim, 0,patch_lid)*patch_normal(dim,0,patch_lid) + patch_normal(dim, 1,patch_lid)*patch_normal(dim,1,patch_lid) + patch_normal(dim, 2,patch_lid)*patch_normal(dim,2,patch_lid) ); 

	    patch_normal(dim, 1, patch_lid) = patch_normal(dim, 1, patch_lid)/sqrt(patch_normal(dim, 0,patch_lid)*patch_normal(dim,0,patch_lid) + patch_normal(dim, 1,patch_lid)*patch_normal(dim,1,patch_lid) + patch_normal(dim, 2,patch_lid)*patch_normal(dim,2,patch_lid) );  

	    patch_normal(dim, 2, patch_lid) = patch_normal(dim, 2, patch_lid)/sqrt(patch_normal(dim, 0,patch_lid)*patch_normal(dim,0,patch_lid) + patch_normal(dim, 1,patch_lid)*patch_normal(dim,1,patch_lid) + patch_normal(dim, 2,patch_lid)*patch_normal(dim,2,patch_lid) );  
	  }
	}

        //compute surface determinant //
	int surface_jacobian_size = num_dim*num_dim*num_dim*mesh.num_patches_in_elem();
	real_t surface_jacobian_a[surface_jacobian_size];
	for (int i = 0; i < surface_jacobian_size; i++) surface_jacobian_a[i] = 0.0;
	auto surface_jacobian = ViewCArray <real_t> ( &surface_jacobian_a[0], mesh.num_patches_in_elem(), num_dim, num_dim, num_dim);
        

	for (int dim = 0; dim < num_dim; dim++){
	  for (int patch_lid = 0; patch_lid < mesh.num_patches_in_elem(); patch_lid++){
            int gauss_patch_gid = mesh.gauss_patch_pt_in_elem(elems_in_vert_gid, patch_lid);
            if ( dim == 0 ){
	      for (int j = 1; j < num_dim; j++){
		int k = 0;
		for (int i = 0; i < num_dim; i++){
	          surface_jacobian(dim, patch_lid, i, k) = mesh.gauss_patch_pt_jacobian(gauss_patch_gid, i,j);    
		}
		k++;
	      }
	      for (int i = 0; i < num_dim; i++){
	        surface_jacobian(dim, patch_lid, i, 2) = patch_normal(dim, i, patch_lid);
	      }
	    }
	    else if (dim == 1){
	      for (int j = 0; j < num_dim; j = j+2){
		int k = 0;
		for (int i = 0; i < num_dim; i++){
	          surface_jacobian(dim, patch_lid, i, k) = mesh.gauss_patch_pt_jacobian(gauss_patch_gid, i,j);   
		}
		k++;
	      }
	      for (int i = 0; i < num_dim; i++){
	        surface_jacobian(dim, patch_lid, i, 2) = patch_normal(dim, i, patch_lid);
	      }
	    }
	    else if ( dim == 2){
	      for (int j = 0; j < 2; j++){
		for (int i = 0; i < num_dim; i++){
	          surface_jacobian(dim, patch_lid, i, j) = mesh.gauss_patch_pt_jacobian(gauss_patch_gid, i,j);    
		}
	      }
	      for (int i = 0; i < num_dim; i++){
	        surface_jacobian(dim, patch_lid, i, 2) = patch_normal(dim, i, patch_lid);
	      }
	    };
	  }
	}	
	
	int det_surf_jac_size = num_dim*mesh.num_patches_in_elem();
	real_t det_surf_jacobian_a[det_surf_jac_size];
	for (int i = 0; i < det_surf_jac_size; i++) det_surf_jacobian_a[i] = 0.0;
        auto det_surf_jacobian = ViewCArray <real_t> ( &det_surf_jacobian_a[0], num_dim, mesh.num_patches_in_elem() );
        
        // compute surface jacobian determinant
	for (int dim = 0; dim < num_dim; dim++){
	  for (int patch_lid = 0; patch_lid < mesh.num_patches_in_elem(); patch_lid++){
	    det_surf_jacobian(dim, patch_lid)  = surface_jacobian(dim, patch_lid, 0, 0) * ( surface_jacobian(dim, patch_lid,1, 1) *  surface_jacobian(dim, patch_lid,2, 2) -  surface_jacobian(dim, patch_lid,2, 1) *  surface_jacobian(dim, patch_lid,1, 2)) -
                 surface_jacobian(dim, patch_lid,0, 1) * ( surface_jacobian(dim, patch_lid,1, 0) *  surface_jacobian(dim, patch_lid,2, 2) -  surface_jacobian(dim, patch_lid,1, 2) *  surface_jacobian(dim, patch_lid,2, 0)) +
                 surface_jacobian(dim, patch_lid,0, 2) * ( surface_jacobian(dim, patch_lid,1, 0) *  surface_jacobian(dim, patch_lid,2, 1) -  surface_jacobian(dim, patch_lid,1, 1) *  surface_jacobian(dim, patch_lid,2, 0));
	  }
	}


        for (int dim = 0; dim < num_dim; dim++){
          for(int patch_gauss_lid = 0; patch_gauss_lid < mesh.num_patches_in_elem(); patch_gauss_lid++){

            int patch_gid = mesh.gauss_patch_pt_in_elem(elems_in_vert_gid, patch_gauss_lid);
	    surface_int(dim) += ref_elem.ref_patch_basis(patch_gauss_lid, vertex)
				  * mat_pt.pressure(patch_gauss_lid)
				  * ref_elem.ref_patch_g_weights(patch_gauss_lid)
				  * det_surf_jacobian(dim,patch_gauss_lid);// mesh.gauss_patch_pt_det_j(patch_gid);
	  }// end loop over gauss_lid
        }// end loop over dim
       
        for (int dim = 0; dim < num_dim; dim++){
          force(dim, current) = surface_int(dim) - volume_int(dim);
	  //std::cout << force(dim) << std::endl;
        }
        
	
        for (int dim = 0; dim < num_dim; dim++){
          for (int index = 0; index < num_basis; index++){
            //std::cout << vel_r(index,dim) - vel(index,dim) << std::endl;
	    Mv(dim) += res_mass( index )*(elem_state.BV_vel_coeffs( current, elems_in_vert_gid,index,dim) - elem_state.BV_vel_coeffs( 0, elems_in_vert_gid,index,dim));
          }// end loop over index
        }// end loop over dim
       

       /* 
        for (int dim = 0; dim < num_dim; dim++){
          std::cout<< "Mv/dt at dim " << dim << " and correction_step " << t_step << " is " << Mv(dim) << std::endl;
	  std::cout << "force at dim " << dim << " is  " << dt*force(dim,current) << std::endl;
        }
        */

	for (int dim = 0; dim < num_dim; dim++){
	  // Assign values to nodal res
          nodal_res( elems_in_vert, dim ) = Mv(dim)/dt + 0.5*( force( dim, 0 ) + force( dim, current) );
        }// end loop over dim
	

      }// end loop over elements in vertex


      
      for (int dim = 0;  dim < num_dim; dim++){  
        for (int elems_in_vert = 0; elems_in_vert< num_elems_in_vert; elems_in_vert++){
          sum_nodal_res(dim) += nodal_res( elems_in_vert, dim);
        }
      }
      for (int dim = 0; dim < num_dim; dim++){
        elem_state.BV_vel_coeffs( update, elem_gid, vertex, dim ) = 0.0;
      }

      for (int dim = 0 ; dim < num_dim; dim++){
          elem_state.BV_vel_coeffs( update, elem_gid, vertex, dim ) = elem_state.BV_vel_coeffs(current, elem_gid, vertex, dim) - (dt/lumped_mass)*sum_nodal_res(dim);
       	  //std::cout << "--- dim ---"<< std::endl;
	  //std::cout << "  "<< dim << "  "<< std::endl;
	  //std::cout << " sum_nodal_res = "<< sum_nodal_res( dim ) << std::endl;
	  //std::cout << " vel_coeff = "<< elem_state.BV_vel_coeffs( update, elem_gid, vertex, dim ) << std::endl;
      }

    }// end loop over vertex

    if ( update == num_correction_steps ){
      for (int node_lid = 0; node_lid < mesh.num_nodes_in_elem(); node_lid++){
        int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);

        for (int dim = 0; dim < mesh.num_dim(); dim++){
          node.vel(1, node_gid, dim) = 0.0;
        }// end loop over dim

        for (int dim = 0; dim < mesh.num_dim(); dim++){
          for (int vert = 0; vert < ref_elem.num_basis(); vert++){
            node.vel( 1, node_gid, dim ) += ref_elem.ref_nodal_basis( node_lid, vert ) * elem_state.BV_vel_coeffs( num_correction_steps, elem_gid, vert, dim );
          }// end loop over vertex
        }// end loop over dim

        for (int dim = 0; dim <  mesh.num_dim(); dim++){
          for (int vert = 0; vert < ref_elem.num_basis(); vert++){
	    elem_state.BV_vel_coeffs( 0, elem_gid, vert, dim ) = elem_state.BV_vel_coeffs( num_correction_steps, elem_gid, vert, dim );
	  }
        }
/*
        //// print statements ///
        for (int elem_id = 0; elem_id < mesh.num_elems(); elem_id++){
          std::cout << " ------- elem id ------- " << std::endl;
          std::cout << elem_id << std::endl;
          std::cout << " ----------------------- " << std::endl;
          for (int dim = 0; dim < mesh.num_dim(); dim++){
            std::cout << "-------- dim ------" <<std::endl;
            std::cout << dim << std::endl;
            for (int basis_id = 0; basis_id < ref_elem.num_basis(); basis_id++){
              std::cout << elem_state.BV_vel_coeffs( 0, elem_id, basis_id, dim ) << ", ";  
            }
            std::cout<<std::endl;
          }
          std::cout << " ----------------------- " << std::endl;
          std::cout << " ----------------------- " << std::endl;
        }
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
*/
      }// end loop over node_lid
     }// end if

  }// end loop over elements

}// end get_nodal_res



  /*    
  // Create CArray to store time integral of force and artificial viscosity
  real_t time_int_a[num_dim];
  for (int i = 0; i < num_dim; i++) time_int_a[i] = 0.0;
  auto time_integral = ViewCArray <real_t> (&time_int_a[0], num_dim);  
  */
   
   /*    	
      for (int dim = 0; dim < num_dim; dim++){

        real_t temp_dt = sub_dt;

       for (int prev_times = 0; prev_times <= current; prev_times++){    
         real_t point = ( (time_points[prev_times] + 1)* 0.5 * temp_dt );// <-- scale quadrature point

	  time_integral( dim ) = 0.5*temp_dt // <-- scale weights
		                  * time_weights[prev_times]
              		          *( force( dim, prev_times ) );
				     
          temp_dt += temp_dt;
        }// end loop over prev_times
      }// end loop over dim
   */


/*        
  // Create CArray for vel_bar used in artificial viscosity //
  real_t vel_bar_a[num_dim*t_step];
  for (int i = 0; i < num_dim*t_step; i++) vel_bar_a[i] = 0.0;
  auto vel_bar = ViewCArray <real_t> (&vel_bar_a[0], num_dim, t_step);

  // Create CArray for Q //
  real_t Q_a[num_dim*t_step];
  for (int i = 0; i < num_dim*t_step; i++) Q_a[i] = 0.0;
  auto Q = ViewCArray <real_t> ( &Q_a[0], num_dim, t_step);
*/
  /* 
    //for(int cell_lid = 0; cell_lid < mesh.num_cells_in_elem(); cell_lid++){
    //  int cell_gid = mesh.cells_in_elem(elem_gid, cell_lid);
    //for(int node_lid = 0; node_lid < mesh.num_nodes_in_cell(); node_lid++){
    for ( int vertex = 0; vertex < ref_elem.num_basis(); vertex++){
	int node_lid = ref_elem.vert_node_map(vertex);
        int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);

        // Create a view of vel and vel_n //
        auto vel = ViewCArray <real_t> (&elem_state.BV_vel_coeffs( current, elem_gid, vertex, 0 ), num_dim);
        auto vel_n = ViewCArray <real_t> (&elem_state.BV_vel_coeffs( 0, elem_gid ,vertex, 0 ), num_dim);
        
        // Initialize Q, sigma, vel_bar volume integral of "force" and time_integral //

        for (int prev_times = 0; prev_times <= current; prev_times++){
          for (int dim = 0; dim < num_dim; dim++){
            Q(dim, prev_times) = 0.0;
            vel_bar(dim, prev_times) = 0.0;
            force_cell_volume(dim, prev_times) = 0.0;
            for (int i = 0; i < num_dim; i++){
              sigma(i,dim,prev_times) = 0.0;
            }// end loop over i
            time_integral(dim) = 0.0;
	  }// end loop over dim
        }// end loop over prev_times

        
        // Compute vel_bar //
          
        for (int prev_times = 0; prev_times <= current; prev_times ++){
	  for (int dim_j = 0; dim_j < num_dim; dim_j++){
            // Loop over node_lid in cells to compute vel_bar //
            //for (int node_lid_in_cell = 0; node_lid_in_cell < mesh.num_nodes_in_cell(); node_lid_in_cell++){
	    for (int vert = 0; vert <  ref_elem.num_basis(); vert++){
	      //int node_lid_in_cell = ref_elem.vert_node_map(vert);
              // Get node_gid for each node_lid in cell  //
              //int node_gid_from_cell = mesh.nodes_in_cell(cell_gid, node_lid_in_cell);
              vel_bar(dim_j, prev_times) += vel(prev_times, elem_gid, vertex, dim_j);
              //std::cout << " vel_bar = " << vel_bar(dim_j, prev_times) << std::endl;   
            }// end loop over nodes in cell
          }// end loop over dim_j for vel_bar
        }// end loop over prev_times for vel_bar
        

         real_t mass_a[num_basis*num_basis];
         auto res_mass = ViewCArray <real_t> (&mass_a[0], num_basis);

         for (int basis_n = 0; basis_n < num_basis; basis_n++){
	   for (int basis_m = 0; basis_m < num_basis; basis_m++){
             res_mass(basis_m, basis_n) = 0.0;
	   }// end loop over basis_m
         }// end loop over basis_n 
    
         for(int basis_m = 0; basis_m < num_basis; basis_m++){
           for(int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++){
             int gauss_gid = mesh.gauss_in_elem(elem_gid,gauss_lid);//cell(cells_in_node_gid, gauss_lid);
             mass_vec(basis_m, basis_n) += 1.0//cell_state.density(cell_gid)//cells_in_node_gid)
                                  *ref_elem.ref_nodal_basis(gauss_lid,basis_m)//cell_basis(gauss_lid, basis_m)
                                  *ref_elem.ref_nodal_basis(gauss_lid, vertex)//cell_basis(gauss_lid, node_lid)
                                  *mesh.gauss_pt_det_j(gauss_gid)//cell_pt_det_j(gauss_gid)
                                  *ref_elem.ref_node_g_weights(gauss_lid);//cell_g_weights(gauss_lid);
           } // end loop over gauss in element
           mass += mass_vec(basis_m);
         } // end loop over basis_m


        // artificial viscosity //
        
        for (int dim_j=0; dim_j < num_dim; dim_j++){ 
          for (int prev_times = 0; prev_times <= current; prev_times++){
            real_t alpha_a[3];
            real_t speed = sqrt( node.vel(prev_times, node_gid, 0)*node.vel(prev_times, node_gid, 0) 
			+ node.vel(prev_times, node_gid, 1)*node.vel(prev_times, node_gid, 1)
			+ node.vel(prev_times, node_gid, 2)*node.vel(prev_times, node_gid, 2) );
            alpha_a[0] = speed + cell_state.cs(cell_gid);
            alpha_a[1] = cell_state.cs(cell_gid);
            alpha_a[2] = speed - cell_state.cs(cell_gid);
        
	    real_t alpha = 0.0;


	    for (int i = 0; i < 3; i++){
	      real_t temp1 = alpha_a[0] > alpha_a[1] ? alpha_a[0] : alpha_a[1];
	      alpha = alpha_a[2] > temp1 ? alpha_a[2] : temp1;          	  
            }
        
           // std::cout << "a0 = " << alpha_a[0] << ", a1 = " << alpha_a[1] << ", a2 = " << alpha_a[2] << std::endl;
	   // std::cout << " alpha = " << alpha << std::endl;
		  
	    // Fill Q //
            // Q = alpha_k*(vel - vel_bar) //
            Q(dim_j, prev_times) = alpha
                                   *(node.vel(prev_times,node_gid,dim_j) - vel_bar(dim_j,prev_times));
            //std::cout <<  "Q at dim " << dim_j << " is " << Q(dim_j, prev_times) << std::endl;  
          }// end loop over prev_times for Q
        }// end loop over dim_j for Q 
       
        
 
       // Loop over cells in node //
       for (int cells_in_node_lid = 0; cells_in_node_lid < mesh.num_cells_in_node(node_gid); cells_in_node_lid++){
         // get cell_gid //
         int cells_in_node_gid = mesh.cells_in_node(node_gid, cells_in_node_lid);
       
         for (int prev_times = 0; prev_times <= current; prev_times++){
       	   for (int dim_j=0; dim_j < num_dim; dim_j++){
	     for (int dim_i = 0; dim_i < num_dim; dim_i++){
               // Fill sigma //
               sigma(dim_i, dim_j, prev_times) = cell_state.stress( 0, cells_in_node_gid, dim_i, dim_j);
               //std::cout << " sigma at dim_i = "<< dim_i << " and dim_j = "<< dim_j <<" is equal to " << sigma(dim_i,dim_j,prev_times) << std::endl;
             }//  end dim_j sigma
           } // end dim_i sigma
         } // end loop over prev_times        
         
	 real_t cell_force_a[num_dim*num_basis*t_step];
	 auto cell_force = ViewCArray <real_t> (&cell_force_a[0], num_dim, num_basis, t_step);
       }// end loop over cells_in_node_lid
       
       for (int dim_j=0; dim_j < num_dim; dim_j++){

         real_t temp_dt = sub_dt;

         for (int prev_times = 0; prev_times <= current; prev_times++){    
             // Begin time integration of force_cell_volume and Q //
             // int^{t^m}_{t^n}(Q^{r,m}_p + \int_{V_h}(\grad\varphi\cdot\sigma)dV)dt
           real_t point = ( time_points[prev_times]* 0.5 * temp_dt + 0.5*sub_dt );
     //std::cout << " point = " << point << std::endl;
           time_integral(dim_j) += 0.5 *temp_dt * time_weights[prev_times]
              		           *(force_cell_volume(dim_j,prev_times) + 0.0*Q(dim_j,prev_times))
				   *std::legendre(t_step, point); // *legendre::eval(t_step,point)<-- depending on the compiler, legendre namespace conflicts with std::legendre
				     
           temp_dt += temp_dt;
             
//	     std::cout << " pseudo_dt at step " << prev_times << " is " << temp_dt << std::endl;
//	     std::cout<< " force at time_step " << prev_times <<" and dim "<< dim_j  << " is = "<< force_cell_volume(dim_j,prev_times)*std::legendre(t_step,point) << std::endl;	
//	     std::cout<< " Q at time_step " << prev_times <<" and dim "<< dim_j  << " is = "<< Q(dim_j,prev_times) << std::endl;                    
//           std::cout << "time integral is = " << time_integral(dim_j) << std::endl;
         }// end loop over prev times for time_integral
       } //end loop over dim_j for time_integral

          
       // Store lo_res with node_gid and cell_gid //
       for (int dim = 0; dim < num_dim; dim++){
         elem_state.nodal_res( elem_gid, vertex, dim ) = mass*(vel(dim) - vel_n(dim)) + time_integral(dim);// cells_in_node_gid, dim) = mass*(vel(dim) - vel_n(dim)) + time_integral(dim);
       }// end loop over dim  

    }// end loop over vertex
  */

       /*  	 
       for (int prev_times = 0; prev_times <= current; prev_times++){
         for (int basis_id = 0; basis_id < ref_elem.num_basis(); basis_id++){
           for (int dim_k = 0; dim_k < num_dim; dim_k++){
             for (int dim_j=0; dim_j < num_dim; dim_j++){
               for (int dim_i = 0; dim_i < num_dim; dim_i++){
	         for (int cells_in_node_lid = 0; cells_in_nodes_lid < mesh.num_cells_in_node(node_gid); cell_in_node_lid++){
		   int cells_in_node_gid = mesh.cells_in_node(node_gid, cells_in_node_lid);
                   for (int gauss_cell_lid = 0; gauss_cell_lid < mesh.num_gauss_in_cell(); gauss_cell_lid++){
                 
                     force_cell_volume(dim_k, basis_id, prev_times) += mesh.gauss_cell_pt_jacobian_inverse(cells_in_node_gid, dim_i, dim_j)
                                                           *sigma(dim_j, dim_k, prev_times)
			                                   *ref_elem.ref_cell_gradient(gauss_cell_lid, basis_id, dim_k)
                                                           *ref_elem.ref_cell_g_weights(gauss_cell_lid)
                                                           *mesh.gauss_cell_pt_det_j(cells_in_node_gid);
                   
                   }// end loop over gauss_cell_lid
		 }// end loop over cells_in_node_lid
                 force_cell_volume(dim_k, basis_id, prev_times) += cell_force(dim_k, basis_id, prev_times);
               }// end dim_i for force_cell_volume
             }// end dim_j for force_cell_volume
               //std::cout << "volume integral of force = "<< force_cell_volume(dim_k,prev_times) << std::endl;
           }// end loop over dim_k
	 }// end loop over basis_id
       }// end loop over prev_times for force_cell_volume
       */ 





/*
  // Initialize nodal_res to zero //
#pragma omp simd 
  for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
    for (int vertex = 0; vertex < ref_elem.num_basis(); vertex++){
      int node_lid = elem.vert_node_map(vertex);
      int node_gid = mesh.nodes_in_elem(elem_gid, node_lid);

      for (int cell_lid = 0; cell_lid < mesh.num_cells_in_node(node_gid);  cell_lid++){
        int cell_gid = mesh.cells_in_node(node_gid, cell_lid);
        for (int dim = 0; dim < num_dim; dim++){
          elem_state.nodal_res( elem_gid, vertex, cell_gid, dim)  = 0.0;
        }// end loop over dim
      }// end loop over cell_lid
    }// end loop over vertex
  }// end loop over elem_gid

  real_t time_weights[num_correction_steps];
  for (int i = 0; i < num_correction_steps; i++) time_weights[i] = 0.0;
  
  real_t time_points[num_correction_steps];
  for (int i = 0; i < num_correction_steps; i++) time_points[i] = 0.0;

  if (current == 0){

    time_weights[0] = 0.5;

    time_points[0] = 0;

  }
  else if (current == 1){

    time_weights[0] = 1;
    time_weights[1] = 1;

    time_points[0] = -0.5773502691896257;
    time_points[1] = 0.5773502691896257;

  }
  else if (current == 2){

    time_weights[0] = 0.5555555555555556;
    time_weights[1] = 0.8888888888888888;
    time_weights[2] = 0.5555555555555556;
 
    time_points[0] = -0.7745966692414834;
    time_points[1] = 0.0000000000000000;
    time_points[2] = 0.7745966692414834; 
  }
  else if (t_step == 3){

    time_weights[0] = 0.3478548451374538;    
    time_weights[1] = 0.6521451548625461;
    time_weights[2] = 0.6521451548625461;
    time_weights[3] = 0.3478548451374538;

    time_points[0] = -0.8611363115940526;
    time_points[1] = -0.3399810435848563;
    time_points[2] = 0.3399810435848563; 
    time_points[3] = 0.8611363115940526;

  };

*/
