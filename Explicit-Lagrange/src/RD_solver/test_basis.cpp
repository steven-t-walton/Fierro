#include <iostream>

#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"

using namespace utils;

void test_basis(){

	std::cout <<  "  " << std::endl;
	std::cout << " Basis Values " << std::endl;
  	for (int dof_i = 0; dof_i < ref_elem.num_basis(); dof_i++){
            for (int dof_j = 0; dof_j < ref_elem.num_basis(); dof_j++){
	    	int node_lid = ref_elem.vert_node_map(dof_j);
		double test_fcn = ref_elem.ref_nodal_basis(node_lid, dof_i);
		std::cout << test_fcn <<" , ";
	    }
	    std::cout<<std::endl;
	}

	std::cout <<  "  " << std::endl;
	std::cout << " Dual Basis Values " << std::endl;
	for (int dof_i = 0; dof_i < ref_elem.num_dual_basis(); dof_i++){
            for (int dof_j = 0; dof_j < ref_elem.num_dual_basis(); dof_j++){
	    	int node_lid = ref_elem.ref_dual_vert_node_map(dof_j);
		double test_fcn = ref_elem.BV_basis(node_lid, dof_i);
		std::cout << test_fcn <<" , ";
	    }
	    std::cout<<std::endl;
	}

/*	
	std::cout <<  "  " << std::endl;
	std::cout << " K mass matrix values " << std::endl;
  	for (int dof_i = 0; dof_i < ref_elem.num_basis(); dof_i++){
            for (int dof_j = 0; dof_j < ref_elem.num_basis(); dof_j++){
            	double test_fcn =0.0;
            	for (int dof_k = 0; dof_k < ref_elem.num_basis(); dof_k++){
			int node_lid = elem.vert_node_map(dof_k);
                        int gauss_gid = mesh.gauss_in_elem(4, node_lid);	    	
		        test_fcn += ref_elem.ref_nodal_basis(node_lid, dof_i)
				           *ref_elem.ref_nodal_basis(node_lid, dof_j)
					   *mesh.gauss_pt_det_j(gauss_gid)
					   *ref_elem.ref_node_g_weights(node_lid);
	    	}
			std::cout << test_fcn <<" , ";
	     }
	    	std::cout<<std::endl;
	}
*/

/*	
     for (int elem_gid =0; elem_gid < mesh.num_elems(); elem_gid++){
	std::cout <<  "  " << std::endl;
	std::cout << " T mass matrix values in elem "<< elem_gid <<" are "<< std::endl;
  	for (int dof_i = 0; dof_i < ref_elem.num_dual_basis(); dof_i++){
            for (int dof_j = 0; dof_j < ref_elem.num_dual_basis(); dof_j++){
            	double test_fcn =0.0;
            	for (int dof_k = 0; dof_k < ref_elem.num_dual_basis(); dof_k++){
			int node_lid = elem.dual_vert_node_map(dof_k);
                        int gauss_gid = mesh.gauss_in_elem(elem_gid, node_lid);	    	
		        test_fcn += ref_elem.ref_nodal_dual_basis(node_lid, dof_i)
				           *ref_elem.ref_nodal_dual_basis(node_lid, dof_j)
					   *mesh.gauss_pt_det_j(gauss_gid)
					   *ref_elem.ref_node_g_weights(node_lid);
	    	}
			std::cout << test_fcn <<" , ";
	     }
	    	std::cout<<std::endl;
	}
     }
*/
  
};
