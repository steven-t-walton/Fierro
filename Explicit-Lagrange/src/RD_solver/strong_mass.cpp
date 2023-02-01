#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"
#include<iostream>

using namespace utils;

void get_strong_mass(){
  for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
    for(int gauss_lid = 0; gauss_lid < mesh.num_gauss_in_elem(); gauss_lid++){
      
      int gauss_gid = mesh.gauss_in_elem(elem_gid, gauss_lid);
      real_t vol_gauss = ref_elem.ref_node_g_weights(gauss_lid)*mesh.gauss_pt_det_j(gauss_gid);
      mat_pt.density(gauss_gid) = mat_pt.mass(gauss_gid)/vol_gauss;
      //std::cout << mat_pt.density(gauss_gid) << std::endl;
    
    } // end loop gauss
  }// end loop over elem_gid
}// end get_strong_mass
