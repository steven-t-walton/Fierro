#include<iostream>

#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"

using namespace utils;

void get_position_rdh(int correction_step){

  int update = correction_step;

// Update position with v(t^{n+1}) //
#pragma omp simd  
  for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){

    for (int vert = 0; vert < ref_elem.num_basis(); vert++){
      
      // View of vel coeffs //
      auto vel_update = ViewCArray <real_t> ( &elem_state.BV_vel_coeffs( update, elem_gid, vert, 0), mesh.num_dim() );
      auto vel_n = ViewCArray <real_t> ( &elem_state.BV_vel_coeffs( 0, elem_gid, vert, 0), mesh.num_dim() );

      // View of pos coeffs //
      auto pos_update = ViewCArray <real_t> ( &elem_state.BV_pos_coeffs( update, elem_gid, vert, 0), mesh.num_dim() );
      auto pos_n = ViewCArray <real_t> ( &elem_state.BV_pos_coeffs( 0, elem_gid, vert, 0), mesh.num_dim() );
      
      for (int dim = 0; dim < mesh.num_dim(); dim++){
        pos_update( dim) = pos_n(dim) + 0.5*dt * ( vel_update( dim ) + vel_n(dim) );
      }// end loop over dim

    }// end loop over node_lid
  }// end loop over elem_gid
}// end get_position_rdh()
