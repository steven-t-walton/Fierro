/* state.cpp */

#include<iostream>

#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"

using namespace utils;

// -----------------------------------------------------------------------------
// This function calculates the cell pressure, density, sound speed
//------------------------------------------------------------------------------
void get_state(int cycle, int correction_step){

  real_t det_J0_a[mesh.num_cells()];
  auto det_J0 = ViewCArray(&det_J0_a[0], mesh.num_cells());
  
  real_t rho_0_a[mesh.num_cells()];
  auto rho_0 = ViewCArray(&rho_0_a[0], mesh.num_cells());
  for (int cell_gid = 0; cell_gid < mesh.num_cells(); cell_gid++){
    det_J0( cell_gid ) = 0.0;
    rho_0( cell_gid ) = 0.0;
  }
#pragma omp simd
  if ( cycle == 1 and correction_step == 0){
    for(int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
      for(int cell_lid = 0; cell_lid < mesh.num_cells_in_elem(); cell_lid++){
	int cell_gid = mesh.cells_in_elem(elem_gid, cell_lid);

        det_J0(cell_gid) = mesh.gauss_cell_pt_det_j(cell_gid);
        //std::cout << " det_J0 = " << det_J0(cell_gid);

        rho_0(cell_gid) = cell_state.density(cell_gid);
        //std::cout << ";  rho_0 = " << rho_0(cell_gid) << std::endl;
      }
    }
  };

#pragma omp simd
  for (int cell_gid = 0; cell_gid < mesh.num_cells(); cell_gid++) {

    // calculate the density  
    if ( cycle == 1 and correction_step == 0){
      cell_state.density(cell_gid) = rho_0(cell_gid) * det_J0(cell_gid);
    }
    else {
      cell_state.density(cell_gid) = rho_0(cell_gid) * det_J0(cell_gid)/ mesh.gauss_cell_pt_det_j(cell_gid);
    }; 
    
// calculate the pressure
    cell_state.pressure(cell_gid) =
        material[cell_state.mat_id(cell_gid)].eos_func(
            p_of_de,
            cell_state.mat_id(cell_gid),
            cell_state.density(cell_gid),
            cell_state.ie(correction_step, cell_gid)
            );

    // calculate the sound speed
    cell_state.cs(cell_gid) =
        material[cell_state.mat_id(cell_gid)].eos_func(
            sspd,
            cell_state.mat_id(cell_gid),
            cell_state.density(cell_gid),
            cell_state.ie(correction_step, cell_gid)
            );
  }
}


