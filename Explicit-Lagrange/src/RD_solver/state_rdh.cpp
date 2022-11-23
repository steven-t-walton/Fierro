/* state.cpp */

#include<iostream>
#include<math.h>

#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"

#define PI 3.14159265

using namespace utils;

void get_state(){
// /*

  for( int cell_gid = 0; cell_gid < mesh.num_cells(); cell_gid++){
    //cell_properties(cell_gid);
    real_t elem_coords_x = 0.0;
    real_t elem_coords_y = 0.0;
    elem_coords_x = mesh.cell_coords(cell_gid,0);
    elem_coords_y = mesh.cell_coords(cell_gid,1);

    real_t source = 0.0;

    source = 3.141592653589/(4.0*(0.66666667))* ( cos(3.0*3.141592653589 * elem_coords_x) * cos( 3.141592653589 * elem_coords_y) - cos( 3.141592653589 * elem_coords_x ) * cos( 3.0*3.141592653589 * elem_coords_y ) ); 
    
    cell_state.density(cell_gid) = 1.0;
      
    cell_state.pressure(cell_gid) = 0.25*( cos(2.0*3.141592653589 * elem_coords_x )  + cos(2.0*3.141592653589 * elem_coords_y ) ) + 1.0;
   
    cell_state.ie(1, cell_gid) = source + cell_state.pressure(cell_gid)/(cell_state.density(cell_gid)*(0.66666667));

    real_t ie = 0.0;
    real_t ke = 0.0;
    track_rdh(ie,ke);
    real_t energy = ke;//+ie
    cell_state.total_energy(1,cell_gid) = cell_state.ie(1,cell_gid) + energy;
    
    cell_state.cs(cell_gid) =
        material[cell_state.mat_id(cell_gid)].eos_func(
            sspd,
            cell_state.mat_id(cell_gid),
            cell_state.density(cell_gid),
            cell_state.ie(1, cell_gid)
            );
  }
}
// */  
 /*
  real_t det_J0_a[mesh.num_cells()];
  auto det_J0 = ViewCArray(&det_J0_a[0], mesh.num_cells());
  
  real_t rho_0_a[mesh.num_cells()];
  auto rho_0 = ViewCArray(&rho_0_a[0], mesh.num_cells());
  for (int cell_gid = 0; cell_gid < mesh.num_cells(); cell_gid++){
    det_J0( cell_gid ) = 0.0;
    rho_0( cell_gid ) = 0.0;
  }

#pragma omp simd
  if ( cycle ==1 ){
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
    if ( cycle == 1 ){
      cell_state.density(cell_gid) = 1;//rho_0(cell_gid);
    }
    else {
      cell_state.density(cell_gid) = 1; //rho_0(cell_gid) * det_J0(cell_gid)/ mesh.gauss_cell_pt_det_j(cell_gid);
      //std::cout << " det J(t) = " << mesh.gauss_cell_pt_det_j(cell_gid) <<  std::endl;
      //std::cout <<  " density = " << cell_state.density(cell_gid) << std::endl;
    }; 
  }// end loop over cell_gid
  
  for( int cell_gid = 0; cell_gid  < mesh.num_cells(); cell_gid++){
      real_t elem_coords_x = mesh.cell_coords(cell_gid,0);
      real_t elem_coords_y = mesh.cell_coords(cell_gid,1);
      // pressure //
      cell_state.pressure(cell_gid) = 0.25*( cos(2.0*3.141592653589 * elem_coords_x ) ) + cos(2.0*3.141592653589 * elem_coords_y ) + 1.0;

      // internal energy //
      cell_state.ie(1, cell_gid) = cell_state.ie( 0 , cell_gid) + dt/cell_state.mass(cell_gid) * 3.141592653589/(4.0*7.0/5.0 - 4.0)* ( cos(3.0*3.141592653589 * elem_coords_x) * cos( 3.141592653589 * elem_coords_y) - cos( 3.141592653589 * elem_coords_x ) * cos( 3.0*3.141592653589 * elem_coords_y ) );	
     // }// end loop over cell_lid	      
   // }// end loop over node_lid
  //}// end loop over elem_gid
  }// end loop over cell_gid
  

  // calculate the sound speed
#pragma omp simd
  for (int cell_gid = 0; cell_gid < mesh.num_cells(); cell_gid++){
    cell_state.cs(cell_gid) =
        material[cell_state.mat_id(cell_gid)].eos_func(
            sspd,
            cell_state.mat_id(cell_gid),
            cell_state.density(cell_gid),
            cell_state.ie(1, cell_gid)
            );
  }
  */




