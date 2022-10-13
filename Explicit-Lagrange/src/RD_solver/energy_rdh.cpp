// -----------------------------------------------------------------------------
//  \rho de/dt + \tau : div( u ) = 0 
//
//  currently for
//
//  de/dt + p div( u ) = 0, 
//------------------------------------------------------------------------------
#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"

using namespace utils;

void get_energy_rdh( real_t sub_dt){
  
  // step 1: compute residual of energy equation
  // step 2: update energy
  // step 3: interpolate with order=p-1 basis functions
  
  

  for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
  
    
  
  }// end loop over elem_gid


} // end subroutine
