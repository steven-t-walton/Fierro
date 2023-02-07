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
  for( int cell_gid = 0; cell_gid < mesh.num_cells(); cell_gid++){
    cell_properties(cell_gid);
  }
}



