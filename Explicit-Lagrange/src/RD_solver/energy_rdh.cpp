// -----------------------------------------------------------------------------
//  \rho de/dt + \tau : div( u ) = 0  
//------------------------------------------------------------------------------
#include "utilities.h"
#include "state.h"
#include "geometry.h"
#include "variables.h"

using namespace utils;

void get_energy_rdh( real_t sub_dt){

    // loop over the elements in the mesh
#pragma omp simd
    for (int elem_gid = 0; elem_gid < mesh.num_elems(); elem_gid++){
        for(int cell_lid = 0; cell_lid < mesh.num_cells_in_elem(); cell_lid++){

            // get the global ID for this cell
            int cell_gid = mesh.cells_in_elem(elem_gid, cell_lid);

            // Loop over the nodes in the cell
            for (int node_lid = 0; node_lid < mesh.num_nodes_in_cell(); node_lid++){

                // Get node global id for the local vertex id
                int node_gid = mesh.nodes_in_cell(elem_gid, node_lid); 

                // create view of velocities at current and previous step
                auto vel   = ViewCArray <real_t> (&node.vel( num_correction_steps, node_gid, 0), num_dim);
                auto vel_n = ViewCArray <real_t> (&node.vel(0, node_gid, 0), num_dim);
		
		auto vel_mid = CArray <real_t> (num_dim);

		for (int dim = 0; dim <  num_dim; dim++){
		  vel_mid(dim) = 0.5*(vel(dim) + vel_n(dim));
		}

           }

            real_t x = mesh.cell_coords(cell_gid, 0);
            real_t y = mesh.cell_coords(cell_gid, 1);

            real_t front = 3.14159265/(4.0*((7.0/5.0) - 1.0));  

           
            // Source term for TG problem 
            real_t source = front * (cos(3.0*3.14159265*x)*cos(3.14159265*y) 
                                    - cos(3.0*3.14159265*y)*(cos(3.14159265*x))); 


            source = source*mesh.cell_vol(cell_gid);


            // update the specific energy
            cell_state.ie(num_correction_steps, cell_gid) = 
                cell_state.ie(0, cell_gid) + sub_dt / cell_state.mass(cell_gid) * source;
                
            // loop over all points to find total kinetic energy 
            real_t ke = 0.0;
            for (int node_lid = 0; node_lid < mesh.num_nodes_in_cell(); node_lid++) {
                
                int node_gid = mesh.nodes_in_cell(cell_gid, node_lid);
                // create view into vertex velocity
                auto vel = ViewCArray <real_t> (&node.vel(num_correction_steps, node_gid, 0), num_dim);

                ke += 0.5 * node.mass(node_gid) * 
                    (vel(0)*vel(0) + vel(1)*vel(1) + vel(2)*vel(2));

            }


            cell_state.ke(num_correction_steps,cell_gid) = ke;


            cell_state.total_energy(num_correction_steps, cell_gid) = ke + cell_state.ie(num_correction_steps, cell_gid) + source;

           
        } // end loop over cells in element
    } // end loop over the elements
} // end subroutine
