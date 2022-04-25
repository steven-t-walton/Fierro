// -----------------------------------------------------------------------------
// This code reads the mesh in different formats
//------------------------------------------------------------------------------
#include "mesh.h"
#include "state.h"
#include "variables.h"


// -----------------------------------------------------------------------------
// Reads an ensight .geo mesh file
//------------------------------------------------------------------------------
void read_mesh_ensight(char* MESH,
                       mesh_t &mesh,
                       node_t &node,
                       elem_t &elem,
                       corner_t &corner,
                       const size_t num_dims,
                       const size_t rk_num_bins){

    const size_t rk_level = 0;

	FILE *in;
    char ch;


    //read the mesh    WARNING: assumes a .geo file
    in = fopen(MESH,"r");  
    
    //skip 8 lines
    for (int j=1; j<=8;j++) {
        int i=0;
        while ((ch=(char)fgetc(in))!='\n') {
            i++;
            //printf("%c",ch);
        }
        //printf("\n");
    }  

    

    // --- Read in the nodes in the mesh ---
    
    int num_nodes = 0;
    
    fscanf(in,"%d",&num_nodes);
    printf("Num nodes read in %d\n" , num_nodes);
    
    // intialize node variables
    mesh.initialize_nodes(num_nodes);
    node.initialize(rk_num_bins, num_nodes, num_dims);
    

    // read the initial mesh coordinates
    // x-coords
    for (int node_id = 0; node_id < mesh.num_nodes; node_id++) {
        fscanf(in,"%le",&node.coords(rk_level,node_id, 0));
    }

    // y-coords
    for (int node_id = 0; node_id < mesh.num_nodes; node_id++) {
        fscanf(in,"%le",&node.coords(rk_level,node_id, 1));
    }  

    // z-coords
    for (int node_id = 0; node_id < mesh.num_nodes; node_id++) {
        fscanf(in,"%le",&node.coords(rk_level,node_id, 2));
    }

    
    ch = (char)fgetc(in);
    //printf("%c",ch);

    //skip 1 line
    for (int j=1; j<=1; j++) {
        int i=0;
        while ((ch=(char)fgetc(in))!='\n') {
            i++;
            //printf("%c",ch);
        }
        //printf("\n");
    }
    

    // --- read in the elements in the mesh ---
    int num_elem = 0;
    
    fscanf(in,"%d",&num_elem);
    printf("Num elements read in %d\n" , num_elem);

    // intialize elem variables
    mesh.initialize_elems(num_elem, num_dims);
    elem.initialize(rk_num_bins, num_nodes, num_dims);

    // for each cell read the list of associated nodes
    for (int elem_gid = 0; elem_gid < num_elem; elem_gid++) {
        for (int node_lid = 0; node_lid < 8; node_lid++){
            
            fscanf(in,"%d",&mesh.nodes_in_elem.host(elem_gid, node_lid));  // %d vs zu

            // shift to start node index space at 0
            mesh.nodes_in_elem.host(elem_gid, node_lid) -= 1;
        }
    }
    // update device side
    mesh.nodes_in_elem.update_device();
    
    
    // intialize corner variables
    int num_corners = num_elem*mesh.num_nodes_in_elem;
    mesh.initialize_corners(num_corners);
    corner.initialize(num_corners, num_dims);






    // Close mesh input file
    fclose(in);

}