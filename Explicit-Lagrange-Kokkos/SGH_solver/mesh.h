#ifndef MESH_H
#define MESH_H


#include "matar.h"
#include "state.h"

#define PI 3.141592653589793

/*
==========================
Nodal indexing convention
==========================

              K
              ^         J
              |        /
              |       /
              |      /
      7------------------6
     /|                 /|
    / |                / |
   /  |               /  |
  /   |              /   |
 /    |             /    |
4------------------5     |
|     |            |     | ----> I
|     |            |     |
|     |            |     |
|     |            |     |
|     3------------|-----2
|    /             |    /
|   /              |   /
|  /               |  /
| /                | /
|/                 |/
0------------------1

nodes are ordered for outward normal
patch 0: [0,4,7,3]  xi-minus dir
patch 1: [1,2,6,5]  xi-plus  dir
patch 2: [0,1,5,4]  eta-minus dir
patch 3: [2,3,7,6]  eta-plus  dir
patch 4: [0,3,2,1]  zeta-minus dir
patch 6: [4,5,6,7]  zeta-plus  dir

*/

// sort in ascending order using bubble sort
KOKKOS_INLINE_FUNCTION
void bubble_sort(size_t arr[], size_t num){
    
    for (size_t i=0; i<(num-1); i++){
        for (size_t j=0; j<(num-i-1); j++){
            
            if (arr[j]>arr[j+1]){
                size_t temp = arr[j];
                arr[j] = arr[j+1];
                arr[j+1] = temp;
            } // end if
            
        } // end for j
    } // end for i
} // end function


// mesh sizes and connectivity data structures
struct mesh_t {
    
    size_t num_dims;
    
    size_t num_nodes;
    
    size_t num_elems;
    size_t num_nodes_in_elem;
    size_t num_patches_in_elem;

    size_t num_corners;
    
    size_t num_patches;
    size_t num_bdy_patches;
    size_t num_bdy_sets;
    size_t num_nodes_in_patch;

    
    // corner ids in node
    RaggedRightArrayKokkos <size_t> corners_in_node;
    CArrayKokkos <size_t> num_corners_in_node;
    
    // elem ids in node
    RaggedRightArrayKokkos <size_t> elems_in_node;
    
    // node ids in elem
    DCArrayKokkos <size_t> nodes_in_elem;
    
    // corner ids in elem
    CArrayKokkos <size_t> corners_in_elem;
    
    // elem ids in elem
    RaggedRightArrayKokkos <size_t> elems_in_elem;
    CArrayKokkos <size_t> num_elems_in_elem;
    
    // patch ids in elem
    CArrayKokkos <size_t> patches_in_elem;
    
    // node ids in a patch
    CArrayKokkos <size_t> nodes_in_patch;
    
    // patch ids in bdy set
    DynamicRaggedRightArrayKokkos <size_t> bdy_patches_in_set;
    
    // bdy_patches
    CArrayKokkos <size_t> bdy_patches;
    
    // element ids in a patch
    CArrayKokkos <size_t> elems_in_patch;
    
    
    // initialization methods
    void initialize_nodes(const size_t num_nodes_inp)
    {
        num_nodes = num_nodes_inp;
    }; // end method
    
    
    // initialization methods
    void initialize_elems(const size_t num_elems_inp, const size_t num_dims_inp)
    {
        num_dims = num_dims_inp;
        num_nodes_in_elem = 1;
        for (int dim=0; dim<num_dims; dim++){
            num_nodes_in_elem *= 2;
        }
        num_elems = num_elems_inp;
        nodes_in_elem = DCArrayKokkos <size_t> (num_elems, num_nodes_in_elem);
        corners_in_elem = CArrayKokkos <size_t> (num_elems, num_nodes_in_elem);
    }; // end method
    
    
    // initialization methods
    void initialize_corners(const size_t num_corners_inp)
    {
        num_corners = num_corners_inp;
    }; // end method
    
    
    // build the corner mesh connectivity arrays
    void build_corner_connectivity(){
        
        num_corners_in_node = CArrayKokkos <size_t> (num_nodes); // stride sizes
        
        // initializing the number of corners (node-cell pair) to be zero
        FOR_ALL(node_gid, 0, num_nodes, {
            num_corners_in_node(node_gid) = 0;
        });
        
        
        for (size_t elem_gid = 0; elem_gid < num_elems; elem_gid++){
            FOR_ALL(node_lid, 0, num_nodes_in_elem,{
                
                // get the global_id of the node
                size_t node_gid = nodes_in_elem(elem_gid, node_lid);

                // increment the number of corners attached to this point
                num_corners_in_node(node_gid) = num_corners_in_node(node_gid) + 1;
                
            });  // end FOR_ALL over nodes in element
        } // end for elem_gid
        
        
        // the stride sizes are the num_corners_in_node at the node
        corners_in_node = RaggedRightArrayKokkos <size_t> (num_corners_in_node);

        CArrayKokkos <size_t> count_saved_corners_in_node(num_nodes);

        // reset num_corners to zero
        FOR_ALL(node_gid, 0, num_nodes, {
            count_saved_corners_in_node(node_gid) = 0;
        });
        
        
        // he elems_in_elem data type
        elems_in_node = RaggedRightArrayKokkos <size_t> (num_corners_in_node);
        
        // populate the elems connected to a node list and corners in a node
        for (size_t elem_gid = 0; elem_gid < num_elems; elem_gid++){
            FOR_ALL(node_lid, 0, num_nodes_in_elem, {
                
                // get the global_id of the node
                size_t node_gid = nodes_in_elem(elem_gid, node_lid);
                
                // the column index is the num corners saved
                size_t j = count_saved_corners_in_node(node_gid);

                // Save corner index to this node_gid
                size_t corner_gid = node_lid + elem_gid*num_nodes_in_elem;
                corners_in_node(node_gid, j) = corner_gid;
                
                elems_in_node(node_gid, j) = elem_gid; // save the elem_gid

                // Save corner index to element
                size_t corner_lid = node_lid;
                corners_in_elem(elem_gid, corner_lid) = corner_gid;

                // increment the number of corners saved to this node_gid
                count_saved_corners_in_node(node_gid) = count_saved_corners_in_node(node_gid) + 1;

            });  // end FOR_ALL over nodes in element
        } // end for elem_gid
     
    } // end of build_corner_connectivity
    
    
    // build elem connectivity arrays
    void build_elem_elem_connectivity(){
        
        // find the max number of elems around a node
        size_t max_num_elems_in_node;
        size_t max_num_lcl;
        REDUCE_MAX(node_gid, 0, num_nodes, max_num_lcl, {
            
            // num_corners_in_node = num_elems_in_node
            size_t max_num = num_corners_in_node(node_gid);
            
            if (max_num > max_num_lcl) max_num_lcl = max_num;
                    
        }, max_num_elems_in_node); // end parallel reduction on max
        Kokkos::fence();
        
        // a temporary ragged array to save the elems around an elem
        DynamicRaggedRightArrayKokkos <size_t> temp_elems_in_elem(num_nodes, num_nodes_in_elem*max_num_elems_in_node);
        
        num_elems_in_elem = CArrayKokkos <size_t> (num_elems);
        FOR_ALL(elem_gid, 0, num_elems, {
            num_elems_in_elem(elem_gid) = 0;
        });
        Kokkos::fence();
        
        // find and save neighboring elem_gids of an elem
        FOR_ALL (elem_gid, 0, num_elems, {
            for (int node_lid=0; node_lid<num_nodes_in_elem; node_lid++){
                
                // get the gid for the node
                size_t node_id = nodes_in_elem(elem_gid, node_lid);
                
                // loop over all elems connected to node_gid
                for (int elem_lid = 0; elem_lid < num_corners_in_node(node_id); elem_lid++){
                    
                    // get the global id for the neighboring elem
                    size_t neighbor_elem_gid = elems_in_node(node_id, elem_lid);
                    
                    // a flag to save (=1) or not (=0)
                    size_t save = 1;
                    
                    // a true neighbor_elem_id is not equal to elem_gid
                    if (neighbor_elem_gid == elem_gid ){
                        save = 0;  // don't save
                    } // end if
                    
                    // check to see if the neighbor_elem_gid has been saved already
                    size_t num_saved = temp_elems_in_elem.stride(elem_gid);
                    for (size_t i=0; i<num_saved; i++){
                        
                        if (neighbor_elem_gid == temp_elems_in_elem(elem_gid,i)){
                            save=0;   // don't save, it has been saved already
                        } // end if
                        
                    } // end for i
                    
                    if (save==1){
                        // save the neighboring elem_gid
                        temp_elems_in_elem(elem_gid, num_saved) = neighbor_elem_gid;
                        
                        // increment the number of neighboring elements saved
                        temp_elems_in_elem.stride(elem_gid)++;
                    } // end if save
                    
                } // end for elem_lid in a node

            }  // end for node_lid in an elem
            
            // save the actial stride size
            num_elems_in_elem(elem_gid) = temp_elems_in_elem.stride(elem_gid);
            
        }); // end FOR_ALL elems
        Kokkos::fence();
        
        
        // compress out the extra space in the temp_elems_in_elem
        elems_in_elem = RaggedRightArrayKokkos <size_t> (num_elems_in_elem);
        
        FOR_ALL (elem_gid, 0, num_elems, {
            for (size_t i=0; i<num_elems_in_elem(elem_gid); i++){
                elems_in_elem(elem_gid, i) = temp_elems_in_elem(elem_gid, i);
            } // end for i
        });  // end FOR_ALL elems
        Kokkos::fence();
        
    } // end of build_elem_elem_connectivity
        
    
    // build the patches
    void build_patch_connectivity(){
        
        // building patches
        //CArrayKokkos <size_t> hash_patch(num_nodes, num_nodes, num_nodes,2);  // very slow
        DViewCArrayKokkos <size_t> node_ordering_in_cell;  // node lids in a patch
        
        num_nodes_in_patch = 2*(num_dims-1);  // 2 (2D) or 4 (3D)
        num_patches_in_elem = 2*num_dims; // 4 (2D) or 6 (3D)
        
        
        if(num_dims == 3) {
            size_t node_lids_in_patch_in_elem_3D[24] =
               {0,4,7,3,
                1,2,6,5,
                0,1,5,4,
                2,3,7,6,
                0,3,2,1,
                4,5,6,7};
            node_ordering_in_cell = DViewCArrayKokkos <size_t> (&node_lids_in_patch_in_elem_3D[0],6,4);
        }
        else {
            //   y
            //   |
            // 4---3
            // |   |  -- x
            // 1---2
            //
            size_t node_lids_in_patch_in_elem_2D[8] =
               {1,4,
                3,2,
                1,2,
                4,3};
            node_ordering_in_cell = DViewCArrayKokkos <size_t> (&node_lids_in_patch_in_elem_2D[0],4,2);
        } // end if
        
        
        

        // hash array
        CArrayKokkos <long int> hash_arr;
        if(num_dims==2){
            hash_arr = CArrayKokkos <long int> (num_nodes + num_nodes*num_nodes);
        }
        else
        {
            hash_arr = CArrayKokkos <long int> (num_nodes + num_nodes*num_nodes +
                                                num_nodes*num_nodes*num_nodes);
        } // end if
        
        // save the hash keys
        CArrayKokkos <size_t> hash_keys_in_elem (num_elems, num_patches_in_elem);
        
        
        // step 1) initialize the hash_arr = -1 at the patch hash_key values
        for (size_t elem_gid = 0; elem_gid<num_elems; elem_gid++){
            
            FOR_ALL(patch_lid, 0, num_patches_in_elem, {
                
                size_t sorted_patch_nodes[num_nodes_in_patch];
                
                // first save the patch nodes
                for (size_t node_lid = 0; node_lid<num_nodes_in_patch; node_lid++){
                    sorted_patch_nodes[node_lid] = node_ordering_in_cell(patch_lid,node_lid);
                }  // end for node_lid
                
                // sort nodes from smallest to largest
                bubble_sort(sorted_patch_nodes, num_nodes_in_patch);
                
                size_t hash_key;
                if(num_dims==2) {
                    hash_key = sorted_patch_nodes[0] + num_nodes*sorted_patch_nodes[1];
                }
                else {
                    hash_key = sorted_patch_nodes[1] + num_nodes*sorted_patch_nodes[2] +
                               num_nodes*num_nodes*sorted_patch_nodes[3];  // 3 largest node values
                } // end if on dims
                
                // save hash_keys in the this elem
                hash_keys_in_elem(elem_gid,patch_lid) = hash_key;
                
                hash_arr(hash_key) = -1; // tag the hash_array as having a patch

            }); // end for patch_lid
            
        } // end for elem_gid
        
        
        // step 2) count the number of patches in the mesh
        RUN({
            // serial execution on a GPU
            size_t patch_gid = 0;
            for (size_t elem_gid=0; elem_gid<num_elems; elem_gid++){
                for (size_t patch_lid=0; patch_lid<num_patches_in_elem; patch_lid++) {
                    
                    size_t hash_key = hash_keys_in_elem(elem_gid,patch_lid);
                    
                    // check to see if it is a new patch
                    if (hash_arr(hash_key) == -1){
                        hash_arr(hash_key) = patch_gid; // save the patch_gid
                        patch_gid++;
                    } // end if a new patch
                    
                } // end for patch_lid
            } // end for elem_gid
            
            num_patches = patch_gid;
        }); // end RUN
        
        printf("num_patches = %zu \n", num_patches);
        
        // allocate memory for the patch structures in mesh_t
        patches_in_elem = CArrayKokkos <size_t> (num_elems, num_patches_in_elem);
        elems_in_patch =  CArrayKokkos <size_t> (num_patches, 2);
        nodes_in_patch = CArrayKokkos <size_t> (num_patches, num_nodes_in_patch);
        
        // a temporary variable to help populate patch structures
        CArrayKokkos <size_t> num_elems_in_patch_saved (num_patches);
        
        // initialize the number of elems in a patch saved to zero
        FOR_ALL(patch_gid, 0, num_patches, {
            num_elems_in_patch_saved(patch_gid) = 0;
        });
        
        
        // step 3) populate the patch data structures with indices
        FOR_ALL(elem_gid, 0, num_elems, {
            for (size_t patch_lid=0; patch_lid<num_patches_in_elem; patch_lid++) {
                
                // get the patch_gid
                size_t hash_key = hash_keys_in_elem(elem_gid,patch_lid);
                size_t patch_gid = hash_arr(hash_key);
                
                patches_in_elem(elem_gid, patch_lid) = patch_gid;
                
                // save the elem_gid to the patch
                size_t num_saved = num_elems_in_patch_saved(patch_gid);
                elems_in_patch(patch_gid, num_saved) = elem_gid;
                
                // save the nodes in the patch if it is a new patch
                if (num_saved == 0){
                    // save the patch nodes
                    for (size_t node_lid = 0; node_lid<num_nodes_in_patch; node_lid++){
                        nodes_in_patch(patch_gid,node_lid) = node_ordering_in_cell(patch_lid,node_lid);
                    }  // end for node_lid
                } // end if num_saved
                
            } // end for patch_lid
            
        }); // end FOR_ALL elem_gid
        
        
    } // end patch connectivity method
    
    
    void init_bdy_sets (size_t num_sets){
        
        if(num_sets == 0){
            printf("ERROR: number of boundary sets = 0, setting it = 1");
            num_sets = 1;
        }
        num_bdy_sets = num_sets;
        bdy_patches_in_set = DynamicRaggedRightArrayKokkos <size_t> (num_sets, num_bdy_patches);
    } // end of init_bdy_sets method


}; // end mesh_t


namespace region
{

    // for tagging boundary faces
    enum vol_tag
    {
        global = 0,     // tag every cell in the mesh
        box = 1,        // tag all cells inside a box
        cylinder = 2,   // tag all cells inside a cylinder
        sphere = 3      // tag all cells inside a sphere
    };

} // end of namespace


namespace init_conds
{
    
    // applying initial conditions
    enum init_velocity_conds
    {
        // uniform
        cartesian = 0,   // cart velocity
        radial = 1,      // radial in the (x,y) plane where x=r*cos(theta) and y=r*sin(theta)
        spherical = 2,   // spherical
    
        // linear variation
        radial_linear = 3,     // linear variation from 0,0,0
        spherical_linear = 4,   // linear variation from 0,0,0
    
        // vortical initial conditions
        tg_vortex = 5
    };
    
} // end of initial conditions namespace


// fill instructions
struct mat_fill_t {
    
    // type
    region::vol_tag volume; // 1 is global, 2 are planes, 3 is a sphere
    
    // material id
    size_t mat_id;
    
    // planes
    double x1;
    double x2;
    double y1;
    double y2;
    double z1;
    double z2;
    
    // radius
    double radius1;
    double radius2;

    
    // initial conditions
    init_conds::init_velocity_conds velocity;
    
    // velocity coefficients by component
    double u,v,w;
    
    // velocity magnitude for radial velocity initialization
    double speed;
    
    double sie;  // specific internal energy
    double den;  // density
};


namespace bdy
{
    
    // for tagging boundary faces
    enum bdy_tag
    {
        x_plane  = 0,   // tag an x-plane
        y_plane  = 1,   // tag an y-plane
        z_plane  = 2,   // tag an z-plane
        cylinder = 3,   // tag an cylindrical surface
        sphere   = 4,   // tag a spherical surface
        readFile = 5    // read from a file
    };
    
    
    
    // for enforcing boundary conditions
    enum bdy_hydro_conds
    {
        fixed = 0,          // zero velocity
        reflected = 1,      // reflected or wall condition
        velocity = 2,       // constant velocity
        pressure = 3,       // constant pressure
        acceleration = 4,   // constant acceleration
        contact = 5         // contact surface
    };

} // end of bdy namespace


// tag mesh points on bdy's and set the BC type
struct boundary_t {

    // tag surface type
    bdy::bdy_tag surface;    // 0=xplane, 1=yplane, 2=zplane, 3=cylinder, 4=sphere, 5=read file
    
    // tag surface value or radius
    real_t value;
    
    // BC type
    bdy::bdy_hydro_conds hydro_bc;
    
};


void read_mesh_ensight(char* MESH,
                       mesh_t &mesh,
                       node_t &node,
                       elem_t &elem,
                       corner_t &corner,
                       const size_t num_dims,
                       const size_t rk_num_bins);


void input(CArrayKokkos <material_t> &material,
           CArrayKokkos <mat_fill_t> &mat_fill,
           CArrayKokkos <boundary_t> &boundary,
           CArrayKokkos <double> &state_vars,
           size_t &num_materials,
           size_t &num_fills,
           size_t &num_boundaries,
           size_t &num_dims,
           size_t &num_state_vars);


KOKKOS_FUNCTION
void get_vol_hex(const DViewCArrayKokkos <double> &elem_vol,
                 const size_t elem_gid,
                 const DViewCArrayKokkos <double> &node_coords,
                 const mesh_t &mesh);


KOKKOS_FUNCTION
void get_bmatrix(const ViewCArrayKokkos <double> &B_matrix,
                 const size_t elem_gid,
                 const DViewCArrayKokkos <double> &node_coords,
                 const mesh_t &mesh);


void setup( const CArrayKokkos <material_t> &material,
            const CArrayKokkos <mat_fill_t> &mat_fill,
            const CArrayKokkos <boundary_t> &boundary,
            const mesh_t &mesh,
            const DViewCArrayKokkos <double> &node_coords,
            const DViewCArrayKokkos <double> &node_vel,
            const DViewCArrayKokkos <double> &node_mass,      
            const DViewCArrayKokkos <double> &elem_den,
            const DViewCArrayKokkos <double> &elem_pres,
            const DViewCArrayKokkos <double> &elem_stress,
            const DViewCArrayKokkos <double> &elem_sspd,       
            const DViewCArrayKokkos <double> &elem_sie,
            const DViewCArrayKokkos <double> &elem_vol,
            const DViewCArrayKokkos <double> &elem_mass,
            const DViewCArrayKokkos <size_t> &elem_mat_id,
            const DViewCArrayKokkos <double> &elem_statev,
            const CArrayKokkos <double> &state_vars,
            const size_t num_fills,
            const size_t rk_num_bins,
            const size_t num_bdy_sets);

void ensight( mesh_t &mesh,
              DViewCArrayKokkos <double> &node_coords,
              DViewCArrayKokkos <double> &node_vel,
              DViewCArrayKokkos <double> &node_mass,
              DViewCArrayKokkos <double> &elem_den,
              DViewCArrayKokkos <double> &elem_pres,
              DViewCArrayKokkos <double> &elem_stress,
              DViewCArrayKokkos <double> &elem_sspd, 
              DViewCArrayKokkos <double> &elem_sie,
              DViewCArrayKokkos <double> &elem_vol,
              DViewCArrayKokkos <double> &elem_mass,
              DViewCArrayKokkos <size_t> &elem_mat_id,
              CArray <double> &graphics_times,
              size_t graphics_id,
              double time_value);
#endif 