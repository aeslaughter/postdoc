/*! \file my_mesh.cpp 
\brief Custom mesh features developed for libMesh Mesh class.
*/

// My FEM includes
#include "fem/common/my_mesh.h"     // header for this source code
using namespace SlaughterFEM;

// Standard library includes
#include <string.h>     
// Function for adding boundary ID based on nodal position
void MyMesh :: add_boundary_id(short int id, int dir_idx, Real value){  

    // Test that the test coordinate was entered correctly
    if (dir_idx < -1 || dir_idx > 2){
        printf("ERROR: direction index is out-of-bounds, specify 0, 1, or 2. (see MyMesh :: add_boundary_id in MyMesh.cpp)\n");
        exit(0);
    }

    // Define the local elements to loop through
    MeshBase::const_element_iterator       el     = this->active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = this->active_local_elements_end();
    
    // Loop through all of the elements using the element iterator
    for ( ; el != end_el ; ++el)
    {
        const Elem* elem = *el;	// a pointer to the current element

        // Test if the element is on the boundary
        if(elem->on_boundary() == true){
            
            // Loop through each of the sides of the boundary element
            for (unsigned int i = 0; i < elem->n_sides(); i++){
        
                // If the side has no neighbors (NULL), it is on the boundary
                if (elem->neighbor(i) == NULL){
                    AutoPtr<Elem> this_side = elem->side(i); // pointer to the current side
                
                    // Special case that adds ID to boundary if it doesn't already have one.
                    if(dir_idx == -1){
                        short int boundary_id = (this->boundary_info)->boundary_id(elem, i);
                        short int invalid_id = (this->boundary_info)->invalid_id;
                        
                        if(boundary_id == invalid_id){
                            (this->boundary_info)->add_side(elem, i, id);   
                        }
                        
                   
                    // General case where ID is added based on node position 
                    } else {
                
                        // Variable for storing nodal position test results
                        unsigned int testsum = 0; 
                    
                        // Loop through all the nodes on the current side
                        for(unsigned int j = 0; j < this_side->n_nodes(); j++){
                            Point p = this_side->point(j); 			   	// the current point

                            // Increment testsum if the point is located along the desired nodal coordinate
                            if(p(dir_idx) == value) ++testsum;			  
                        }
                    
                        // Add the boundary ID to this side of all the points of the current side were located along the desired nodal coordinate
                        if(testsum == this_side->n_nodes()){
                            (this->boundary_info)->add_side(elem, i, id);
                        } 	// end if testsum
                        
                    }   // end -1 if special case statement
                } 	// end if test of boundary element side
            }	// end loop of element sides
        }	// end test of boundary elements
    } 	// end loop of elements function
}	// end function

// Overloaded version that allows user to specify coordinate with text
void MyMesh :: add_boundary_id(short int id, const char* dir, Real value){

    // Storage location for numeric version of coordinate test direction
    int dir_idx; 

    // Test the text version of the coordinate directions and return its numeric value
    if(strcmp("x", dir) == 0 || strcmp("X", dir) == 0){
        dir_idx = 0;
    } else if (strcmp("y", dir) == 0 || strcmp("Y", dir) == 0){
        dir_idx = 1;
    } else if (strcmp("z", dir) == 0 || strcmp("Z", dir) == 0){
        dir_idx = 2;
    } else {
        printf("ERROR: the direction index was not recongnized, specify x, y, or z. (see MyMesh :: add_boundary_id in MyMesh.cpp)\n");
        exit(0);
    }

    // Call the main version of the function with the numeric coordinate
    this->add_boundary_id(id, dir_idx, value);	
}

// Function for assigning ID to un-identified boundaries
void MyMesh :: add_boundary_id(short int id){
    this->add_boundary_id(id, -1, 0);  
}
