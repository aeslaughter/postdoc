/*! \file my_mesh.h
\brief Header file for MyMesh.cpp.
*/

// Avoid multiple includes
#ifndef my_mesh_h
#define my_mesh_h

// Include the necessary libMesh headers
#include "libmesh.h"
#include "mesh.h"
#include "boundary_info.h"
#include "side.h"
#include "elem.h"

/*! A namespace for my FEM related functions and classes
 * 
 * The FEM library contains tools that integrate with libMesh, which is 
 * the basis of the FEM code associated with the programs
 * and classes listed in this documentation.
 */
namespace SlaughterFEM {
       
/*! \class MyMesh my_mesh.h "fem/my_mesh/my_mesh.h"
 * \brief Adds additional boundary identification behavior to libMesh::Mesh class 
 * 
 * add_boundary_id: 
 * This function searches every element for sides that 
 * are not shared, if they are not shared then they must lie on a boundary. 
 * These boundary sides are then tested to see if all of the nodes have 
 * the specified value for the desired coordinate. If so, then the 
 * boundary ID is added to the side object.
 * 
 * There are three ways to use this feature, example code for each is
 * included in the member function documentation. 
 * 
 * \ingroup FEMcommon
 */  
class MyMesh : public Mesh {
public:

    //! A function for identifing a boundary based on the nodal position of the sides.
    /*!
        \param id an integer ID to identify the boundary
        \param dir_idx A integer that represents the space dimension to test: 0, 1, and 2 for x, y, and z directions respectively
        \param value A libmesh Real number to test the nodal coordinate against

    This example first finds the elements with neighbors and then clears any existing boundary IDs. Then all boundary sides that lie completely on x = 0.0 are given the boundary ID of 1.
    \code
        mesh.find_neighbors();
        mesh.boundary_info->clear();
        mesh.add_boundary_id(1, 0, 0.0);
    \endcode
    
    Specifing dir_idx = -1 is a special case that inserts the boundary ID to all unset boundaries. In this case, value is not used.
    */
    void add_boundary_id(short int id, int dir_idx, Real value);
    
    //! An alternative form of add_boundary_id, that allows a text input of the coordinate dimension.  
    /*! 
        \param id an integer ID to identify the boundary
        \param dir A char that represents the space dimension to test: x, y, or z
        \param value A libmesh Real number to test the nodal coordinate against
        
    Example code:
    \code
        mesh.find_neighbors();
        mesh.boundary_info->clear();
        mesh.add_boundary_id(1, "x", 0.0);
    \endcode   
    */ 
    void add_boundary_id(short int id, const char* dir, Real value);  
    
    //! An alternative form of add_boundary_id, that adds the specified ID to all boundaries without an ID
    /*! 
    This version essentially acts to redifine all boundaries with a value of invalid_id to the specified value.
        \param id an integer ID to identify the boundary
        
    \code
        mesh.find_neighbors();
        mesh.boundary_info->clear();
        mesh.add_boundary_id(1);
    \endcode    
    */ 
    void add_boundary_id(short int id); 
    
}; // end MyMesh class

}  // end SlaughterFEM namespace

#endif
