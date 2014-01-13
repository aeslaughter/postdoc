/*! \file front_velocity_eq.cpp
  * \brief Source code for FrontVelocityEq class.
  */
  
// Include the header file 
#include "fem/volume_average/front_velocity_eq.h"  
using namespace SlaughterFEM;

// Class constructor
FrontVelocityEq :: FrontVelocityEq(EquationSystems& es, const Order order, const FEFamily family) : 
	EquationBase(es, "FrontVelocityEquation"){
	
	
	// Determine the dimension
	unsigned int dim = this->system().ndim();
	
	// Only 2D and 3D supported 
	if (dim < 2 || dim > 3){
		printf("ERROR: Only 2D and 3D support exixts for FrontVelocityEq class.\n");
		libmesh_error();
	}
	
	// Add x and y velocity variables
    this->system().add_variable("front_velocity_x", order, family);
    this->system().add_variable("front_velocity_y", order, family);
    
    // Add the z-direction velocity variable
    if (dim == 3){
		this->system().add_variable("front_velocity_z", order, family);
	}
}
