//! \file boundary_base.cpp

#include "fem/common/boundary_base.h"
using namespace SlaughterFEM;

// Base class constructor
BoundaryBase :: BoundaryBase(std::string the_type) : type(the_type){
	// Sets the function pointer to NULL, this is used by the EqBase
	// class to determine what to do, use this pointer if it is not
	// NULL or use the value member.
	fptr = NULL; 
}
