//! \file heat_eq_boundary.cpp

#include "fem/heat_eq/heat_eq_boundary.h"
using namespace SlaughterFEM;

// Base class constructor
HeatEqBoundaryBase :: HeatEqBoundaryBase(std::string the_type) : BoundaryBase(the_type){
	// Sets the function pointer to NULL, this is used by the HeatEqSystem
	// class to determine what to do, use this pointer if it is not
	// NULL or use the value member function.
	fptr = NULL; 
}

// Dirichlet condtion class constructor
HeatEqBoundaryDirichlet :: HeatEqBoundaryDirichlet() : HeatEqBoundaryBase("dirichlet"){
	T_constant = 0;
}

// The default behavior of the T member is to return a constant value
Number HeatEqBoundaryDirichlet :: T(const Point&, const Real){
	return T_constant;
}

// Neumann condition class constructor
HeatEqBoundaryNeumann :: HeatEqBoundaryNeumann() : HeatEqBoundaryBase("neumann"){
	q_constant = 0;
}

// The default behavior of the q member is to return a constant value
Number HeatEqBoundaryNeumann :: q(const Point&, const Real){
	return q_constant;
}

// Convection condition class constructor
HeatEqBoundaryConvection :: HeatEqBoundaryConvection() : HeatEqBoundaryBase("convection"){
	h_constant = 0;
	Tinf_constant = 0; 
}

// The default behavior is to return a constant value
Number HeatEqBoundaryConvection :: h(const Point&, const Real){ 
	return h_constant; 
}

// The default behavior is to return a constant value
Number HeatEqBoundaryConvection :: Tinf(const Point&, const Real){ 
	return Tinf_constant; 
}
