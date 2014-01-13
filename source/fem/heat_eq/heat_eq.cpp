/*! \file heat_eq.cpp
  * \brief Source code for HeatEq class.
  */
  
// Include the header file 
#include "fem/heat_eq/heat_eq.h"  
using namespace SlaughterFEM;

// Class constructor
HeatEq :: HeatEq(EquationSystems& es, const Order order, const FEFamily family) : 
	EquationBase(es, "TransientHeatEquation"){
	
	// Add the unknown variable "u" to the the heat equation
    this->system().add_variable("u", order, family);

    // Define the default values for the parameters
	this->system().set_constant<Real>("theta" ,0.5);
    this->system().set_constant<Real>("k", 1);
    this->system().set_constant<Real>("rho", 1);
    this->system().set_constant<Real>("cp", 1);
}

// Class destruction
HeatEq :: ~HeatEq(){}
