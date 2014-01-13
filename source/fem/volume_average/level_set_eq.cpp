/*! \file level_set_eq.cpp
  * \brief Source code for LevelSetEq class.
  */
  
// Include the header file 
#include "fem/volume_average/level_set_eq.h"  
using namespace SlaughterFEM;

// Class constructor
LevelSetEq :: LevelSetEq(EquationSystems& es, const Order order, const FEFamily family) : 
	EquationBase(es, "LevelSetEquation"){
	
	// Add the unknown variable
    this->system().add_variable("phi", order, family);
}
