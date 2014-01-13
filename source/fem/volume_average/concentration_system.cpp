 /*! \file concentration_system.cpp
* \brief Source file for the volume averaging ConcentrationEq class.
*/

// Include file for the ConcentrationEq class
#include "fem/volume_average/concentration_system.h"
using namespace SlaughterFEM;

// The constructor  
ConcentrationSystem :: ConcentrationSystem(EquationSystems& es, const string& name, const unsigned int number) : ImplicitSystemBase(es, name, number){
	
	// Define the concentration variable in libMesh
	//vector<string> var;
	//var.push_back("C");
	//add_variable(var[0], order);
} 

// Assembly function
void ConcentrationSystem :: assemble(){}
