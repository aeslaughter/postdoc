 /*! \file momentum_system.cpp
* \brief Source file for the volume averaging MomentumSystem class.
*/

// Include file for the MomentumSystem class
#include "fem/volume_average/momentum_system.h"
using namespace SlaughterFEM;

// The constructor  
MomentumSystem :: MomentumSystem(EquationSystems& es, const string& name, const unsigned int number) : ImplicitSystemBase(es, name, number){

	/*
	// The number of dimensions
	const MeshBase& mesh = get_mesh();
	const unsigned int dim = mesh.mesh_dimension();

	// Velocity variables
	vector<string> var;

	// Build a vector of velocity components
	if (dim == 2 || dim == 3){
		add_variable("vx", order);
		add_variable("vy", order);		
	} else {
		printf("ERROR: Only 2 and 3 dimensions are supported\n");
		exit(1);
	} 
	
	if (dim == 3){
		add_variable("vz", order);
	}
	*/ 
} 

VectorValue<Number> MomentumSystem :: velocity(const Point& p){
	// The number of dimensions
	const MeshBase& mesh = get_mesh();
	const unsigned int dim = mesh.mesh_dimension();
	
	// Collect the velocity vector components
	VectorValue<Number> v;
	for (unsigned int i = 0; i < dim; i++){
		v(i) = System::point_value(i, p);
	}
	
	// Return the vector containing desired variables
	return v;
}

// Initialization function
void MomentumSystem :: initialize(){
	
	// Call the base initialization function
	SystemBase::initialize();
	
	// Initially the old and current solutions are the same
	*this->old_local_solution = *this->current_local_solution;

	// This ThermoEq system is now initialized
	_initialized = true;
}

// Assembly function
void MomentumSystem :: assemble(){}
