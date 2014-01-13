/*! \file front_velocity_system.cpp
  * \brief Source code for FrontVelocitySystem class.
  */
  
// Include the header file 
#include "fem/volume_average/front_velocity_system.h"
using namespace SlaughterFEM;

// Class constructor
FrontVelocitySystem :: FrontVelocitySystem(EquationSystems& es, const string& name, const unsigned int number) : SystemBase<TransientExplicitSystem>(es, name, number){
	
	// Set the function pointer to null
	//_velocity_ptr = NULL;
}


void FrontVelocitySystem :: add_velocity_function(boost::function<void (DenseVector<Number>&, const Point&, const Real)> func){
	_velocity_ptr = func;
}
		

void FrontVelocitySystem :: initialize(){
	
	// Test that the user has supplied a pointer to the velocity equation
	if (_velocity_ptr == NULL){
		printf("ERROR: A function pointer to the velocity equation must be supplied.\n");
		libmesh_error();
	}
	
	// Project the velocity equation
	this->project(0.0);

	// Update the initialization state
	_initialized = true;	
}

void FrontVelocitySystem :: update_solution(Real t){
	
	this->time = t;
	this->project(t);

}

void FrontVelocitySystem :: project(Real t){
	
	// Project the velocity equation
	MyAnalyticFunction<Number> func_object(_velocity_ptr);
	this->project_solution(&func_object);
}
