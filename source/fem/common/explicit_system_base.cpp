/*! \file explicit_system_base.cpp
  * \brief Source code for the ExplicitSystemBase class.
  */
   
// Include the header file for this source   
#include "fem/common/explicit_system_base.h"
using namespace SlaughterFEM;

// The class constructor
ExplicitSystemBase :: ExplicitSystemBase(EquationSystems& es, const string& name, const unsigned int number) : 
			SystemBase<TransientExplicitSystem>(es, name, number){
	
	// Attach an initilization function for the heat equation
	//this->attach_init_object(*this);	
}

/*
// By default the operator() calls component
void ExplicitSystemBase :: operator()(const Point& p, const Real t, DenseVector<Number>& output){	
	// Re-size the output vector
	output.resize(this->n_vars());

	// Compute the values
	for (unsigned int i = 0; i < output.size(); i++){
		output(i) = component(i,p);
	}
}	


// Do not use the scalar output version of the operator()
Number ExplicitSystemBase :: operator() (const Point&, const Real){
	libmesh_not_implemented();
}

// init() function
void ExplicitSystemBase :: init(Real t_initial){		
			
	// Add a time variable and set system time
	set_constant<Real>("time", t_initial);

	// Initialize
	SystemBase::init();
}

// The initilization function
void ExplicitSystemBase :: initialize(){
	this->project_solution(this);
}
*/
	
/*	
// A function for projecting the nodal data at prescribed time	
void ExplicitSystemBase :: update_solution(Real t){

	// Update the system type
	this->time = t;

	// Step the system with time
	*this->old_local_solution = *this->current_local_solution;

	// Re-initilize the system that stores the data
	this->reinit();
}
*/
