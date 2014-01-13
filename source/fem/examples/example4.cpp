/** \example fem/examples/example4.cpp 
 *
 */

//! @cond example4 Doxygen ignores the following code

// Boost includes
#include <boost/function.hpp>
#include <boost/bind.hpp> 
#include <boost/math/constants/constants.hpp>
 
// My includes 
#include "fem/include.h"
using namespace SlaughterFEM;
 
// Standard library includes  
#include <iostream>
#include <algorithm>
#include <math.h>
#include <stdio.h>
#include <string>
using std::string;

// libMesh includes
#include "libmesh.h"
#include "mesh.h"
#include "mesh_generation.h"
#include "equation_systems.h"
#include "transient_system.h"
#include "explicit_system.h"
#include "analytic_function.h"
#include "exodusII_io.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

// Initialization function for enthalpy, h
void initial_enthalpy(DenseVector<Number>& output, const Point&, const Real){
	output(0) = 0;
}

void initial_velocity(DenseVector<Number>& output, const Point&, const Real){
	output(0) = 0;	// x-direction
	output(1) = 0;	// y-direction
}

// The main function
int main (int argc, char** argv){
	
	// Initialize libraries
    LibMeshInit init (argc, argv);
        
	// Generate a mesh
	Mesh mesh;
	MeshTools::Generation::build_square(mesh, 10, 10, 0., 1, 0., 1, QUAD8);
	mesh.all_second_order();
	
	// Create an equation system
	EquationSystems eq_sys(mesh);

/*
	// Define the system constants
	eq_sys.parameters.set<Real>("rho_s") = 1078;
	eq_sys.parameters.set<Real>("rho_l") = 1078;
	eq_sys.parameters.set<Real>("k_s")   = 0.393;
	eq_sys.parameters.set<Real>("k_l")   = 0.468;
	eq_sys.parameters.set<Real>("c_s")   = 1870;
	eq_sys.parameters.set<Real>("c_l")   = 3249;
	eq_sys.parameters.set<Real>("D")     = 4.8e-9;
	eq_sys.parameters.set<Real>("mu")    = 1.3e-3;
	eq_sys.parameters.set<Real>("hf")    = 3.138e5;
	eq_sys.parameters.set<Real>("K0")    = 5.56e-11;
	eq_sys.parameters.set<Real>("beta_T")= 3.832e-4;
	eq_sys.parameters.set<Real>("beta_c")= 0.257;
	eq_sys.parameters.set<Real>("Te")    = 257.75;
	eq_sys.parameters.set<Real>("Ce")    = 0.803;
	eq_sys.parameters.set<Real>("Tm")    = 633.59;
	eq_sys.parameters.set<Real>("k_p")   = 0.30;
	 
	// Define dimensionless groups (should be computed from above)
	eq_sys.parameters.set<Real>("Pr")    = 9.025;
	eq_sys.parameters.set<Real>("Le")    = 27.84;
	eq_sys.parameters.set<Real>("Da")    = 8.896e-8;
	eq_sys.parameters.set<Real>("Ra_T")  = 1.938e7;
	eq_sys.parameters.set<Real>("Ra_c")  = -2.514e7;
	eq_sys.parameters.set<Real>("m")  	 = 0.905;
	eq_sys.parameters.set<Real>("R_k")   = 0.840;
	eq_sys.parameters.set<Real>("R_c")   = 0.576;
*/
	eq_sys.parameters.set<Real>("dt") = 0.01; 
	
	// Add the momentum equation
	MomentumEq momentum(eq_sys); 
	
	// Create a variable linker for the velocity solution of momentum eq
	EqVariableLinker<TransientNonlinearImplicitSystem> velocity(eq_sys, momentum.get_system());
	velocity.set_variables(momentum.get_velocity_variables());
	
printf("Name = %s\n", momentum.get_system().name().c_str());
printf("Name = %s\n", velocity.name().c_str());

	// Add the energy equation
	EnergyEq energy(eq_sys, velocity); // Add EqVariableLinker as to constructor
	//energy.add_initial_function(initial_enthalpy);
	 
	// 
	VolumeAverageEqData(eq_sys, velocity);
	
	
	
	
	// Link the velocity vector variables from the MomentumEq
	//energy.set_velocity_system(momentum.get_system(), momentum.velocity_components);
	 	
	 
	 	 







}
//! @endcond
