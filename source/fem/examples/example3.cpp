/** \example fem/examples/example3.cpp 
 * A test function for multiple variable nodal data with libmesh
 * 
 * This example demonstrates a number of methods that can be used to
 * use the libMesh::EquationSystems to create data at the nodal level.
 * 
 * The demo function demonstrates a limitation of libmesh to handle
 * both temporal and spatially changing nodal data with multiple
 * variables. It also presents a solution using the MyAnalyticFunction 
 * class that is based on boost::function. A wrapper class is also 
 * presented that is under development.
 */

//! @cond example3 Doxygen ignores the following code

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

// A class for demonstrating the use and inheritance of the EqDataBase class
class HeatEqData : public EqDataBase{
	public:
		HeatEqData(EquationSystems& sys) : EqDataBase(sys, "heat_eq_data"){
			add_variable("k");		// add conductivity
			add_variable("cp");		// add specific heat
			get_system().init();	// initilize the data storage system
		}
		
		Number k(const Point& p, const Real t){
			return 1 + exp(-t) * sin(pi()*p(0)) * sin(pi()*p(1));
		}

		Number cp(const Point& p, const Real t){
			return 10 + exp(-t) * sin(pi()*p(0)) * sin(pi()*p(1));

		}
				
	protected:		
		double pi(){
			return boost::math::constants::pi<double>();
		}

		void value(DenseVector<Number>& output, const Point& p, const Real){
			
			// Gather the time from the system
			Real t = get_system().time;
			
			// Return the values of k and cp
			output.resize(2);
			output(0) = k(p,t);
			output(1) = cp(p,t);
		}		
};


// A function for space dependent conductivity; time is ignored by libmesh
void func1(DenseVector<Number>& output, const Point& p, const Real t){
	
	// Display the time, it does not get updated (only show at 0,0)
	if(p(0) == 0 && p(1) == 0){
		printf("\nTime = %f\n", t);
	}
	
	// Define pi
	const double pi = boost::math::constants::pi<double>(); 
	
	// Update the output vector
	output.resize(2); 
	output(0) = 1 + exp(-t) * sin(pi*p(0)) * sin(pi*p(1));
	output(1) = 10 + exp(-t) * sin(pi*p(0)) * sin(pi*p(1));
}

// A function for initializing thermal conductivity (boost::function method)
#include <boost/math/constants/constants.hpp>
void func2(DenseVector<Number>& output, const Point& p, const Real, const Parameters& parameters){
		
	// Extract the time	
	Real t = parameters.get<Real>("time");

	// Display the time, only show at 0,0
	if(p(0) == 0 && p(1) == 0){
		printf("\nTime = %f\n", t);
	}
	
	// Define pi
	const double pi = boost::math::constants::pi<double>(); 

	// Update the output vector
	output.resize(2); 
	output(0) = 1 + exp(-t) * sin(pi*p(0)) * sin(pi*p(1));
	output(1) = 10 + exp(-t) * sin(pi*p(0)) * sin(pi*p(1));
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

// METHOD 1: 
// The libmesh way, which seems to not allow multiple variables to 
// update with time. Example 9 does what is required with a wrapper 
// function, but a wrapper function with a vector output is not accepted
// in the project_solution member.
//
// Thus, this first example allows to project multiple variables but 
// can not use time as an input.

	// Create a system and add two variables
	eq_sys.add_system<TransientExplicitSystem>("data");
	eq_sys.get_system("data").add_variable("k", SECOND, MONOMIAL);
	eq_sys.get_system("data").add_variable("cp", SECOND, MONOMIAL);
	
	// Update the system time, this is not used and included here 
	// for demonstration purposes
	eq_sys.get_system("data").time = 1; 

	// Initilize the equation system storing the data
	eq_sys.get_system("data").init();
	
	// Create a AnalyticFuction and use this to project the solution
	AnalyticFunction<Number> fobj1(func1);
	eq_sys.get_system("data").project_solution(&fobj1);
	
	// Export the solution
	ExodusII_IO(mesh).write_equation_systems("example3_1.ex2", eq_sys);

// METHOD 2: 
// This method utilizes MyAnalyticFunction to bind the EquationSystem
// parameters to the function pointer. This then allows for additional
// variables to be passed to the function in similar fashion to the 
// wrapper function in Example 9. Of course, anything could be bound
// using the boost::bind. Thus, this method is infinitaly expandable and
// support multiple variables.

	// Clear the system
	eq_sys.clear();

	// Add two variables
	eq_sys.add_system<TransientExplicitSystem>("data");
	eq_sys.get_system("data").add_variable("k", SECOND, MONOMIAL);
	eq_sys.get_system("data").add_variable("cp", SECOND, MONOMIAL);
	
	// Set the time and intilize the system
	eq_sys.parameters.set<Real>("time") = 1;
	eq_sys.init();

	// Generate the boost::function pointers
	boost::function< void (DenseVector<Number>&, const Point&, const Real, const Parameters&)> fptr_full;
	boost::function< void (DenseVector<Number>&, const Point&, const Real)> fptr_shrt;
	
	// Link function to the full pointer
	fptr_full = func2;

	// Bind the parameters
	fptr_shrt = boost::bind(fptr_full, _1, _2, _3, eq_sys.parameters);

	// Create a MyAnalyticFuction and use this to project the solution
	MyAnalyticFunction<Number> fobj2(fptr_shrt);
	eq_sys.get_system("data").project_solution(&fobj2);
	ExodusII_IO(mesh).write_equation_systems("example3_2.ex2", eq_sys);
	
// METHOD 3:
// A wrapper class that is still under development	
	
	// Clear the system
	eq_sys.clear();
	
	// Create the data class
	HeatEqData data(eq_sys);
	
	// Project the solution at time t
	data.update_solution(1);
	
	// Output the solution to file
	ExodusII_IO(mesh).write_equation_systems("example3_3.ex2", eq_sys);

	return 0;
}
//! @endcond
