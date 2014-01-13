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
#include "system.h"
#include "mesh_generation.h"
#include "equation_systems.h"
#include "nonlinear_implicit_system.h"
#include "transient_system.h"
#include "explicit_system.h"
#include "analytic_function.h"
#include "function_base.h"
#include "exodusII_io.h"
#include "point.h"
#include "point_locator_tree.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

// A class for initilizing velocity
template <typename Output = Number>
class EqInit : public FunctionBase<Output>{
	
	public:
		
		// Storage location for the EquationSystems
		EquationSystems& eq_sys;
		
		// Constructors that creates some variables
		EqInit(EquationSystems& sys) : eq_sys(sys){
			eq_sys.add_system<TransientExplicitSystem>("data");
			eq_sys.get_system("data").add_variable("x");
			eq_sys.get_system("data").add_variable("y");
		};
		
		// Returns the velocity by its components
		virtual Output component(unsigned int index, const Point& p, Real t = 0){
		
			printf("Index: %i\n", index);
			printf("\tp(0) = %g; p(1) = %g\n", p(0), p(1));
		
			// Locate the element for the current point
			PointLocatorTree locate(eq_sys.get_mesh());
			const Elem* elem = locate(p);

			if (index == 0){
				printf("\tx = %g\n", x(p));
				return x(p); // add some intensive calculation here
			} else {
				printf("\ty = %g\n", y(elem));
				return y(elem); // add some intensive calculation here
			}		
			
			
		} 
		
		// Do not use the scalar output version of the operator()
		Output operator() (const Point& p, const Real t = 0.){
			libmesh_not_implemented();
		}

		// Use the component member to build the vector output version of operator()
		void operator() (const Point& p, const Real t, DenseVector<Output>& output){
			
			// Re-size the output to the number of variables
			output.resize(eq_sys.get_system("data").n_vars());
			
			// Use the component function to compute the values
			for (unsigned int i = 0; i < output.size(); i++){
				output(i) = component(i,p,t);
			}
		}
		
		// Create a copy of the class
		virtual AutoPtr<FunctionBase<Output> > clone () const{
			return AutoPtr<FunctionBase<Output> > (new EqInit<Output>(eq_sys));
		}	
		
		// A parameter that has some complicated calculation
		Output x(const Point& p){ 
			return (Output) p(0); 
		}
		
		// Another parameter that takes some work
		Output y(const Elem* elem){ 
			return (Output) elem->volume(); 
			printf("\ty = %g\n",elem->volume());
		}
};

// The main function
int main (int argc, char** argv){

	// Initialize libraries
    LibMeshInit init (argc, argv);
        
	// Generate a mesh
	Mesh mesh;
	MeshTools::Generation::build_square(mesh, 1, 1, -1., 1, -1, 1, QUAD4);
	mesh.all_first_order();
	
	// Create an equation system
	EquationSystems eq_sys(mesh);

	// Project the data using the EqInit class
	EqInit<Number> data(eq_sys);
	eq_sys.get_system("data").init();
	eq_sys.get_system("data").project_solution(&data);
	
	// Test that the projection is working
	//Point p(1,1);
	//Number x = eq_sys.get_system("data").point_value(0,p);
	//Number y = eq_sys.get_system("data").point_value(1,p);
	

	// Done
	return 0;
}
