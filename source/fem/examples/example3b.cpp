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

#include "fem/common/eq_core.h"
using namespace SlaughterFEM;

// Bring in everything from the libMesh namespace
using namespace libMesh;

// A class for initilizing velocity
template <typename Output = Number>
class EqInit : public FunctionBase<Output>{
	
	public:
		
		EquationSystems& eq_sys;
		
		EqInit(EquationSystems& sys) : eq_sys(sys){
			eq_sys.add_system<TransientExplicitSystem>("data");
			eq_sys.get_system("data").add_variable("x");
			eq_sys.get_system("data").add_variable("y");
		};
		
		// Returns the velocity by its components
		virtual Output component(unsigned int index, const Point& p, Real t = 0){
		
			PointLocatorTree locate(eq_sys.get_mesh());

			if (index == 0){
				return x(p); // add some intensive calculation here
			} else {
				return y(p); // add some intensive calculation here
			}		
		} 
		
		// Do not use the scalar output version of the operator()
		Output operator() (const Point& p, const Real t = 0.){
			libmesh_not_implemented();
		}

		// Use the component member to build the vector output version of operator()
		void operator() (const Point& p, const Real t, DenseVector<Output>& output){
				
			output.resize(eq_sys.get_system("data").n_vars());
			
			for (unsigned int i = 0; i < output.size(); i++){
				output(i) = component(i,p,t);
			}
		}
		
		// Create a copy of the class
		virtual AutoPtr<FunctionBase<Output> > clone () const{
			return AutoPtr<FunctionBase<Output> > (new EqInit<Output>(eq_sys));
		}	
		
		Output x(const Point&){
			return 1;
		}
		
		Output y(const Point&){
			return 2;
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

/*
	// Create a system and add two variables (x and y velocity components)
	eq_sys.add_system<TransientNonlinearImplicitSystem>("test");
	eq_sys.get_system("test").add_variable("vx", FIRST);
	eq_sys.get_system("test").add_variable("vy", FIRST);

	// Initilize the equation system
	eq_sys.get_system("test").init();
*/	
	
	// Project the velocity using the EqInit class
	EqInit<Number> data(eq_sys);
	eq_sys.get_system("data").init();
	
	eq_sys.get_system("data").project_solution(&data);
	
	Point p(1,1);
	
	Number x = data.point_value(p);
	Number y = data.point_value(p);
	
	printf("x = %g; y = %g\n", x, y);
/*
	eq_sys.get_system("test").project_solution(&init_velocity);
	
	// Test that values are correct for each point on the only element
	DenseVector<Number> v(2);
	Elem* elem = mesh.elem(0);
	for (unsigned int i = 0; i < elem->n_nodes(); i++){
		Point& p = elem->point(i);
		v(0) = eq_sys.get_system("test").point_value(0,p);
		v(1) = eq_sys.get_system("test").point_value(1,p);
		printf("Point %g,%g; v = %f,%f\n", p(0),p(1),v(0),v(1));
	}
*/	

	// Done
	return 0;
}
