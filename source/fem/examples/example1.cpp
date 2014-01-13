/** \example fem/examples/example1.cpp 
 * A test function for testing the HeatEq class.
 * 
 * A 2-D FEM solution to the following equation on the domain from 0
 * to 1 in both the x and y directions.
 * 
 * \f[
 * 1 + \exp(-t) * \sin(\pi x) \sin(\pi y)
 * \f]
 */
 
// Standard includes  
#include <stdio.h>
#include <math.h>

// BOOST includes
#include <boost/math/constants/constants.hpp>

// Include my FEMcommon and FEMheateq libraries
#include "fem/common/include.h"
#include "fem/heat_eq/include.h"
using namespace SlaughterFEM;

// Include my common library
#include "common/include.h"
using namespace SlaughterCommon;

// libMesh related includes
#include <mesh_generation.h>
#include <vtk_io.h>
#include <exodusII_io.h>
using namespace libMesh;

// This is the exact solution, it is also used for defining boundary
// contidions and the initial condition 
Number exact_solution (const Point& p, const Real t){
	const double pi = boost::math::constants::pi<double>();  
	return 1 + exp(-t) * sin(pi*p(0)) * sin(pi*p(1));
}

void exact_solution_vec(DenseVector<Number>& output, const Point& p, const Real t){
	const double pi = boost::math::constants::pi<double>();  
	output(0) = 1 + exp(-t) * sin(pi*p(0)) * sin(pi*p(1));
}

// Wrapper function for exact solution, used for initilization        
Number initial_function(const Point& p, const Parameters&, const std::string&, const std::string&){
    return exact_solution(p, 0.0);
}

// Wrapper function for boundary function
void boundary_function(DenseVector<Number>& output, const Point& p, const Real t){
	output(0) = exact_solution(p, t);
}

// Begin main function
int main (int argc, char** argv){
  
 /* These fail in libMesh's debug mode
 
	// Add command-line options
	UserOptions opt("Program options");
	opt.add_flag("help","List the available options");
	opt.add_option<int>("nx",10, "Number of elements in x direction");
	opt.add_option<int>("ny",10, "Number of elements in y direction");
	opt.add_option<int>("num-steps,N", 100, "Number of time steps");
	opt.add_option<double>("t-step,t", 0.01, "Time step (sec.)");
	opt.add_option<double>("t-start", 0, "Start time (sec.)");
	opt.apply_options(argc, argv);

	// Collect a few options
	int nx = opt.get<int>("nx");
	int ny = opt.get<int>("ny");
	int N  = opt.get<int>("num-steps");
	double dt = opt.get<double>("t-step");
	double tstart = opt.get<double>("t-start");
 */
 	// Debug version (libmesh debug build doesn't like boost::program_options)
	int nx = 10;
	int ny = 10;
	int N  = 100;
	double dt = 0.01;
	double tstart = 0;

    // Initialize libraries
    LibMeshInit init (argc, argv);
        
    // Create a mesh 
    MyMesh mesh;
	MeshTools::Generation::build_square(mesh, nx, ny, 0., 1., 0., 1., QUAD8);
	mesh.all_second_order();

    // Create a HeatEq object
    EquationSystems eq_sys(mesh); 
    boost::shared_ptr<HeatEq> heateq(new HeatEq(eq_sys));

	// Assign an initilization function
	heateq->system().add_initial_function(initial_function);

	// Define a parameters for this problem
	const double pi = boost::math::constants::pi<double>();  
    heateq->system().set_constant<Real>("k", 1 / (2 * pi * pi));
   
    // Add boundary IDs 
    mesh.find_neighbors();
    mesh.boundary_info->clear();
    mesh.add_boundary_id(0); // all boundaries

    // Dirichlet Boundary uusing the standard function pointer
	boost::shared_ptr<HeatEqBoundaryDirichlet> pBC0_func = 
		heateq->system().add_boundary<HeatEqBoundaryDirichlet>(0);
	pBC0_func->fptr = boundary_function;

	// Dirichlet Boundary using boost::function
	//boost::shared_ptr<HeatEqBoundaryDirichlet> pBC0 = eqObj.add_boundary<HeatEqBoundaryDirichlet>(1);
	//pBC0->T_constant = 1;
	
	// Method of initializing with a boost::function (overrides above)
	boost::function<void (DenseVector<Number>& output, const Point&, const Real)> init_ptr;
	init_ptr = exact_solution_vec;
	heateq->system().add_initial_function(init_ptr);
	
	// Initialize the equation system
	double time = tstart; 
	heateq->system().init(time);
	
	// Define a general filename
    FileParts outfile("../data/fem/examples/output/example1.ex2");
  
    // Export the initial mesh
	ExodusII_IO(mesh).write_equation_systems(outfile.add_tstep(0,3,"_"), eq_sys);
	
	// Loop through time
	for (int t_step = 0; t_step < N; t_step++){
		
		// Incremenet the time counter, set the time and the
		// time step size as parameters in the EquationSystem.
		time += dt;
	
		// Display a progress message
		printf("time = %f; step %d of %d\n", time, t_step, N);

		// Update the old solution vector
		heateq->system().update_solution(time, dt);

		// Assemble and solve the linear system
		heateq->system().solve();
 
		// Output evey 10 timesteps to file.
		if ( (t_step+1)%10 == 0){
			std::string out = outfile.add_tstep(t_step+1,3,"_");
			ExodusII_IO(mesh).write_equation_systems(out, eq_sys);
        }
    }

	// All done. 
	return 0;
}




