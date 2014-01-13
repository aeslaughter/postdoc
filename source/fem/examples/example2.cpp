/** \example fem/examples/example2.cpp 
 * A test function for testing the HeatEq class.
 * 
 * This class implements Example 8-1 from Bhatti (2005; p. 552). Run
 * the program with the \c --help flag for a complete list of run-time
 * options.
 * 
 * This class will fail to compile in debug mode. libMesh and 
 * Boost::program_options do not get along.
 */
 
// Standard C++ includes  
#include <iostream>
#include <algorithm>
#include <math.h>
#include <stdio.h>
#include <functional>

// BOOST includes
#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>
#include <boost/ref.hpp>
#include <boost/bind.hpp>
#include <boost/math/constants/constants.hpp>

// libMesh includes
#include <libmesh.h>
#include <libmesh_common.h>
#include <exodusII_io.h>
#include <vtk_io.h>
#include <gmv_io.h>
#include <mesh_generation.h>
#include <mesh.h>
#include <mesh_triangle_interface.h>
#include <mesh_generation.h>
#include <elem.h>
#include <mesh_tetgen_interface.h>
#include <node.h>
#include <face_tri3.h>
#include <mesh_triangle_holes.h>
using namespace libMesh;

// Include portions of my fem library
#include "fem/common/include.h"
#include "fem/heat_eq/include.h"
using namespace SlaughterFEM;

// Include my common library
#include "common/include.h"
using namespace SlaughterCommon;

// Short-hand for boundary conditions
typedef boost::shared_ptr<HeatEqBoundaryDirichlet> pDirichlet;
typedef boost::shared_ptr<HeatEqBoundaryNeumann> pNeumann;
typedef boost::shared_ptr<HeatEqBoundaryConvection> pConvection;

// My boundary condition function
void dirichlet_function(DenseVector<Number>& output, const Point&, const Real){
	output(0) = 300;
}	

// Wrapper function for exact solution, used for initilization        
Number initial_function(const Point&, const Parameters&, const std::string&, const std::string&){
    return 50;
}

// Function prototype for gathering command-line options with UserOptions
UserOptions gather_command_line(int argc, char** argv);

// Begin main function
int main (int argc, char** argv){
//PerfLog P("Example 2 Program");
//P.push("init","Program Initilization");

	// Gather command-line options
	UserOptions user = gather_command_line(argc, argv);

    // Initialize libraries
    LibMeshInit init (argc, argv);
        
	// Initilize the mesh and order variables
	MyMesh mesh;
	Order order;
	
	// Patch test (Bhatti Example 8-1)
	if(user.get_flag("patch")){
		GMVIO(mesh).read("../data/fem/examples/input/example2.gmv");
		mesh.all_first_order();
		order = FIRST;
	
	// 2-D box
	} else if (user.get_flag("2D")){
		MeshTools::Generation::build_square(mesh, 10, 10, 0., 0.04, 0., 0.04, TRI6);
		mesh.all_second_order();
		order = SECOND;
		
	// 3-D box	
	} else if (user.get_flag("3D")){
		MeshTools::Generation::build_cube(mesh, 10, 10, 10, 0., 0.04, 0., 0.04, 0., 0.04, TET10);
		mesh.all_second_order();
		order = SECOND;
	
	// Multi-element implementation of Bhatti Exmample 8-1
	} else {
		mesh.set_mesh_dimension(2);
		mesh.add_point(Point(0,0));
		mesh.add_point(Point(0.02,0));
		mesh.add_point(Point(0.02,0.04));
		mesh.add_point(Point(0,0.02));
	
		TriangleInterface t(mesh);
		t.desired_area() = 1e-4;
		t.triangulation_type() = TriangleInterface::PSLG;
		t.smooth_after_generating() = true;
		t.triangulate();
		mesh.all_second_order();
		order = SECOND;
	}
	
	// Create an equation system
 	EquationSystems eq_sys(mesh); 
 
	// Create a HeatEq class
	boost::shared_ptr<HeatEq> heateq(new HeatEq(eq_sys, order));
	 	
 	// Define the material constants
	heateq->system().set_constant<Real>("k", user.get<Real>("conductivity"));
    heateq->system().set_constant<Real>("rho", user.get<Real>("density"));
    heateq->system().set_constant<Real>("cp", user.get<Real>("specific-heat"));

	// Link to the initialization function
	heateq->system().add_initial_function(initial_function);

    // Add boundary IDs (this is some custom functionality that I added)
    mesh.find_neighbors();
    mesh.boundary_info->clear();
    mesh.add_boundary_id(0, "y", 0.0);  // bottom
    mesh.add_boundary_id(1, "x", 0.02); // right
    mesh.add_boundary_id(2, "x", 0.0);  // left
    mesh.add_boundary_id(3); // top

	// Convection boundary at bottom (user-specified)
	pConvection pC = heateq->system().add_boundary<HeatEqBoundaryConvection>(0);
	pC->h_constant = user.get<Real>("h-coefficient");
	pC->Tinf_constant = user.get<Real>("Tinf");	

	// Flux boundary at right-side (user-specified)
	pNeumann pN = heateq->system().add_boundary<HeatEqBoundaryNeumann>(1);
	pN->q_constant = user.get<Real>("flux");
	
	// Flux boundary at left-side (symetry; defaults to q = 0)
	heateq->system().add_boundary<HeatEqBoundaryNeumann>(2);
	
	// Top constant temperature boundary
	pDirichlet pD = heateq->system().add_boundary<HeatEqBoundaryDirichlet>(3);
	pD->fptr = dirichlet_function; // links the boundary function

	// Initialize system
	heateq->system().init(0.0);

	// Define a general filename
    FileParts outfile("../data/fem/examples/output/example2.ex2");
       
    // Export the initial mesh
	//ExodusII_IO(mesh).write_equation_systems(outfile.add_tstep(0,3,"_"), eq_sys);
	MyVTKIO vtk("../data/fem/examples/output/example2.vtu", eq_sys);
	vtk.write(0.0);

	// Define time stepping variables
	Real time  = 0.;
	Real dt = user.get<Real>("dt");	
	
//P.pop("init");
//P.push("soln", "Transient Solution");  
    
    // Begin the time loop
    int N = user.get<int>("num-steps");
    int d = user.get<int>("output-div");
	for (int t_step = 0; t_step < N; t_step++){
			
		// Incremenet the time counter, set the time and the
		// time step size as parameters in the EquationSystem.
		time += dt;

		// Display a progress message
		printf("time = %f; step %d of %d\n", time, t_step, N);

		// Update the old solution vector
		heateq->system().update_solution(time, dt);
		
		//! \todo This needs to be incoporated into the HeatEqSystem
		heateq->system().reinit();
		heateq->system().rhs->zero();

		// Solve the system
		heateq->system().solve();
		
		// Output evey 10 timesteps to file.
		if ( (t_step+1)%d == 0){
			std::string out = outfile.add_tstep(t_step+1,3,"_");
			//ExodusII_IO(mesh).write_equation_systems(out, eq_sys);
			vtk.write(time);
		}
    }

//P.pop("soln");
	// All done. 
	return 0;
}

// A subfunction for defining and gathering command-line options
UserOptions gather_command_line(int argc, char** argv){
	
	// The general options and title
	UserOptions user("General Options");
	user.add_title("Example 2: FEM solution of the heat equation\n");
	user.add_flag("help","List the available options");
	
	// Domain related options
	UserOptions type("Domain Options");
	type.add_flag("patch", "Run as 2-element patch test as in Bhatti example 8-1");
	type.add_flag("2D", "Run with a 2-D box domain");
	type.add_flag("3D", "Run with a 3-D cube domain"); 
	
	// Time integration related options
	UserOptions t_opt("Time Intergration Options");
	t_opt.add_option<int>("num-steps,n", 300, "Number of time steps");
	t_opt.add_option<Real>("dt",1,"Time step division (sec.)");
	t_opt.add_option<int>("output-div,d", 10, "About the data after this many time steps");

	// Mesh refinement options
	UserOptions r_opt("Mesh Refinement Options");
	r_opt.add_flag("refine", "Utilize adaptive mesh refinement");
	r_opt.add_option<double>("refine-fraction,r", 0.80, "Max. fraction of elements to refine");
	r_opt.add_option<double>("coarsen-fraction,c", 0.07, "Max. fraction of elements to coarsen");
	r_opt.add_option<int>("h-level,l", 5, "Max. allowed refinement steps for element");

	// Material constant options
	UserOptions m_opt("Material Constants (defaults to Bhatti Example 8-1)");
	m_opt.add_option<Real>("conductivity,k",3, "Thermal conductivity (W/(mK))");
	m_opt.add_option<Real>("density,p", 1600, "Density (kg/m^3)");
	m_opt.add_option<Real>("specific-heat,c", 800, "Specific heat capacity (J/(kgK))");
	
	// Boundary condition options
	UserOptions b_opt("Boundary Options (defaults to Bhatti Example 8-1)");
	b_opt.add_option<Real>("flux,q", 0, "Flux boundary value, right side (W/m^2)");
	b_opt.add_option<Real>("h-coefficient,h", 200, "Convection heat transfer coefficient (W/m^2)");
	b_opt.add_option<Real>("Tinf,i", 50, "Convection boundary layer temperature (C)");
	b_opt.add_option<Real>("T-top,T", 300, "Top boundary temperature (disabled)");

	// Link the groups together
	user.add(type).add(t_opt).add(r_opt).add(m_opt).add(b_opt);

	// Apply the command-line options
	user.apply_options(argc, argv);
	
	// Return the main class
	return user;
}


