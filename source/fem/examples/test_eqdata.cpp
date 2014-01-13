/** \example test_eqdata.cpp 
 * A test function for the VolumeAverageEqData class
 */
 
 #include <iostream>
 
 
// LibMesh related includes
#include "libmesh.h"
#include "mesh.h"
#include "mesh_generation.h"
using namespace libMesh;

// BOOST includes
#include <boost/shared_ptr.hpp>

// My custom fem library includes
#include "fem/include.h"
using namespace SlaughterFEM;


//typedef EqVariableLinker<TransientNonlinearImplicitSystem> Linker;


// Test element

// NODE ID:
// 		3     2
// 		*-----*
//		|	  |		
// 		*-----*
// 		0     1
//
// NODE DATA (id; x,y position; x,y velocity)
// 	0; -1,-1; 1,2
//	1; 1,-1;  1,0
//	2; 1,1;   1,1
//	3; -1,1;  0,0

void initial_velocity(DenseVector<Number>& output, const Point& p, const Real){

 	output.resize(2);
	if (p(0) == -1 && p(1) == -1){
		output(0) = 1; // x-direction
		output(1) = 2; // y-direction
	} else if (p(0) == 1 && p(1) == -1){
		output(0) = 1;
		output(1) = 0;		
	} else if (p(0) == 1 && p(1) == 1){
		output(0) = 1;
		output(1) = 1;	
	} else if (p(0) == -1 && p(1) == 1){
		output(0) = 0;
		output(1) = 0;	
	}

}

void initial_concentration(DenseVector<Number>& output, const Point& p, const Real){
	output(0) = 0.192;
}

void initial_enthalpy(DenseVector<Number>& output, const Point& p, const Real){
	output(0) = 1;
}

void initial_temperature(DenseVector<Number>& output, const Point& p, const Real){
	output(0) = 287;
}

int main(int argc, char** argv){   
	
	// Initilize libMesh
    LibMeshInit init (argc, argv);

	// Create the mesh, single element
	Mesh mesh;
    MeshTools::Generation::build_square(mesh, 1, 1, -1, 1, -1, 1, QUAD4); 

	// Create an equation system
	EquationSystems eq_sys(mesh);

	// Create the equation systems via shared_ptr pointers
	boost::shared_ptr<ThermoEq> thermo(new ThermoEq(eq_sys));
	boost::shared_ptr<MomentumEq> momentum(new MomentumEq(eq_sys, FIRST)); 
	boost::shared_ptr<ConcentrationEq> concentration(new ConcentrationEq(eq_sys, FIRST));
	boost::shared_ptr<EnergyEq> energy(new EnergyEq(eq_sys, FIRST));

	// Link the initialization functions
	momentum->add_initial_function(initial_velocity);
	concentration->add_initial_function(initial_concentration);
	energy->add_initial_function(initial_temperature);

	// Add the necessary cross-linking between the systems via their pointers
	thermo->momentum = momentum;
	thermo->energy = energy;
	thermo->concentration = concentration;
	energy->thermo = thermo;
	energy->momentum = momentum;

	// Initialize the various systems
	momentum->init();
	concentration->init();
	energy->init();
	thermo->init();

	// Get a pointer to the first (and only) element
	Elem* elem = mesh.elem(0);
	
	// Display the values for initial enthalpy; computed from temp.
	printf("\n\nTEMP. TO ENTHALPY CONVERSION:\n");
	for (unsigned int i = 0; i < elem->n_nodes(); i++){
		Point& p = elem->point(i);
		Number h = energy->point_value(0,p);
		printf("\th = %g (%g, %g)\n", h, p(0), p(1));
	}		
	
	
	
	
//	thermo->test(elem);
//	DenseVector<Number> output(1);
	
	//energy.enthalpy(output, p, 0.0);
	
	//energy.init();
	
	// Thermodynamics links

	//data.init();

	//printf("initilized = %s\n", (data.velocity_ptr->initialized())?"true":"false");

/*
	// Print some basic element information
	printf("\n\nELEMENT INFORMATION:\n");	
		
	// Test that values are correct, use first (only) element
	VectorValue<Number> v;
	Elem* elem = mesh.elem(0);
	for (unsigned int i = 0; i < elem->n_nodes(); i++){
		Point& p = elem->point(i);
		v = data.velocity_ptr->point_value(p);
		printf("\tPoint %g, %g\tv = %g, %g\n", p(0),p(1),v(0),v(1));
	}	
 	// Test the element length
	printf("\nTESTING: element data\n");
	Number h = data.element_length(elem);
	printf("\tElement length (%f): %f\n", 1.065004, h);

	// Test the tau_1 value
	printf("\nTESTING: Tau_1\n"); 
	printf("\tTau_1: %g\n", data.tau_1(elem->point(0), data.element_length(elem)));
	printf("\tTau_1: %g\n", data.tau_1(elem->point(1), data.element_length(elem)));
	printf("\tTau_1: %g\n", data.tau_1(elem->point(2), data.element_length(elem)));
	printf("\tTau_1: %g\n", data.tau_1(elem->point(3), data.element_length(elem)));
 
 */
 
 
 
 
	// Test the EnergyEq assembly
 	printf("\nTESTING: energy eq. assembly\n");
	energy->assemble();
}
