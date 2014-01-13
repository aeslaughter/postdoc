/*! \example test_thermo.cpp 
 * A test function for the thermodynamic relations in the ThermoEq class,
 * the results from this function should be compared with the MATLAB 
 * function test_thermo.m.
 */
 
// My includes
//#include "fem/test/test_functions.h"
#include "fem/common/system_base.h"

#include "fem/include.h"
using namespace SlaughterFEM;

// Use the libMesh namespace
using namespace libMesh;
#include "dof_map.h"
#include "exodusII_io.h"





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
	output.resize(1);
	if (p(0) == -1 && p(1) == -1){
		output(0) = 0.192; // x-direction
	} else if (p(0) == 1 && p(1) == -1){
		output(0) = 0.192;
	} else if (p(0) == 1 && p(1) == 1){
		output(0) = 0.192;
	} else if (p(0) == -1 && p(1) == 1){
		output(0) = 0.192;
	}
}

void initial_temperature(DenseVector<Number>& output, const Point& p, const Real){
	output.resize(3);
	if (p(0) == -1 && p(1) == -1){
		output(0) = 287; 
	} else if (p(0) == 1 && p(1) == -1){
		output(0) = 287;
	} else if (p(0) == 1 && p(1) == 1){
		output(0) = 287;
	} else if (p(0) == -1 && p(1) == 1){
		output(0) = 287;
	}
	
	// h_dot and delta_h_dot
	//! \todo this needs to be handled internally
	output(1) = 0;
	output(2) = 0;
	
}

void new_enthalpy(DenseVector<Number>& output, const Point& p, const Real){
	output.resize(1);
	double x = 0.9;
	if (p(0) == -1 && p(1) == -1){
		output(0) = x*3.8002; // x-direction
	} else if (p(0) == 1 && p(1) == -1){
		output(0) = x*2.8499;
	} else if (p(0) == 1 && p(1) == 1){
		output(0) = x*3.5154;
	} else if (p(0) == -1 && p(1) == 1){
		output(0) = x*4.2149;
	} 
}
	
int main(int argc, char** argv){   

	// Initilize libMesh
    LibMeshInit init (argc, argv);

	// Create the mesh, single element
	Mesh mesh;
    MeshTools::Generation::build_square(mesh, 1, 1, -1, 1, -1, 1, QUAD4); 

	// Create an equation system
	EquationSystems eq_sys(mesh);

	//boost::shared_ptr<EnergySystem> energy(new EnergySystem(eq_sys, FIRST));


	// Create the equation systems via shared_ptr pointers
	boost::shared_ptr<ThermoEq> thermo(new ThermoEq(eq_sys));
	//boost::shared_ptr<MomentumEq> momentum(new MomentumEq(eq_sys, FIRST)); 
	//boost::shared_ptr<ConcentrationEq> concentration(new ConcentrationEq(eq_sys, FIRST));
	//boost::shared_ptr<EnergyEq> energy(new EnergyEq(eq_sys, FIRST));
/*
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
	printf("\nTEMP. TO ENTHALPY CONVERSION:\n");
	for (unsigned int i = 0; i < elem->n_nodes(); i++){
		Point& p = elem->point(i);
		Number h = energy->point_value(0,p);
		//printf("\th = %g (%g, %g)\n", h, p(0), p(1));
	}		
	
	// Display the variables values for Eqs. 17 to 21
	//thermo->test(elem);
	
	//unsigned int idx = thermo->get_system().variable_number("density");
	//printf("idx = %d\n", idx);
	
	//const DofMap& dof_map = thermo->get_system.get_dof_map();
	
	
	
	// Test the EnergyEq assembly
 	printf("\nTESTING: energy eq. assembly\n");
	energy->solve();
	*/
	
	
}
