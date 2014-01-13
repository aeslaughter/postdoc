/** \example fem/examples/example6.cpp 
 */
 
// Standard includes  
#include <stdio.h>
#include <math.h>
#include <string>
using std::string;

// Include my common library
#include "common/include.h"
using namespace SlaughterCommon;

// Include my FEMcommon and FEMvolavg libraries
#include "fem/common/include.h"
#include "fem/volume_average/include.h"
using namespace SlaughterFEM;

// libMesh related includes
#include <libmesh.h>
#include <mesh.h>
#include <mesh_generation.h>
#include <vtk_io.h>
#include <exodusII_io.h>
using namespace libMesh;	
	

// Initial function
void phi_init(DenseVector<Number>& output, const Point& x, const Real){	
	const double pi = boost::math::constants::pi<double>();  
	output(0) = std::pow(x(0) - 0.5, 2) + std::pow(x(1) - 0.75, 2) - std::pow(0.15, 2);
}

// This needs to be outside of the System
void velocity_func(DenseVector<Number>& output, const Point& x, const Real t){
	output.resize(2);
	const double pi = boost::math::constants::pi<double>();  
	output(0) = cos(pi*t/8)*sin(2*pi*x(1))*std::pow(sin(pi*x(0)),2);
	output(1) = cos(pi*t/8)*sin(2*pi*x(0))*std::pow(sin(pi*x(1)),2);
}

int main (int argc, char** argv){

    // Initialize libraries
    LibMeshInit init (argc, argv);
  
    // Create a mesh 
    MyMesh mesh;
	MeshTools::Generation::build_square(mesh, 20, 20, 0., 1., 0., 1., QUAD4);
	mesh.all_first_order();
	//MeshTools::Generation::build_square(mesh, 2, 2, 0., 1., 0., 1., QUAD8);
	//mesh.all_second_order();

    // Create a LevelSetSystem
    EquationSystems eq_sys(mesh); 	
	boost::shared_ptr<FrontVelocityEq> velocity(new FrontVelocityEq(eq_sys, SECOND));
	velocity->system().add_velocity_function(velocity_func);
	velocity->system().init();

    LevelSetSystem& levelset = eq_sys.add_system<LevelSetSystem>("LevelSetEquation");	
	levelset.velocity = velocity;
    levelset.add_initial_function(phi_init);
    levelset.init();
    
	FileParts outfile("../data/fem/examples/output/example6.ex2");
	ExodusII_IO(mesh).write_equation_systems(outfile.add_tstep(0,5,"_").c_str(), eq_sys);
	
	Number t_stop = 8;
	Number t = 0;
	unsigned int cnt = 1;

	while (t < t_stop){
		
		// Solve
		levelset.solve();
		
		//velocity.solve();
		
		velocity->system().update_solution();
	
		t = levelset.time;
		velocity->system().time = t;
	
		printf("Time = %g (%d)\n", t, cnt);
		ExodusII_IO(mesh).write_equation_systems(outfile.add_tstep(cnt, 5, "_").c_str(), eq_sys);
		cnt++;
		
		//levelset.update_solution();
	}
	
	//levelset.reinit();
	//printf("M Global:\n");
	//levelset.get_matrix("mass_matrix_inverse").print(std::cout);
	//levelset.matrix->print(std::cout);	
  
}
