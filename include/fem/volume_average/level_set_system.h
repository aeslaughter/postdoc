 /*! \file level_set_system.h
  * \brief Header file for LevelSetSystem class.
  */
  
// Limit multiple includes
#ifndef level_set_system_h
#define level_set_system_h

// Standard library includes
#include <vector>
#include <string>
using std::vector;
using std::string;

// BOOST includes
#include <boost/shared_ptr.hpp>
#include <boost/math/constants/constants.hpp>

// LibMesh: general includes
#include <libmesh.h>
#include <libmesh_common.h>
#include <point.h>
#include <string_to_enum.h>
#include <getpot.h>
#include <mesh.h>
#include <boundary_info.h>
#include <dirichlet_boundaries.h>
#include <analytic_function.h>
#include <mesh_refinement.h>
#include <parmetis_partitioner.h>
#include <error_vector.h>
#include <kelly_error_estimator.h>
#include <linear_implicit_system.h>
#include <equation_systems.h>
#include <transient_system.h>
#include <vector_value.h>
#include <fe.h>
#include <fe_interface.h>
#include <fe_type.h>
#include <dof_map.h>
#include <quadrature_gauss.h>
#include <sparse_matrix.h>
#include <numeric_vector.h>
#include <dense_matrix.h>
#include <dense_vector.h>
#include <elem.h>
#include <perf_log.h>
#include <quadrature_simpson.h>
using namespace libMesh;

// Includes from this library
#include "fem/common/my_analytic_function.h"
#include "fem/common/implicit_system_base.h"
#include "fem/common/boundary_base.h"
#include "fem/volume_average/front_velocity_eq.h"

// Add to FEM and HeatEq namespaces
namespace SlaughterFEM{

/*! \class LevelSetSystem level_set_system.h "fem/volume_average/level_set_system.h"
 * \brief A class for solving the level set equation

 * 
 * \ingroup FEMvolavg
 */  	
class LevelSetSystem : public ImplicitSystemBase<TransientLinearImplicitSystem>{

	public:	
		//! Class constructor
		/*!
		 */
		LevelSetSystem(EquationSystems& es, const string& name, const unsigned int number);
		
		~LevelSetSystem();
		
		boost::shared_ptr<FrontVelocityEq> velocity;
		
		//! DG-SSPRK Solution
		/*!
		 */
		void solve();
		
		//! Pointer to the front velocity equation
		//boost::shared_ptr<FrontVelocityEq> 
		
		void update_solution();
		
		void update_solution(Real time, Real dt);


		
	//	SparseMatrix<Number>* mass_matrix_inverse;
		
		
	//protected:	
		
		//!
		void initialize();
	
		//! Time step calculatin
		Number time_step();
	
		//! libMesh assembly function
		void assemble();
		
		void assemble_mass();
		
		void assemble_stiffness();
		
		// make this a private wrapper to a function pointer supplied by user
		//VectorValue<Number> velocity(const Point& x, const Real t);
		
		
		Number h;
		
		unsigned int _count;
		

}; // class

} // namespaces
#endif
