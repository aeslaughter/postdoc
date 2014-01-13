 /*! \file front_velocity_system.h
  * \brief Header file for FrontVelocitySystem class.
  */
  
// Limit multiple includes
#ifndef front_velocity_system_h
#define front_velocity_system_h

// Standard library includes
#include <vector>
#include <string>
using std::vector;
using std::string;

// BOOST includes
#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp>

// libMesh includes
#include <libmesh.h>
#include <libmesh_common.h>
#include <point.h>
#include <transient_system.h>
#include <system.h>
#include <equation_systems.h>
#include <explicit_system.h>
using namespace libMesh;

// Includes from this library
#include "fem/common/my_analytic_function.h"
#include "fem/common/system_base.h"

// Add to FEM and HeatEq namespaces
namespace SlaughterFEM{

/*! \class FrontVelocitySystem front_velocity_system.h "fem/volume_average/front_velocity_system.h"
 * \brief A class for explicitly defining the front velocity equation for
 * use in the level set solution.
 * 
 * \ingroup FEMvolavg
 */  	
class FrontVelocitySystem : public SystemBase<TransientExplicitSystem>{

	public:	
		//! Class constructor
		/*!
		 */
		FrontVelocitySystem(EquationSystems& es, const string& name, const unsigned int number);

		void add_velocity_function(boost::function<void (DenseVector<Number>&, const Point&p, const Real t)> func);
		
		virtual void update_solution(Real t = 0);
		
		//DenseVector<Number> velocity(const Point& p, const Real t);
		
	protected:
		boost::function<void (DenseVector<Number>&, const Point&p, const Real t)> _velocity_ptr;
		
		virtual void initialize();
		
		void project(Real t = 0);

}; // class
}  // namespace
#endif
