 /*! \file front_velocity_eq.h
  * \brief Header file for FrontVelocityEqclass.
  */
  
// Limit multiple includes
#ifndef front_velocity_eq_h
#define front_velocity_eq_h

// Standard library includes
#include <vector>
#include <string>
using std::vector;
using std::string;

// libMesh includes
#include <libmesh.h>
#include <libmesh_common.h>
using namespace libMesh;

// Includes from this library
#include "fem/common/equation_base.h"
#include "fem/volume_average/front_velocity_system.h"

// Add to FEM and HeatEq namespaces
namespace SlaughterFEM{

/*! \class LevelSetSystem level_set_system.h "fem/volume_average/level_set_system.h"
 * \brief A class for solving the level set equation

 * 
 * \ingroup FEMvolavg
 */  	
class FrontVelocityEq : public EquationBase<FrontVelocitySystem>{

	public:	
		//! Class constructor
		/*!
		 */
		FrontVelocityEq(EquationSystems& es, const Order order = FIRST, const FEFamily family = MONOMIAL);


}; // class
}  // namespace
#endif
