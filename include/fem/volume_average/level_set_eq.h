/*! \file level_set_eq.h
  * \brief Header file for LevelSetEq class.
  */

// Avoid multiple includes
#ifndef level_set_eq_h
#define level_set_eq_h

// Standard C++ includes
#include <string>
using std::string;

// libMesh includes 
#include <libmesh.h>
#include <libmesh_common.h>
#include <system.h>

// My includes  
#include "fem/volume_average/level_set_system.h"
#include "fem/common/equation_base.h"
  
namespace SlaughterFEM{
	
/*! \class LevelSetEq level_set_eq.h "fem/volume_average/level_set_eq.h"
 * \brief An equation class for the level set equation.
 * 
 * \ingroup FEMvolavg
 */  
class LevelSetEq : public EquationBase<LevelSetSystem>{
	
	public:
		LevelSetEq(EquationSystems& es, const Order order = FIRST, const FEFamily family = LAGRANGE); 
		
}; // class
}  // namespace
#endif
