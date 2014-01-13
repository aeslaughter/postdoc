/*! \file heat_eq.h
  * \brief Header file for HeatEq class.
  */

// Avoid multiple includes
#ifndef heat_eq_h
#define heat_eq_h

// Standard C++ includes
#include <string>
using std::string;

// libMesh includes 
#include <libmesh.h>
#include <libmesh_common.h>
#include <system.h>

// My includes  
#include "fem/heat_eq/heat_eq_system.h"
#include "fem/common/equation_base.h"
  
namespace SlaughterFEM{
	
/*! \class HeatEq heat_eq.h "fem/heat_eq/heat_eq.h"
 * \brief An equation class for the heat equation.
 * 
 * \see HeatEqSystem
 * \see example1.cpp
 * \see example2.cpp
 * \ingroup FEMTransientHeat
 */  
class HeatEq : public EquationBase<HeatEqSystem>{
	
	public:
		HeatEq(EquationSystems& es, const Order order = FIRST, const FEFamily family = LAGRANGE); 
		
		~HeatEq();
	
}; // class
}  // namespace
#endif
