 /*! \file concentration_system.h
  * \brief Header file for the volume averaging ConcentrationEq class.
  */
  
// Limit multiple includes
#ifndef volume_average_concentration_system_h
#define volume_average_concentration_system_h

// LibMesh includes
#include <libmesh.h>
#include <libmesh_common.h>
#include <point.h>
#include <elem.h>
#include <mesh.h>
#include <nonlinear_implicit_system.h>
#include <equation_systems.h>
#include <transient_system.h>
#include <fe.h>
#include <fe_type.h>
#include <dof_map.h>
#include <quadrature_gauss.h>
#include <sparse_matrix.h>
#include <numeric_vector.h>
#include <dense_matrix.h>
#include <dense_vector.h>
#include <perf_log.h>

// Use the standard LibMesh names space
using namespace libMesh;

// Standard library includes
#include <vector>
#include <string>
using std::string;
using std::vector;

// Boost
#include <boost/shared_ptr.hpp>

// Includes from this library
#include "fem/common/implicit_system_base.h"
#include "fem/common/my_analytic_function.h"
#include "fem/common/boundary_base.h"

// Add to FEM namespace
namespace SlaughterFEM{	

/*! \class MomentumEq momentum_eq.h "fem/volume_average/momentum_eq.h"
 * \brief A class for solving the volume averaged momentum equation with libMesh
 * \ingroup FEMVolumeAverage
 */  	
class ConcentrationSystem : public ImplicitSystemBase<TransientNonlinearImplicitSystem>{

	public:	
		//! Class constructor
		/* \param sys A reference to the EquationSystems class to be used
		 * \param order The finite element order to use
		 * 
		 * \see momentum_eq.cpp
		 */
		ConcentrationSystem(EquationSystems& es, const string& name, const unsigned int number);

	private:	
		//! libMesh assembly function
		void assemble();
};

} // namespaces
#endif
