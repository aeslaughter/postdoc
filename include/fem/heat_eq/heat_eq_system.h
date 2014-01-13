 /*! \file heat_eq_system.h
  * \brief Header file for HeatEqSystem class.
  */
  
// Limit multiple includes
#ifndef heat_eq_system_h
#define heat_eq_system_h

// BOOST includes
#include <boost/shared_ptr.hpp>

// LibMesh: general includes
#include <libmesh.h>
#include <libmesh_common.h>
#include <point.h>
#include <string_to_enum.h>
#include <getpot.h>

// LibMesh: mesh includes
#include <mesh.h>
#include <boundary_info.h>
#include <dirichlet_boundaries.h>
#include <analytic_function.h>
#include <mesh_refinement.h>
#include <parmetis_partitioner.h>
#include <error_vector.h>
#include <kelly_error_estimator.h>

// LibMesh: system includes
#include <linear_implicit_system.h>
#include <equation_systems.h>
#include <transient_system.h>
#include <vector_value.h>

// libMesh: Define finite element object and degree of freedom object
#include <fe.h>
#include <fe_type.h>
#include <dof_map.h>

// libMesh: Define Gauss quadrature rules
#include <quadrature_gauss.h>

// libMesh: Datatypes for finite element matrix and vector components.
#include <sparse_matrix.h>
#include <numeric_vector.h>
#include <dense_matrix.h>
#include <dense_vector.h>
#include <elem.h>

// libMesh: Allows for monitoring performance
#include <perf_log.h>

// Use the standard LibMesh names space
using namespace libMesh;

// Standard library includes
#include <vector>
#include <string>
//#include <memory>

// Includes from this library
#include "fem/common/my_analytic_function.h"
#include "fem/common/implicit_system_base.h"
#include "fem/common/boundary_base.h"
#include "fem/heat_eq/heat_eq_boundary.h"

// Add to FEM and HeatEq namespaces
namespace SlaughterFEM{

/*! \class HeatEqSystem heat_eq_system.h "fem/heat_eq/heat_eq_system.h"
 * \brief A class for solving the heat equation with libMesh
 * 
 * Provides mechanims for defining a libMesh equation 
 * system of the heat equation including, adding dirichlet, neumann, and convection
 * boundary conditions, including a heat source term, and
 * applying adaptive mesh refinement.
 * 
 * \b Strong \b Form:
 * \f[
 * \rho c_p \frac{\partial T}{\partial t} - \nabla\cdot\mathbf{q} + s = 0,
 * \f]
 * where \f$ t \f$ is time, \f$ \rho \f$ is density, \f$ c_p \f$ is specific heat, \f$ T \f$
 * is temperature, \f$ \mathbf{q} \f$ is heat flux vector, and \f$ s \f$ is the heat
 * source.
 * 
 * \b Weak \b Form:
 * \f[
 * \int_\Omega w^T \rho c_p \frac{\partial T}{\partial t} d\Omega -
 * \int_\Omega \nabla w^T \cdot\mathbf{q} d\Omega +
 * \int_\Omega w^T s d\Omega + 
 * \int_{\Gamma_q} w^T \overline{\mathbf{q}} d\Gamma +
 * \int_{\Gamma_h} w^T h T d\Gamma - 
 * \int_{\Gamma_h} w^T h T_{\infty} d\Gamma = 0,
 * \f]
 * where \f$ \Omega \f$ defines the entire doman and \f$ \Gamma \f$ defines
 * the boundaries. The subscripts \f$ q \f$ and \f$ h \f$ for the boundary
 * intergrals indicate the known heat flux (\f$ \overline{\mathbf{q}} \f$)
 * boundary and the convective boundary (\f$ q = h(T - T_\infty) \f$), respectively.
 * \f$ w \f$ is the test function.
 * 
 * \ingroup FEMTransientHeat
 */  	
class HeatEqSystem : public ImplicitSystemBase<TransientLinearImplicitSystem, HeatEqBoundaryBase>{

	public:	
		//! Class constructor
		/*!
		 * \param es A reference to the EquationSystems class to be used
		 * \param name The name of the system
		 * \param number The number of the system
		 * 
		 * Should be added using EquationSystems::add_system
		 *
		 * \see heat_eq.cpp
		 */
		HeatEqSystem(EquationSystems& es, const string& name, const unsigned int number);
		
		//! Class destructor
		~HeatEqSystem();

	protected:	
		//! libMesh assembly function
		void assemble();
		
		//! Method for constructing the general time integration stiffness matrix and RHS vector
		void build_stiffness_and_rhs(DenseMatrix<Number>& K_hat, 
			DenseVector<Number>& F_hat, 
			DenseMatrix<Number> Me, 
			DenseMatrix<Number> Ke, 
			DenseVector<Number> Fe_old, 
			DenseVector<Number> Fe, 
			DenseVector<Number> u_old, 
			Real dt, 
			Real theta);
				
}; // Heat_eq class

} // namespaces
#endif
