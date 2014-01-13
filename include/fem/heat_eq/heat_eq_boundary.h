/*! \file heat_eq_boundary.h
  * \brief Header file for Heat_eq_boundary classes
  */

// Avoid multiple includes
#ifndef heat_eq_boundary_h
#define heat_eq_boundary_h

// Standard library includes
#include <vector>
#include <string>

// Boost
#include <boost/function.hpp>

// LibMesh includes
#include <libmesh.h>
#include <dense_vector.h>
#include <point.h>

// My includes
#include "fem/common/boundary_base.h"

// Use the libMesh namespace
using namespace libMesh;

// Include this in my namespace
namespace SlaughterFEM{

/*! \class HeatEqBoundaryBase heat_eq_boundary.h "fem/heat_eq/heat_eq_boundary.h"
 * \brief A base class for heat equation boundary conditions
 * 
 * Provides basic format for defining a boundary conditions for application
 * to the HeatEq class. This class is used to provide uniform access
 * for all types of boundaries that are dervied from the base, as is 
 * in the case for the HeatEq class.
 * 
 * \see HeatEqBoundaryDirichlet
 * \see HeatEqBoundaryNeumann
 * \see HeatEqBoundaryConvection
 * 
 * \ingroup FEMTransientHeat
 */  	
class HeatEqBoundaryBase : public BoundaryBase{
	public:		
		//! The constructor for the base class of heat equation boundaries
		/*! \param the_type A string that indicates the type of boundary
		 *
		 * Three types are currently supported:
		 * - "dirichlet"
		 * - "convection"
		 * - "neumann"
		 */ 
		HeatEqBoundaryBase(std::string the_type);
		
		//! A function for returning temperature for dirichlet condition
		/*! 
		 * \param p A libMesh Point
		 * \param t The current time
		 */ 
		virtual Number T(const Point&, const Real){
			return 0;
		};
		
		//! A function for returning convection term, h
		/*! 
		 * \param p A libMesh Point
		 * \param t The current time
		 */ 
		virtual Number h(const Point&, const Real){
			return 0;
		};

		//! A function for returning convection term, Tinf
		/*! 
		 * \param p A libMesh Point
		 * \param t The current time
		 */ 
		virtual Number Tinf(const Point&, const Real){
			return 0;
		};
		
		//! A function for returning neumann condition term, q
		/*! 
		 * \param p A libMesh Point
		 * \param t The current time
		 */ 
		virtual Real q(const Point&, const Real){ 
			return 0;
		};
};

/*! \class HeatEqBoundaryDirichlet heat_eq_boundary.h "fem/heat_eq/heat_eq_boundary.h"
 * \brief A class for dirichlet boundary conditions
 * \ingroup FEMTransientHeat
 */ 
class HeatEqBoundaryDirichlet : public HeatEqBoundaryBase{
	public:	
		//! Class constructor
		/*! Creates a dirichlet boundary conditions class, the constructor
		 * specifies the type as "dirichlet."
		 */ 
		HeatEqBoundaryDirichlet();
		
		//! A virtual function for returning the temperature
		/*! By default this function returns a constant value defined
		 * in the T_constant attribute. If constant temperature is to be 
		 * used then only T_constant must be changed (zero by default).
		 * To create a non-constant flux a derived class should be 
		 * created that inherts this class. This T() function should
		 * then be defined in this new class to return the desired 
		 * temperature value.
		 */ 
		virtual Real T(const Point& p, const Real t);
		
		//! A constant temperature value, used in default operation
		Number T_constant;		
	
	protected:	
	
		//! Computes and returns the boundary temperature
		/*!
		 * This returns the temperature value that is computed by the T function
		 * that may be customized. This member is used by the HeatEq
		 * class thus is protected.
		 * 
		 * \param output A reference to the vector that stores the
		 * stores the computed flux
		 * \param p A libMesh point
		 * \param t The time
		 * 
		 * \see T
		 */ 
		void value(DenseVector<Number>& output, const Point& p, const Real t){
			output(0) = T(p,t);
		}
};

/*! \class HeatEqBoundaryNeumann heat_eq_boundary.h "fem/heat_eq/heat_eq_boundary.h"
 * \brief A class for flux boundary conditions
 * \ingroup FEMTransientHeat
 */  	
class HeatEqBoundaryNeumann : public HeatEqBoundaryBase{
	public:
		//! Class constructor
		/*! Creates a neumann boundary conditions class, the constructor
		 * specifies the type as "neumann". It also sets q_constant 
		 * to zero.
		 */
		HeatEqBoundaryNeumann();
		
		//! A virtual function for returning the flux
		/*! By default this function returns a constant value defined
		 * in the q_constant attribute. If constant flux is to be used
		 * then only q_constant must be changed (it is zero by default).
		 * To create a non-constant flux a derived class should be 
		 * created that inherts this class. This q() function should
		 * then be defined in this new class to return the desired 
		 * flux value.
		 */ 
		virtual Real q(const Point& p, const Real t);
		
		//! A constant flux value, used in default operation
		Number q_constant;	
	
	protected:
	
		//! Computes and returns the heat flux
		/*!
		 * This returns the flux value that is computed by the q function
		 * that may be customized. This member is used by the HeatEq
		 * class thus is protected.
		 * 
		 * \param output A reference to the vector that stores the
		 * stores the computed flux
		 * \param p A libMesh point
		 * \param t The time
		 * 
		 * \see q
		 */ 
		void value(DenseVector<Number>& output, const Point& p, const Real t){
			output(0) = q(p,t);
		}
};

/*! \class HeatEqBoundaryConvection heat_eq_boundary.h "fem/heat_eq/heat_eq_boundary.h"
 * \brief A class for convection boundary conditions
 * \ingroup FEMTransientHeat
 */ 
class HeatEqBoundaryConvection : public HeatEqBoundaryBase{
	public:
		//! Class constructor
		/*! Creates a convection boundary conditions class, the constructor
		 * specifies the type as "convection". It also sets h_constant 
		 * and Tinf_constant to zero.
		 */
		HeatEqBoundaryConvection();
		
		//! A virtual function for returning the convection term, h
		/*! By default this function returns a constant value defined
		 * in the h_constant attribute. If constant h is to be used
		 * then only h_constant must be changed (it is zero by default).
		 * To create a non-constant h a derived class should be 
		 * created that inherts this class. This h() function should
		 * then be defined in this new class to return the desired 
		 * value.
		 */ 
		virtual Number h(const Point& p, const Real t);
		
		//! A virtual function for returning the convection term, Tinf
		/*! By default this function returns a constant value defined
		 * in the Tinf_constant attribute. If constant Tinf is to be used
		 * then only Tinf_constant must be changed (it is zero by default).
		 * To create a non-constant Tinf a derived class should be 
		 * created that inherts this class. This Tinf() function should
		 * then be defined in this new class to return the desired 
		 * value.
		 */ 
		virtual Number Tinf(const Point& p, const Real t);
		
		//! A constant h term, used in the default operation
		Number h_constant;
		
		//! A constant Tinf term, used in the default operation		
		Number Tinf_constant;
		
	protected:
	
		//! Computes and returns the convective term
		/*!
		 * This returns the convective term (\f$h \cdot T_{\infty}\f$) value 
		 * that is computed by the h and Tinf functions that may be 
		 * customized. This member is used by the HeatEq
		 * class thus is protected.
		 * 
		 * \param output A reference to the vector that stores the
		 * stores the computed flux
		 * \param p A libMesh point
		 * \param t The time
		 * 
		 * \see Tinf
		 * \see h
		 */ 
		virtual void value(DenseVector<Number>& output, const Point& p, const Real t){
			output(0) = h(p,t) * Tinf(p,t);
		}
};

} //namespaces
#endif
