/*! \file boundary_base.h
  * \brief A base class for my boundary condition class
  */ 
  
// Avoid multiple includes
#ifndef boundary_base_h
#define boundary_base_h  
  
// C++ Standard includes
#include <vector>
using std::vector;  
  
// Boost includes
#include <boost/function.hpp>

// LibMesh includes
#include <libmesh.h>
#include <dense_vector.h>
#include <point.h>

// Include this in my namespace
namespace SlaughterFEM{

/*! \class BoundaryBase boundary_base.h "fem/equation_systems/boundary_base.h"
 * \brief A base class for equation boundary conditions for integration with EqBase class.
 * 
 * Provides basic format for defining a boundary conditions for application
 * to the ImplicitSystemBase class. This class is used to provide uniform access
 * for all types of boundaries that are derived from the base.
 * 
 * \see ImplicitSystemBase
 * \see HeatEqBoundaryBase
 * 
 * \ingroup FEMcommon
 */  	
class BoundaryBase{
	public:		
		//! The constructor for the base class
		/*! 
		 * \param the_type A string that indicates the type of boundary
		 * 
		 * Note that if the inherited class is to be recognized by the 
		 * EqBase class as a dirichlet condition this value MUST be
		 * "dirichlet" as done with the HeatEqBoundaryDirichlet class.
		 * 
		 * \see HeatEqBoundaryDirichlet
		 */ 
		BoundaryBase(std::string the_type);
		
		//! The boundary id
		int id;	
		
		//! Vector containing the variable indices to apply the condition to
		vector<unsigned int> variables;
		
		//! A name that specifies the boundary type
		const std::string type;		
	
		//! A function pointer (libMesh required format)
		void (*fptr)(DenseVector<Number>&, const Point&, const Real);	

		//! A pure virtual function for returning the desired value	
		/*!
		 * For any class derived from this base class the value function
		 * must be defined. The value returned through the DenseVector
		 * reference is used by the Eq_base class when applying the 
		 * boundary conditions.
		 * 
		 * \see Eq_base
		 */ 
		virtual void value(DenseVector<Number>& output, const Point&, const Real) = 0;
		
}; // class
}  // namespace
#endif
