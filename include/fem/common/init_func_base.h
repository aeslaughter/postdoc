/*! \file init_func_base.h
  * \brief A base class for system initilization
  */ 

// Limit multiple includes
#ifndef init_func_base_h
#define init_func_base_h

// Standard includes
#include <string>
using std::string;
  
// libMesh includes
#include <libmesh.h>
#include <libmesh_common.h> 
#include <auto_ptr.h>
#include <dense_vector.h>
#include <point.h>
#include <function_base.h>

// Use libMesh naming
using namespace libMesh;  
  
// Add this to my FEM library namespace
namespace SlaughterFEM{

/*! \class InitFuncBase init_func_base.h "fem/common/init_func_base.h"
 * \brief A base class defining initilization functions.
 * \tparam Output The desired type of output (defaults to a libMesh::Number)
 * \ingroup FEMcommon
 */  	
template <typename Output = Number>
class InitFuncBase : public FunctionBase<Output>{
	public:
	
	//! Class constructor
	InitFuncBase(){}
	
	//! Class destructor
	~InitFuncBase(){}
	
	//! Returns the desired value for the given variable at the specified point
	/*!
	 * \param index The numeric reference to the variable of interest
	 * \param p The current point at which to compute the value of the desired variable
	 * \return The value of interest for the given variable and point
	 */ 
	virtual Number component(unsigned int index, const Point& p, Real time = 0) = 0;
		
	//! Create a copy of the class
	//! \todo Test if this is adequate for derived classes
	virtual AutoPtr<FunctionBase<Output> > clone () const{
		return AutoPtr<FunctionBase<Output> > (new InitFuncBase<Output>());
	}	
	
	//! Not implemented
	Output operator() (const Point& p, const Real t = 0.){
		libmesh_not_implemented();
	}

	//! Not implemented
	void operator() (const Point& p, const Real t, DenseVector<Output>& output){	
		libmesh_not_implemented();		
	}
};
}
#endif
