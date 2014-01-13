/*! \file explicit_system_base.h
  * \brief Header file for ExplicitSystemBase class.
  */
  
// Limit multiple includes
#ifndef explicit_system_base_h
#define explicit_system_base_h

// Standard include
#include <string>
using std::string;

// libMesh includes
#include <libmesh.h>
#include <equation_systems.h>
#include <transient_system.h>
#include <explicit_system.h>
#include <numeric_vector.h>
#include <auto_ptr.h>

// Includes from this library
#include "fem/common/system_base.h"
#include "fem/common/my_analytic_function.h"
 
// Add this class to my FEM namespace
namespace SlaughterFEM{
	
/*! \class ExplicitSystemBase explicit_system_base.h "fem/common/explicit_system_base.h"
 * \brief A base class for defining nodal data
 * 
 * Provides mechanims for defining nodal data using a libmesh
 * system. The data may very spatially and temporally. In the inherited
 * class the pure virtual functions dictate the behavior.
 *
 * \ingroup FEMcommon
 */  		
class ExplicitSystemBase : public SystemBase<TransientExplicitSystem>{
	public:
		// The init() functin defined here is for the system
		using SystemBase<TransientExplicitSystem>::init;

		//! Pure virtual function for returning variable values at a point
		/*!
		 * This function must be defined in the inhereted class. It should
		 * return the value of each of the variables in the System.
		 * \param index The index of the variable to return
		 * \param p The Point at which the value is computed
		 * \param t The time, this should not be used because libMesh 
		 * 		forces time to zero in projections. If time is needed extract
		 * 		it from the System variable.
		 */ 
	//	virtual Number component(unsigned int index, const Point& p, Real t = 0) = 0;
			
		//! Pure virtual function for creating a clone of this class
		/*!
		 * This must be defined in the inherting class, it should return
		 * a fully function copy of the class itself, see VolumeAverageEqData
		 * for an example.
		 * 
		 * \see VolumeAverageEqData
		 */
	//	virtual AutoPtr<FunctionBase<Number> > clone () const = 0;
		
		//! Function for returing a complete vector of variables
		/*!
		 * As with the component function, this must be defined to return
		 * a vector (using output reference) of all of the variables in
		 * the system.
		 * 
		 * By default this loops through each variable and calls component.
		 *
		 * \param p The Point at which the value is computed
		 * \param t The time, this should not be used because libMesh 
		 * 		forces time to zero in projections. If time is needed extract
		 * 		it from the System variable.
		 * \param output A reference to the vector for storing the values
		 * 
		 * \see component
		 * \see VolumeAverageEqData
		 */ 
	//	virtual void operator()(const Point& p, const Real t, DenseVector<Number>& output);		
		
		//! A scalar output version of the operator()
		/*!
		 * This function should not be necessary, the component and 
		 * vector version of the operator() are used instead of this
		 * when data is projected.
		 * 
		 * \param p The Point at which the value is computed
		 * \param t The time, this should not be used because libMesh 
		 * 		forces time to zero in projections. If time is needed extract
		 * 		it from the System variable.
		 * 
		 * \see component
		 * \see operator()(const Point& p, const Real t, DenseVector<Number>& output) 
		 */ 
	//	virtual Number operator() (const Point& p, const Real t = 0.);

		//! Initializes the equation system 
		/*! 
		 * Basic initilization class, by default it simply calls the 
		 * System intialization function.
		 * 
		 * \param t_initial The start time, defaults to zero.
		 * \todo There is an init() for both the System and FunctionBase, does this cause problems
		 */ 
	//	virtual void init(Real t_initial = 0);

		//! Projects the nodal data at the prescribed time
		/*!
		 * This function creates a boost::function pointer to the 
		 * value member of this class, packages that pointer into a 
		 * libmesh acceptable format with MyAnalyticFunction, and then
		 * uses it to project the nodal data.
		 * 
		 * \param t A libmesh::Real containing the current time, this 
		 * time is used to redefine the TransientExplicitSystem time
		 * that contains the nodal data, thus it may be used for 
		 * temporal data or it may just be ignored by letting it default
		 * to zero.
		 */ 
		virtual void update_solution(Real t) = 0;
		
	protected:
	
		//! Class constructor
		/*!
		 * This class is meant to be inherited, as such the constructor
		 * is protected.
		 * 
		 * The class requires that an existing EquationSystem object
		 * be passed in, the system with the name given by \c data_name
		 * is added to the EquationSystem object.
		 * 
		 * \param sys libmesh::EquationSystems object that the new system
		 * for storing data will be added
		 * \param data_name The name given the the system for storing data
		 * 
		 * \see example3.cpp
		 */ 
		ExplicitSystemBase(EquationSystems& es, const string& name, const unsigned int number);
		
		//! Class destructor
		~ExplicitSystemBase(){}
		
		//! Initilization function
		/*!
		 * This function is called by libMesh when initilizing the 
		 * System via the init() function.
		 */
		virtual void initialize() = 0;
		
}; // ends class
}  //ends namespace
#endif
