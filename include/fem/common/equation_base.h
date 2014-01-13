/*! \file equation_base.h
  * \brief Header file for EquationBase class.
  */
  
// Limit multiple includes
#ifndef equation_base_h
#define equation_base_h 
  
// Standard C++ includes
#include <string>
using std::string;
 
// BOOST includes
//#include <boost/enable_shared_from_this.hpp>

// libMesh includes
#include <libmesh.h>
#include <equation_systems.h>
#include <system.h> 
 
// Add this class to my FEM namespace 
namespace SlaughterFEM{
 
/*! \class EquationBase equation_base.h "fem/common/equation_base.h"
 * \brief A base class for wrapping libMesh systems.
 * 
 * This class offers a base for creating a class that acts like a
 * libMesh::System. It automatically adds the system to the EquationSystems
 * variable provided. It then provides access to the system using the
 * system() member function or the operator(). 
 * 
 * This class was designed to offer the ability to manipulate a libMesh::System
 * without accessing it through the EquationSystems object, which requires
 * the name. 
 * 
 * This may be of little use for simple, single equation systems. But, becomes
 * useful for multi-equation systems that depend on each other. This then 
 * allows a pointer to be established that access the system, rather than
 * a reference.
 * 
 * \see HeatEq
 * \see ThermoEq
 * \see EnergyEq
 * \see example1.cpp
 * 
 * \ingroup FEMcommon
 */  	
template <class SystemType> class EquationBase{

	public:
		//! Class constructor
		/*!
		 * \param es The EquationSystems that the System will be added
		 * \param name The name of the system being added
		 */ 
		EquationBase(EquationSystems& es, const string name) : _name(name), _es(es){
			_es.add_system<SystemType>(name);
		}

		//! Class destructor
		~EquationBase(){};
		
		//! Access to the system via the () operator.
		/*!
		 * \see example1.cpp
		 */ 
		SystemType& operator()(){
			return _es.get_system<SystemType>(_name);
		}
		
		//! Access to the system via a member function
		SystemType& system(){
			return _es.get_system<SystemType>(_name);
		}
		
		//! Initialization function
		/*!
		 * Calls the base systems init() function by default, but 
		 * it is virtual allowing for the user to customize the behavior.
		 * \see FrontVelocityEq
		 */
		virtual void init(){
			this->system().init();
		}  

	protected:
	
		//! Storage location for the system name
		const string _name;
		
		//! Storage for the EquationSystems reference
		EquationSystems& _es;
		
}; // class
} // namespace
#endif
