 /*! \file energy_system.h
  * \brief Header file for the volume averaging EnergySystem class.
  */
    
 /*! \page volavg_energy Energy Equation
  * \dontinclude energy_system.cpp
  * \tableofcontents
  * 
  * The following some details of the volume averaged energy 
  * equation, which defined in the EnergyEq class.
  * The code is an implementation of the work presented in Samanta and 
  * Zabaras (2005). This equation as a single unknown, the enthalpy 
  * (\f$h\f$), which is created in the class constructor of EnergyEq class.
  * 
  * \section Constructor
  * The constructor does two simple tasks:
  * -# It defines the default
  * behavior of the initialization function, which assumes that the user
  * supplied function for initialization supplies a temperature. 
  * -# The enthalpy variable is defined. This variable must be named "h" or 
  * the program will not function correctly.
  * 
  * \par
  * \snippet energy_eq.cpp constructor 
  * 
  * \section Initialization
  * \snippet energy_system.cpp initialize 
  * 
  * \section Assembly
  * \subsection Preparation
  * This section explains the definition of the numerous variables that
  * must be established before the actual assembly can begin.
  * 
  * \skip void EnergyEq :: assemble()
  * \until QGauss
  * 
  */  
  
// Limit multiple includes
#ifndef volume_average_energy_system_h
#define volume_average_energy_system_h

// BOOST includes
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp> 

// LibMesh includes
#include <libmesh.h>
#include <libmesh_common.h>
#include <point.h>
#include <elem.h>
#include <mesh.h>
#include <system.h>
#include <nonlinear_implicit_system.h>
#include <equation_systems.h>
#include <transient_system.h>
#include <fe.h>
#include <fe_type.h>
#include <dof_map.h>
#include <quadrature_gauss.h>
#include <sparse_matrix.h>
#include <numeric_vector.h>
#include <type_vector.h>
#include <dense_matrix.h>
#include <dense_vector.h>
#include <perf_log.h>

#include <petsc_linear_solver.h>

// Use the standard LibMesh names space
using namespace libMesh;

// Standard library includes
#include <vector>
#include <string>
using std::string;
using std::vector;

// Includes from this library
#include "fem/common/system_base.h"
#include "fem/volume_average/thermo_system.h"
#include "fem/volume_average/momentum_system.h"

//#include "fem/common/eq_base.h"
#include "fem/common/my_analytic_function.h"
#include "fem/common/boundary_base.h"


// Add to FEM namespace
namespace SlaughterFEM{

//! A short-hand for the EqVariableLinker class
//typedef EqLink<TransientNonlinearImplicitSystem> Linker;

/*! \class EnergySystem energy_system.h "fem/volume_average/energy_system.h"
 * \brief A class for solving the volume averaged energy equation with libMesh
 * \ingroup FEMVolumeAverage
 */  	
class EnergySystem : public SystemBase<TransientNonlinearImplicitSystem>{
	
public:	

		// Calls to initialized() are from EqCore class
		//using EqBaseLink::initialized;
		
		// List inherited properties this class uses internally
		//using EqBase<TransientNonlinearImplicitSystem>::get_variable;
	
		//! Class constructor
		/*!
		 * \param es A reference to the EquationSystems class to be used
		 * \param name The name of the system
		 * \param number The number of the system
		 * 
		 * Should be added using EquationSystems::add_system
		 * 
		 * \see energy_system.cpp
		 */
		EnergySystem(EquationSystems& es, const string& name, const unsigned int number);
	
		//! Utilize a enthalpy initialization function instead of temperature
		/*!
		 * By default the EnergyEq class expects an initial temperature
		 * function add using the add_initial_function member. From the
		 * temperature the enthalpy is computed. Calling
		 * this function changes this behavior to bypass the conversion
		 * and use an enthalpy initilization function directly.
		 * 
		 * This can not be reversed.
		 */ 
		void enthalpy_initialization();

		void solve();

		//! A smart pointer to the ThermoEq class
		boost::shared_ptr<ThermoSystem> thermo;
		
		//! A smart pointer to the MomentumEq class
		boost::shared_ptr<MomentumSystem> momentum;

	private:	
		//! Initialization flag (false = temperature; true = enthalpy)
		bool enthalpy_init_;
		
		//! Converts temperature to enthalpy, Eq. 12
		/*!
		 */ 
		void temperature_to_enthalpy(DenseVector<Number>& output, const Point& p, const Real);
		
		//! Compute the alpha term of Eq. 69
		/*!
		 * 
		 */
		Number alpha(Gradient grad_T, Gradient grad_h, Number f);
		
		//! Intialization function
		virtual void initialize();
		
		//! libMesh assembly function
		void assemble();
				
		void update_rhs();
		
		
}; // EnergyEq class

} // namespaces
#endif 

