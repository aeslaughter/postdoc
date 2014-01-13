/*! \file system_base.h
  * \brief Header file for SystemBase class.
  */
  
// Limit multiple includes
#ifndef system_base_h
#define system_base_h 
  
#include <string>
using std::string;

// BOOST includes
#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp> 

// LibMesh includes
#include <libmesh.h>
#include <libmesh_common.h>
#include <vector_value.h>
#include <point.h>
#include <parameters.h>
#include <equation_systems.h>
#include <transient_system.h>
#include <dof_map.h>

// My FEM includes
#include "fem/common/my_analytic_function.h"
#include "fem/common/init_func_base.h"

// Define namespaces used
using namespace libMesh;

namespace SlaughterFEM{

/*! A short-hand method for defining a function pointer for the initialization function (libmesh method)
 * \ingroup FEMcommon
 */
typedef libMesh::Number (*system_base_init_func_libmesh)(const Point& p, const Parameters& parameters, const std::string& sys_name, const std::string& unknown_name);

/*! A short-hand method for defining a function pointer for the initialization function (my boost::function method)
 * \ingroup FEMcommon
 */
typedef libMesh::Number (*system_base_init_func_boost)(DenseVector<Number>&, const Point&, const Real);

//! \todo This should be named TransientSystemBase
/*! \class SystemBase SystemBase.h "fem/common/system_base.h"
 * \brief My base class for libMesh transient systems
 * 
 * Provides mechanism for defining a libMesh transient systems. 
 * This class must be inherited and the pure virtual functions must be 
 * defined in the inherited class.
 * 
 * The public members are designed to be used by the use in their main
 * program. The protected members, although some must be defined, should
 * not be needed by the user once defined.
 * 
 * \ingroup FEMcommon
 */  	
template <class TheSystem> class SystemBase: 
	public TheSystem,					// the type of System to inherit (e.g., ExplicitSystem)
	public System::Initialization{		// allows for libMesh initialization class{
	
	public:	
	
		//! Template function for extracting an EquationSystem parameter
		/*!
		 * \tparam ParamType The type of parameter being extracted (e.g., double)
		 * \param name The name of the parameter (e.g., "density")
		 * \return The value for the specified parameter
		 * 
		 * The primary reason for this is for future expansion to 
		 * allow for parameters that vary with space and time.
		 */ 
		template<typename ParamType> ParamType get_constant(string name){  
			ParamType output = this->get_equation_systems().parameters. template get<ParamType>(name);
			return output;
		}

		//! Template function for setting an EquationSystem parameter
		/*!
		 * \tparam ParamType The type of parameter being added (e.g., double)
		 * \param name The name of the parameter (e.g., "density")
		 * \param var The value of the parameter (e.g., 123) 
		 * \return The specified value for the parameter
		 * 
		 * The primary reason for this is for future expansion to 
		 * allow for parameters that vary with space and time.
		 */ 
		template<typename ParamType> ParamType set_constant(string name, ParamType var){  
			this->get_equation_systems().parameters. template set<ParamType>(name) = var;
			return var;
		}
	
		//! Add a function pointer to the initialization function
		/*! 
		 * This version uses the libmesh documented method for adding 
		 * initial conditions for an equation system.  It must 
		 * be specified in a function with the format given below.
		 * (see libMesh help for input format).
		 * 
		 * \param func A function pointer to the initialization function
		 */ 
		void add_initial_function(system_base_init_func_libmesh func){
			init_func = func;		// using traditional function pointer
			bst_init_func = NULL;
		}
		
		//! Add function pointer to act as initialization function using boost::function indirectly
		/*! 
		 * The initial conditions for the the equation must 
		 * be specified, this version takes the function pointer and
		 * converts it to a boost::function and uses the overloaded 
		 * version of the \c add_initial_function to apply the function
		 * to the equation system.
		 * 
		 * \param func A function pointer to the initialization function
		 */ 
		void add_initial_function(system_base_init_func_boost func){
			
			// Create a boost::function pointer
			boost::function<void (DenseVector<Number>& output, const Point&, const Real)> init_ptr;
			
			// Apply the supplied function to the boost::function pointer
			init_ptr = func;
			
			// Call the version of this function that uses boost::function
			add_initial_function(init_ptr);
		}
			
		//! Add a boost::function pointer act as initialization function
		/*! 
		 * The initial conditions for equation must 
		 * be specified, this function may be passed to the EqBase
		 * class using a boost::function which allows for extreme 
		 * flexibility when binding is used.
		 * 
		 * \param bst_fptr A boost::function that points to the 
		 * initialization function or class
		 */ 
		void add_initial_function(boost::function<void (DenseVector<Number>& output, const Point&, const Real)> bst_fptr){
			bst_init_func = bst_fptr;	// using boost::function
			init_func = NULL;
		}

		//! Add initial conditions using the EqInitBase class
		/*!
		 * Using the EqInitBase class allows for a more effiencent
		 * handling of initilizing a system that has multiple variables,
		 * such as velocity. When libMesh projects a solution it calls
		 * the initilization function for each variable, thus if a
		 * traditional function is called that sets multiple variables
		 * it will be called for each variables.
		 * 
		 * This calling is done via the compoment member of the libMesh::
		 * FunctionBase class. EqInitBase allows the user to reimplement
		 * this component member so that it only performs the calculation
		 * one piece at a time.
		 *
		 * \param ptr A boost::shared_ptr an EqInitBase derived class.
		 */  
		
		void add_initial_object(boost::shared_ptr<InitFuncBase<Number> > ptr){
			init_ptr = ptr;
			bst_init_func = NULL;	
			init_func = NULL;
		}
		
		//! Updates the solution with time
		/*! 
		 * Passes the current solution to the old solution
		 * 
		 * \param time The current time
		 * \param dt The change in time
		 */ 
		//virtual void update_solution(Real time, Real dt){
		//	*this->old_local_solution = *this->current_local_solution;
		//}	
		
		// A function for returning a value at a point given the name
		/*!
		 * \param name The name of the variable to extract
		 * \param p The Point at which the values are desired
		 */
		Number point_value(const string name, const Point& p){
			unsigned int idx = System::variable_number(name);
			return System::point_value(idx, p);
		}
		
		//! A function for getting a vector of the variables at a point
		/*!
		 * This function simply employs the libMesh::System point_value
		 * member function for each of the variables. This function is 
		 * virtual so that derived class can return specific variables.
		 *
		 * \param p The Point at which the values are desired
		 * 
		 * \see MomentumEq
		 */  	
		virtual VectorValue<Number> point_value(const Point& p){
			
			// Get a reference to the system object
			unsigned int n = this->n_vars();
			
			// Collect the velocity vector components
			VectorValue<Number> var;
			for (unsigned int d = 0; d < n; d++){
				var(d) = TheSystem::point_value(d, p);
			}
			
			// Return the vector containing desired variables
			return var;
		}
	
	//! \todo remove this, it is not really needed
		//! Return the number of dimensions
		/*!
		 * A function for returning the number of dimensions for the 
		 * system. This is simply a short-cut that by passes calling
		 * getting the mesh and then dimensions from the system.
		 * 
		 * \return The dimensionality of the system
		 */
		unsigned int ndim(){
			return _ndim;
		}  

		//! Return the initilization status
		/*!
		* Derived classes can control the initilized_ flag, this function
		* is a public method to return the value of this flag.
		*/
		virtual bool initialized(){
		  return _initialized;
		} 
	
	protected:
	
		SystemBase(EquationSystems& es, const string& name, const unsigned int number) : 
			TheSystem(es, name, number), _ndim(es.get_mesh().mesh_dimension()){
			
			// Test that the template class has System as a base
			try{
				dynamic_cast<System*>(this);
			} catch (int e){
				std::cout << "The template class supplied to the SystemBase class must be derived from libMesh::System." << e << std::endl;
				libmesh_error();
			}			
			
			// Attach itself as initialization object
			attach_init_object(*this);
			
			// The system has not been initilized
			_initialized = false;
		}
		
		~SystemBase(){}
		
		//! libMesh initilize function
		virtual void initialize(){

			// Project the solution using boost::function
			if(bst_init_func != NULL){
				MyAnalyticFunction<Number> func_object(bst_init_func);
				this->project_solution(&func_object);

			// Project the solution using a function pointer directly
			} else if (init_func != NULL){
				this->project_solution(init_func, NULL, TheSystem::get_equation_systems().parameters);
				
			// Use the EqInitBase class (uses component() member)	
			} else if (init_ptr){
				this->project_solution(init_ptr.get());
				
			} else {
				printf("ERROR: No initilization function exists\n");
				exit(1);
			}
			
			// Update the initialization state
			_initialized = true;	
		}
		
		//! The number of dimensions
		const unsigned int _ndim;			
					
		//! Initiliation status
		bool _initialized;	
		
		//! A shared pointer to the EqInitBase class
		boost::shared_ptr<InitFuncBase<Number> > init_ptr;
		
		//! A function pointer to the initialization function		 
		system_base_init_func_libmesh init_func;
		
		//! A boost::function for initialization 
		boost::function<void (DenseVector<Number>& output, const Point&, const Real)> bst_init_func;
		
}; // class
}  // namespace
#endif
