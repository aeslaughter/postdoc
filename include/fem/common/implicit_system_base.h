/*! \file implicit_system_base.h
  * \brief A base class for my ImplicitSystemBase wrappers
  */ 

// Avoid multiple includes
#ifndef implicit_system_base_h
#define implicit_system_base_h

// C++ Standard includes
#include <string>
#include <vector>
using std::string;
using std::vector;

// Boost includes
#include <boost/shared_ptr.hpp>
#include <boost/bind.hpp>

// LibMesh: general includes
#include <libmesh.h>
#include <libmesh_common.h>
#include <parameters.h>
#include <point.h>
#include <string_to_enum.h>
#include <getpot.h>

// LibMesh: mesh includes
#include <mesh.h>
#include <mesh_generation.h>
#include <boundary_info.h>
#include <dirichlet_boundaries.h>
#include <analytic_function.h>
#include <mesh_refinement.h>
#include <parmetis_partitioner.h>
#include <error_vector.h>
#include <kelly_error_estimator.h>

// LibMesh: system includes
#include <linear_implicit_system.h>
#include <nonlinear_implicit_system.h>
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

// Use the libMesh namespace
using namespace libMesh;

// Includes from this library
#include "fem/common/system_base.h"
#include "fem/common/my_analytic_function.h"
#include "fem/common/boundary_base.h"
#include "fem/common/init_func_base.h"

// Add this to the FEM namespace
namespace SlaughterFEM{

/*! \class ImplicitSystemBase system_base.h "fem/common/system_base.h"
 * \brief A template base class for using libmesh to solve equations
 * 
 * Provides mechanism for defining a libMesh equation system. It utilizes
 * the EqBoundaryBase class for implementing boundary conditions.
 * 
 * This class must be inherited and the pure virtual function must
 * be defined in the inherited class.
 * 
 * The public members are designed to be used by the use in their main
 * program. The protected members, although some must be defined, should
 * not be needed by the user once defined.
 * 
 * \see HeatEq
 * \see example1.cpp
 * \see example2.cpp
 * 
 * \ingroup FEMcommon
 */  	
template <typename Type, class TypeBoundaryBase = BoundaryBase> class ImplicitSystemBase :
	public SystemBase<Type>,		// inherit the SystemBase class
	public System::Assembly{		// allows for libMesh initialization class
	
	public:

		 /* 
		 * Note, when using set_constant/get_constant from this class
		 * the template term must be used:
		 * 		EqCore<Type>::template set_constant<Real>("time", 0);
		 * This is not required in classes who inherit ImplicitSystemBase and 
		 * specify a template parameter as in HeatEq and EnergyEq.
		 */ 
	
		//! Class constructor
		/*!
		 * The ImplicitSystemBase class is designed to be inherited, thus the 
		 * constructor is protected.
		 * 
		 * When inheriting this class it is import to explicitly call
		 * this constructor in the constructor of the inheriting class, 
		 * otherwise the class will not function properly. See the 
		 * HeatEq class for an example.
		 * 
		 * \tparam Type The type of libmesh system to be added, e.g.,
		 * TransientLinearImplicitSystem as in the HeatEq class.
		 * \tparam TypeBoundaryBase The base class for the boundary conditions.
		 * By default it uses EqBoundaryBase but the user might what to
		 * change this to use an inherited class of EqBoundaryBase, 
		 * again as done in the HeatEq class.
		 * 
		 * \param sys A libmesh EquationSystems to work from
		 * \param name The name of the system being created
		 * \tparam Type The type of libMesh::System (e.g. TrasientLinearImplicitSystem)
		 * \tparam TypeBoundaryBase The base class for boundary condition class,
		 * defaults to BoundaryBase
		 * \see HeatEq
		 */ 
		ImplicitSystemBase(EquationSystems& es, const string& name, const unsigned int number) : 
			SystemBase<Type>(es, name, number){
			
			// By default mesh refinement is not used
		//	using_mesh_refinement = false;
			   
			// Attach an initilization function for the heat equation
			//this->attach_init_object(*this);	
			
			// Attach the assembly member function
			//System::attach_assemble_object(*this);
			
			// System initilization state
			//SystemBase<Type>::_initialized = false;
		}
		
		//! Initializes the equation system 
		/*! 
		 * This function must be called before the equation is solved,
		 * but after all of the boundaries are defined. It does two 
		 * things, adds the dirichlet conditions and then calls the 
		 * libMesh equation systems init() function.
		 */ 
		virtual void init(Real t_initial = 0){		

			// Add a time variable and set system time
			SystemBase<Type>::template set_constant<Real>("time", t_initial);
			this->time = t_initial;

			// Apply Dirichlet boundary
			apply_dirichlet();

			// Initialize
			System::init();		
		}
		
		//! Updates the solution with time
		/*! 
		 * Passes the current solution to the old solution and reapplies
		 * the dirichlet boundary constraints 
		 * 
		 * \param time The current time
		 * \param dt The change in time
		 */ 
		virtual void update_solution(Real time, Real dt){
		
			// Apply the new time and time step
			this->time = time;
			SystemBase<Type>::template set_constant<Real>("dt", dt);
			SystemBase<Type>::template set_constant<Real>("time", time);

			// Step the system with time
			*SystemBase<Type>::old_local_solution = *SystemBase<Type>::current_local_solution;
			
			// Update the contraints (dirichlet boundaries)
			const MeshBase& mesh = this->get_mesh();
			DofMap& dof_map = this->get_dof_map();
			dof_map.create_dof_constraints(mesh, time);
		}	

		//! Adds a boundary object
		/*! 
		 * A template class for adding a boundary condition
		 * 
		 * \tparam TypeNewBoundary The type of boundary being added. The
		 * type added must be derived from EqBoundaryBase. See example1.cpp
		 * for an example.
		 * 
		 * \param id The integer identification for the boundary
		 * \param var A vector of integers indentifiy which variables to apply the boundary condition
		 * 
		 * \see EqBoundaryBase
		 * \see example1.cpp
		 */ 
		template<class TypeNewBoundary> boost::shared_ptr<TypeNewBoundary> add_boundary(int id, vector<unsigned int> var = vector<unsigned int>(1,0)){
			
			// Check that the enough variable exists
			for (unsigned int i = 0; i < var.size(); i++){
				if (var[i] > this->n_vars()){
					printf("ERROR: The index specified for application to the boundary is out of range of the existing variables for the current system.\n");
				}
			}
			
			// Check that the desired boundary id is unique
			for(unsigned int i = 0; i < bc_ptrs.size(); ++i){
				if(bc_ptrs[i]->id == id){
					printf("ERROR: The boundary id of %d was already used, id's must be unique!\n", id);
					exit(1);	
				}
			}
			
			// Create the boundary via a shared_ptr
			boost::shared_ptr<TypeNewBoundary> ptr(new TypeNewBoundary);
					
			// Assign the boundary id		
			ptr->id = id;
			
			// Assign the variable vector
			ptr->variables = var;

			// Store the pointer
			bc_ptrs.push_back(ptr);

			// Return the pointer
			return ptr;	
		}
		 	
		//! Adds a boundary object (string vector input)
		/*! 
		 * A template class for adding a boundary condition
		 * 
		 * \tparam TypeNewBoundary The type of boundary being added. The
		 * type added must be derived from EqBoundaryBase. See example1.cpp
		 * for an example.
		 * 
		 * \param id The integer identification for the boundary
		 * \param str A vector of strings indentifiy which variables to apply the boundary condition
		 *
		 * \see EqBoundaryBase
		 * \see example1.cpp
		 */ 
		template<class TypeNewBoundary> boost::shared_ptr<TypeNewBoundary> add_boundary(int id, vector<string> str){
			
			// Storage location for the variable indices
			vector<unsigned int> var;
			
			// Loop through each of the variables specified
			for (unsigned int i = 0; i < str.size(); i++){
				
				// Check that the variable exists
				if (!this->has_variable(str[i])){
					printf("ERROR: The variable, %s, does not exist (see ImplicitSystemBase::add_boundary(int id, vector<string> var)\n", str[i].c_str());
					exit(1);
				}
				
				// Add the variable index to the storage vector
				var.push_back(this->variable_number(str[i]));
			}
			
			// Call the index vector version of add_boundary
			add_boundary(id, var);
		}		 	

	protected:
				
		//! libMesh assembly function
		virtual void assemble() = 0;
		
		//! A shared pointer to the InitFuncBase class
		boost::shared_ptr<InitFuncBase<Number> > init_ptr;
				
		//! Storage vector for added boundaries
		vector< boost::shared_ptr<TypeBoundaryBase> > bc_ptrs;

		//! Apply the libMesh based Dirichlet boundary condition
		void apply_dirichlet(){
			
			// Loop through each of the boundaries
			for(unsigned int i = 0; i < bc_ptrs.size(); ++i){

				// Test that the boundary is a dirichlet boundary
				if(bc_ptrs[i]->type.compare("dirichlet") == 0){
					
					// Create a storage container for the boundary id
					std::set<boundary_id_type> boundary_id;
					boundary_id.insert(bc_ptrs[i]->id);
					
					// Apply the boundary condition using the function pointer (libmesh method)
					if(bc_ptrs[i]->fptr != NULL){
						// Create an AnalyticFunctoin from the supplied pointer
						AnalyticFunction<Number> func_object(bc_ptrs[i]->fptr);
						
						// Create the Boundary condition object
						DirichletBoundary bc(boundary_id, bc_ptrs[i]->variables, &func_object);
						
						// Apply the boundary condition
						this->get_dof_map().add_dirichlet_boundary(bc);

					// Apply the boundary condition with the value member (my method)
					} else {
						
						// Create a boost::function pointer
						boost::function< void (DenseVector<Number>&, const Point&, const Real)> bst_fptr;
						
						// Bind the pointer with the value member
						bst_fptr = boost::bind(&TypeBoundaryBase::value, bc_ptrs[i], _1, _2, _3);
						
						// Create an MyAnalyticFunction using the binded pointer
						MyAnalyticFunction<Number> func_object(bst_fptr);
						
						// Create the Boundary condition object
						DirichletBoundary bc(boundary_id, bc_ptrs[i]->variables, &func_object);
						
						// Apply the boundary condition
						this->get_dof_map().add_dirichlet_boundary(bc);
					}
				}
			}
		}

		 //! A method for getting a vector of indices for the boundary id
		 /*!
		  * This is a useful function for finding a boundary object
		  * with a specific id or to test if an id has been used.
		  * 
		  * \return A vector containing the index value of the given id.
		  * There should be only one value to this vector, but if empty 
		  * the boundary doesn't exist. The id's are checked for 
		  * uniqueness when they are added.
		  * 
		  * \param id The boundary id to search
		  * 
		  * \see HeatEq
		  */ 
		vector<int> get_boundary_index(const int id){
			
			// Create storage vector
			vector<int> idx;

			// Loop through each boundary
			for(int i = 0; i < (int)bc_ptrs.size(); ++i){
				
				// Extract the id for the current boundary
				int tmp = bc_ptrs[i]->id; 
				
				// Test the boundary and append the vector if they match
				if(tmp == id){
					idx.push_back(i);
				}
			}

			return idx; // Return the vector
		}

}; // end class
} //end namespace
#endif
