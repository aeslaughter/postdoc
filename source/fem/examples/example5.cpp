/** \example fem/examples/example5.cpp 
 * A test function for testing the HeatEq class.
 * 
 * This class implements Example 8-1 from Bhatti (2005; p. 552). Run
 * the program with the \c --help flag for a complete list of run-time
 * options.
 */
 
// Standard library includes  
#include <iostream>
#include <algorithm>
#include <math.h>
#include <stdio.h>
#include <functional>

// BOOST includes
#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp>

// Include my fem library
#include "fem/include.h"
using namespace SlaughterFEM;

// Include my common library
#include "common/include.h"
using namespace SlaughterCommon;

// libMesh
#include <libmesh.h>
#include <libmesh_common.h>

#include <exodusII_io.h>
#include <vtk_io.h>

#include <mesh.h>
#include <mesh_generation.h>

#include <equation_systems.h>
#include <transient_system.h>
#include <nonlinear_implicit_system.h>
#include <explicit_system.h>

#include <fe.h>
#include <fe_type.h>
#include <quadrature_gauss.h>
#include <dof_map.h>
#include <sparse_matrix.h>
#include <numeric_vector.h>
#include <dense_matrix.h>
#include <dense_vector.h>

#include <point_locator_tree.h>
  
// Bring in everything from the libMesh namespace
using namespace libMesh;

// A class for containing the nodal data
class EqVolAvgData : public EqDataBase{
	
	public:
	
		// Class constructor
		EqVolAvgData(EquationSystems& sys) : EqDataBase(sys, "eq_data") : velocity_linker(sys){
			add_variable("epsilon");		// volume fraction
			//add_variable("tau_1");			// first stabilization term
			//add_variable("element_length"); // element length, Eq. 69, h
			//add_variable("P_alpha");		// stablization term, Eq. 75
			
			get_system().init();			// initialize the data storage system
		}
			
		// This is called when the solution is updated
		void value(DenseVector<Number>& output, const Point& p, const Real){
			
			// Gather the time from the system
			Real t = get_system().time;
			
			// Return the values of k and cp
			output.resize(2);
			output(0) = epsilon(p);
			output(1) = tau_1(p);
			output(2) = element_length(p);
		}		

		


		// A function for returning the velocity, from the momentum equation
		DenseVector<Number> velocity(const Point& p){
			
			// Get a reference to the system object for the Momentum equation
			TransientNonlinearImplicitSystem& momentum_system = 
				eq_sys.get_system<TransientNonlinearImplicitSystem>("momentum");
			
			// Get the number of dimensions
			const unsigned int dim = dimension();	
								
			// Collect the velocity vector components
			DenseVector<Number> s(dim);
			for (int d = 0; d < dim; d++){
				s(d) = momentum_system.point_value(d, p);
			}
			
			// Return the vector
			return s;
		}
		
		vector<Number> velocity(const Point& p){
			
			DenseVector<Number> a = velocity(p);
			vector<Number> b;
			for (unsigned int i = 0; i < a.size(); i++){
				b.push_back(a(i));
			}
			
			return b;
		}
		
		// A function for returning the volume fraction
		Number epsilon(const Point& p){
			return 1;
		}

		// Returns the \tau_1 value for the advective stabilization term (Eq. 63)
		Number tau_1(const Point& p){
			
			// Indicate the use of pow and min function from std library
			using std::pow;
			using std::min;

			// Collect the constants
			Number Da = eq_sys.parameters.get<Number>("Da");
			Number Pr = eq_sys.parameters.get<Number>("Pr");			
			
			// Collect the data evaluated at current point
			Number eps = epsilon(p);
			Number h = element_length(p);
			DenseVector<Number> v = velocity(p);

			// Compute the first \tau value
			Number tau_1a = pow(eps, 2) / pow((1-eps), 2) * (Da / Pr);
			
			// Compute the \tau_{SUPG} value, Eq. 64
			Number v_norm = v.l2_norm();			// ||v^h||
			Number Re = v_norm * h / (2 * Pr);		// Re_v, Eq. 67
			
			// z(Re) of Eq. 70
			Number zRe;
			if (Re >= 0 && Re <= 3){
				zRe = Re/3;
			} else {
				zRe = 1;
			}			

			// The \tau_{SUPG} value
			Number tau_1b = eps * h / (2 * v_norm) * zRe;
			
			// Return the minimum of the two values computed
			return min(tau_1a, tau_1b);
		}
		
		// Returns the element length, Eq. 69
		Number element_length(Elem* elem){
				
			// Get the number of dimensions
			const unsigned int dim = dimension();	
				
			// Get the variable index of the element length
			unsigned int idx = get_system().variable_number("element_length");	
				
			// Get the FE type for the element length variable
			FEType fe_type = get_system().variable_type(idx);

			// Build a Finite Element object of the specified type. 
			AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));

			// A Gauss quadrature rule for numerical integration.
			QGauss qrule (dim, fe_type.default_quadrature_order());

			// Tell the finite element object to use our quadrature rule.
			fe->attach_quadrature_rule(&qrule);

			// The element shape functions evaluated at the quadrature points.
			const std::vector<std::vector<RealGradient> >& dN = fe->get_dphi();
						
			// Re-initialize the fe system for this element
			fe->reinit(elem);
			
			// Initialize the element length parameter, h
			Number h = 0;
			
			// Get a vector of the points for the quadrature rule
			vector<Point>& pvec = qrule.get_points();
			
			/* Compute the element length at the Gauss points, this 
			 * differs from Eq. 69 which uses nodes.
			 */ 
			 
			// Loop over Gauss points (nodal summation of Eq. 69)
			for (unsigned int qp = 0; qp < qrule.n_points(); qp++){
				
				// Get the velocity vector
				DenseVector<Number> s = velocity(pvec[qp]);
		
				// Compute the L2-norm of the velocity vector
				Number s_norm = s.l2_norm();
				
				// Loop over shape functions
				for (unsigned int i = 0; i < dN.size(); i++){
				
					// If zero velocity, then element length is set to 0
					if (s_norm == 0 ){
						h = 0;
					
					// Else compute the element length
					} else {
						
						// Dot product, using normalized velocity
						for (int d = 0; d < dim; d++){
							h += 2 * dN[i][qp](d) * s(d)/s_norm;
						}
					}
				}
			}
			
			// Return the value of the element length
			return h;	
		}

		// Returns the P_alpha stabilization term, Eq. 75
		Number P_alpha(const Point& p){
				
			// Get the number of dimensions
			const unsigned int dim = dimension();	
				
			// Get the variable index of the element length
			unsigned int idx = get_system().variable_number("element_length");	
				
			// Get the FE type for the element length variable
			FEType fe_type = get_system().variable_type(idx);

			// Build a Finite Element object of the specified type. 
			AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));

			// A Gauss quadrature rule for numerical integration.
			QGauss qrule (dim, fe_type.default_quadrature_order());

			// Tell the finite element object to use our quadrature rule.
			fe->attach_quadrature_rule(&qrule);

			// The element shape functions evaluated at the quadrature points.
			const std::vector<std::vector<RealGradient> >& dN = fe->get_dphi();
			
			// Get a reference to the mesh
			MeshBase& mesh = eq_sys.get_mesh();
			
			// Determine the element that contains the point
			PointLocatorTree locate(mesh);
			
			// Locate the element
			const Elem* elem = locate(p);
			
			// Re-initialize the fe system for this element
			fe->reinit(elem);
			
			// Initialize the element length parameter, h
			Number P = 0;
			
			// Get a vector of the points for the quadrature rule
			vector<Point>& pvec = qrule.get_points();
			
			// Loop over Gauss points (nodal summation of Eq. 69)
			for (unsigned int qp = 0; qp < qrule.n_points(); qp++){
				
				// Gather the necessary components
				Number eps = epsilon(pvec[qp]);
				Number tau = tau_1(pvec[qp]);
				DenseVector<Number> v = velocity(pvec[qp]);
				
				// Loop over shape functions
				for (unsigned int i = 0; i < dN.size(); i++){
				
					// Dot product, using normalized velocity
					for (int d = 0; d < dim; d++){
						P += (1 / eps) * tau * v(d) * dN[i][qp](d);
					}
					
				}
			}
			
			// Return the value P_{\alpha}^e
			return P;	
		}
	
	private:
		EqVariableLinker velocity_linker;

};

// Function for initializing Momentum Eq.
void momentum_init_velocity(DenseVector<Number>& output, const Point&, const Real){
	output(0) = 0;
	output(1) = 2;
}

// Function for initializing Energy Eq.
void energy_init_enthalpy(DenseVector<Number>& output, const Point&, const Real){
	output(0) = 0;
}

// Function prototype for assembly of Energy Eq.
void energy_assemble(EquationSystems& eq_sys, const std::string& system_name);


// Initialize Momentum Eq.
void init_function(EquationSystems& es, const std::string& system_name){
	
	// Get a reference to the system
	TransientNonlinearImplicitSystem& system =
	  es.get_system<TransientNonlinearImplicitSystem>(system_name);

	// Create a function acceptable for solution projection
	boost::function< void (DenseVector<Number>&, const Point&, const Real)> fptr;
	
	// Link the correct function
	if (system_name.compare("momentum") == 0){
		fptr = momentum_init_velocity;
	
	} else if (system_name.compare("energy") == 0){
		fptr = energy_init_enthalpy;
	}

	// Convert into a libmesh acceptable format
	MyAnalyticFunction<Number> fobj(fptr);

	// Project the solution 
	system.project_solution(&fobj);
}
        

// Begin main function
int main (int argc, char** argv){

    // Initialize libraries
    LibMeshInit init (argc, argv);

	// Generate a mesh object
	MyMesh mesh;
	
	// Create a 2D grid
	MeshTools::Generation::build_square(mesh, 1, 1, 0., 1, 0., 1, QUAD4);
	mesh.all_first_order();
	Order order = FIRST;

	// Create an equation system
 	EquationSystems eq_sys(mesh); 

	// Add constants
	eq_sys.parameters.set<Number>("Da") = 1;
	eq_sys.parameters.set<Number>("Pr") = 1;
	eq_sys.parameters.set<Number>("dt") = 0.01;

	// Create momentum equation
	TransientNonlinearImplicitSystem& momentum =
		eq_sys.add_system<TransientNonlinearImplicitSystem>("momentum");

	// Add 2D velocity variables
	momentum.add_variable("vx", order);
	momentum.add_variable("vy", order);
	
	// Attach initialization function
	momentum.attach_init_function(init_function);
	momentum.init();

	// Create the energy equation
	TransientNonlinearImplicitSystem& energy =
		eq_sys.add_system<TransientNonlinearImplicitSystem>("energy");
	
	// Add the enthalpy
	energy.add_variable("h", order);
	
	// Attach initialization function
	energy.attach_init_function(init_function);
	
	// Attach assembly function
	energy.attach_assemble_function(energy_assemble);

	// Initialize the energy
	energy.init();

	// Create the data class (Initializes upon creation)
	EqData data(eq_sys);
	
	// Project the solution at time t = 0
	data.update_solution(0);
	
	// Output the data
	ExodusII_IO(mesh).write_equation_systems("example5.ex2", eq_sys);

	// Call the assembly function for testing
	energy_assemble(eq_sys, "energy");



} // end main


// Energy equation assembly function
// Equation references are from Zabaras and Samanta, 2004
void energy_assemble(EquationSystems& eq_sys, const std::string& system_name){

	// Check that we are assembling the correct system
	libmesh_assert (system_name == "energy");

	// Get a constant reference to the mesh object.
	const MeshBase& mesh = eq_sys.get_mesh();

	// The dimension that we are running
	const unsigned int dim = mesh.mesh_dimension();

	// Get a reference to the system object for the Energy equation
	TransientNonlinearImplicitSystem& system = 
		eq_sys.get_system<TransientNonlinearImplicitSystem>(system_name);
		
	// Get a reference to the system object for the Momentum equation (velocity)
	TransientNonlinearImplicitSystem& momentum_system = 
		eq_sys.get_system<TransientNonlinearImplicitSystem>("momentum");
			
	// Get a constant reference to the Finite Element type
	// for the first (and only) variable in the system.
	FEType fe_type = system.variable_type(0);

	// Build a Finite Element object of the specified type. 
	AutoPtr<FEBase> fe      (FEBase::build(dim, fe_type));
	AutoPtr<FEBase> fe_face (FEBase::build(dim, fe_type));

	// A Gauss quadrature rule for numerical integration.
	// Let the \p FEType object decide what order rule is appropriate.
	QGauss qrule (dim,   fe_type.default_quadrature_order());
	QGauss qface (dim-1, fe_type.default_quadrature_order());

	// Tell the finite element object to use our quadrature rule.
	fe->attach_quadrature_rule      (&qrule);
	fe_face->attach_quadrature_rule (&qface);

	// Here we define some references to cell-specific data that
	// will be used to assemble the linear system.  We will start
	// with the element Jacobian * quadrature weight at each integration point.   
	const std::vector<Real>& JxW      = fe->get_JxW();
	const std::vector<Real>& JxW_face = fe_face->get_JxW();

	// The element shape functions evaluated at the quadrature points.
	const std::vector<std::vector<Real> >& N = fe->get_phi();
	const std::vector<std::vector<Real> >& N_face = fe_face->get_phi();

	// Element shape function gradients evaluated at quadrature points
	const std::vector<std::vector<RealGradient> >& dN = fe->get_dphi();

	// The XY locations of the quadrature points used for face integration
	const std::vector<Point>& qface_points = fe_face->get_xyz();

	// A reference to the \p DofMap object for this system.
	const DofMap& dof_map = system.get_dof_map();

	// Define data structures to contain the element matrix
	// and right-hand-side vector contribution (Eq. 107)  
	DenseMatrix<Number> Me;			// [\hat{M} + \hat{M}_{\delta}]
	DenseMatrix<Number> Ne;			// [\hat{N} + \hat{N}_{\delta}]
	DenseMatrix<Number> Ke;			// [\hat{K} + \hat{K}_{\delta}]
	DenseVector<Number> Fe;			// [\hat{F} + \hat{F}_{\delta}]
	//DenseVector<Number> Fe_old;		// element force vector (previous time)
	DenseVector<Number> h_old;		// element enthalpy vector (previous time)
	
	DenseVector<Number> v_old;		// element velocity vector (previous time)

	DenseMatrix<Number> Mstar;		// general time integration stiffness matrix (Eq. 125)
	DenseVector<Number> R;			// general time integration force vector (Eq. 126)

	// Storage for the degree of freedom indices
	std::vector<unsigned int> dof_indices;

	// Here we extract the parameters that will be needed
	const Real dt = eq_sys.parameters.get<Real>("dt");		// time step
	Real time = system.time;			// current time

	// Loop over all the elements in the mesh that are on local processor
	MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
	const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end(); 

	for ( ; el != end_el; ++el){    

		// Pointer to the element current element
		const Elem* elem = *el;

		// Get the degree of freedom indices for the current element
		dof_map.dof_indices(elem, dof_indices);

		// Compute the element-specific data for the current element
		fe->reinit(elem);

		// Zero the element matrices and vectors
		Me.resize (dof_indices.size(), dof_indices.size());
		Ne.resize (dof_indices.size(), dof_indices.size());
		Ke.resize (dof_indices.size(), dof_indices.size());
		Fe.resize (dof_indices.size());
		h_old.resize (dof_indices.size());

		// Compute the RHS and mass and stiffness matrix for this element (Me)
		for (unsigned int qp = 0; qp < qrule.n_points(); qp++){
			
			// Get the velocity at this Gauss point (std::vector format)
			vector v(3,0); // Need to get this from my data class
			
			for (unsigned int i = 0; i < N.size(); i++){
				
				// This is not correct
				Fe(i) += 0; 

				for (unsigned int j = 0; j < N.size(); j++){		

					Me(i,j) += JxW[qp]*(N[i][qp]) * N[j][qp]);  
					Ne(i,j) += JxW[qp]*(N[i][qp]*v[j]*dN[j][qp]);	
				}
			}
		}
/*   
	    // BOUNDARY CONDITIONS    
		// Loop through each side of the element for applying boundary conditions
		for (unsigned int s = 0; s < elem->n_sides(); s++){
			
			// Only consider the side if it does not have a neighbor
			if (elem->neighbor(s) == NULL){
				
				// Pointer to current element side
				const AutoPtr<Elem> side = elem->side(s);
					
				// Boundary ID of the current side
				int boundary_id = (mesh.boundary_info)->boundary_id(elem, s);

				// Get index of the boundary class with the same id
				// this vector is empty if there is no match and only
				// contains a single value if there is a match
				std::vector<int> idx = get_boundary_index(boundary_id);	
				
				// Continue of there is a match						
				if(!idx.empty()){							
							
					// Compute the shape function values on the element face
					fe_face->reinit(elem, s);
											
					// Create a shared pointer to the boundary class	
					boost::shared_ptr<HeatEqBoundaryBase> ptr = bc_ptrs[idx[0]];
					
					// Determine the type of boundary considered
					std::string type = ptr->type;

					// Loop through quadrature points
					for (unsigned int qp = 0; qp < qface.n_points(); qp++){

						// DIRICHLET (libMesh version; handled at initialization)
						if(type.compare("dirichlet") == 0){
							// The dirichlet conditions are handled at initlization
							// but I don't want to throw an error if they are
							// encountered, so just do nothing
							
						// NEUMANN condition
						} else if(type.compare("neumann") == 0){

							// Current and past flux values					
							const Number q = ptr->q(qface_points[qp], time);
							const Number q_old = ptr->q(qface_points[qp], time - dt);
									
							// Add values to Fe						
							for (unsigned int i = 0; i < psi.size(); i++){
								Fe(i) 	  += JxW_face[qp] * q * psi[i][qp];
								Fe_old(i) += JxW_face[qp] * q_old * psi[i][qp];		
							}	

						// CONVECTION boundary
						} else if(type.compare("convection") == 0){	

							// Current and past h and T_inf
							const Number h 		  = ptr->h(qface_points[qp], time);
							const Number h_old    = ptr->h(qface_points[qp], time - dt);
							const Number Tinf	  = ptr->Tinf(qface_points[qp], time);
							const Number Tinf_old = ptr->Tinf(qface_points[qp], time - dt);
							
							// Add values to Ke and Fe						
							for (unsigned int i = 0; i < psi.size(); i++){
								Fe(i) 	  += (1) * JxW_face[qp] * h * Tinf * psi[i][qp];
								Fe_old(i) += (1) * JxW_face[qp] * h_old * Tinf_old * psi[i][qp];
								
								for (unsigned int j = 0; j < psi.size(); j++){
									Ke(i,j) += JxW_face[qp] * psi[i][qp] * h * psi[j][qp];

								}		
							}
				
						// Un-registerd type		
						} else {
							printf("WARNING! The boundary type, %s, was not understood!\n", type.c_str());
					
						}	// (end) type.compare(...) statemenst	
					} //(end) for (int qp = 0; qp < qface.n_points(); qp++)	
				} // (end) if(!idx.empty)
			} // (end) if (elem->neighbor(s) == NULL){
		} // (end) for (int s = 0; s < elem->n_sides(); s++)	

		// Zero the pervious time-step temperature vector for this element
		u_old.resize(dof_indices.size());
		
		// Gather the temperatures at the nodes
		for (unsigned int i = 0; i < psi.size(); i++){
			u_old(i) = system.old_solution(dof_indices[i]);
		}

		// Build K_hat and F_hat (appends existing)
		K_hat.resize(dof_indices.size(), dof_indices.size());
		F_hat.resize(dof_indices.size());
		build_stiffness_and_rhs(K_hat, F_hat, Me, Ke, Fe_old, Fe, u_old, dt, theta);
		
		// Applies the dirichlet constraints to K_hat and F_hat
		dof_map.heterogenously_constrain_element_matrix_and_vector(K_hat, F_hat, dof_indices);
		
		// Apply the local components to the global K and F
		system.matrix->add_matrix(K_hat, dof_indices);
		system.rhs->add_vector(F_hat, dof_indices);	
*/		
	} // (end) for ( ; el != end_el; ++el)

} // (end) assemble()




