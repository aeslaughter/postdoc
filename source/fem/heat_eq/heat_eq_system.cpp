// Header file for this source code
#include "fem/heat_eq/heat_eq_system.h"
using namespace SlaughterFEM;
using namespace libMesh;

// The class constructor
HeatEqSystem :: HeatEqSystem(EquationSystems& es, const string& name, const unsigned int number) : 
			ImplicitSystemBase(es, name, number){}
			
// The class descturcor
HeatEqSystem :: ~HeatEqSystem(){}

// Stiffness and RHS assembly
void HeatEqSystem :: assemble(){

	// Get a constant reference to the mesh object.
	const MeshBase& mesh = get_mesh();

	// The dimension that we are running
	const unsigned int dim = mesh.mesh_dimension();

	// Get a constant reference to the Finite Element type
	// for the first (and only) variable in the system.
	FEType fe_type = variable_type(0);

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
	const std::vector<std::vector<Real> >& psi = fe_face->get_phi();

	// The element shape function gradients evaluated at the quadrature
	// points.
	const std::vector<std::vector<RealGradient> >& B = fe->get_dphi();

	// The XY locations of the quadrature points used for face integration
	const std::vector<Point>& qface_points = fe_face->get_xyz();

	// A reference to the \p DofMap object for this system.
	const DofMap& dof_map = get_dof_map();

	// Define data structures to contain the element matrix
	// and right-hand-side vector contribution.  
	DenseMatrix<Number> Me;			// element mass matrix
	DenseMatrix<Number> Ke;			// element stiffness matrix
	DenseVector<Number> Fe;			// element force vector (current time)
	DenseVector<Number> Fe_old;		// element force vector (previous time)
	DenseVector<Number> u_old;		// element temperature vector (previous time)
	
	DenseMatrix<Number> K_hat;		// general time integration stiffness matrix
	DenseVector<Number> F_hat;		// general time integration force vector

	// Storage for the degree of freedom indices
	std::vector<unsigned int> dof_indices;

	// Here we extract the parameters that will be needed
	const Real dt = get_constant<Real> ("dt");		// time step
	Real time = this->time;					 		// current time
	const Real k = get_constant<Real> ("k");		// conductivity
	const Real rho = get_constant<Real> ("rho");	// density
	const Real cp = get_constant<Real> ("cp");		// specific heat
	const Real theta = get_constant<Real>("theta");	// intergration parameter

	// Loop over all the elements in the mesh that are on local processor
	MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
	const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end(); 
	
	for ( ; el != end_el; ++el){    

		// Pointer to the element current element
		const Elem* elem = *el;

		// Get the degree of freedom indices for the current element
		dof_map.dof_indices(elem, dof_indices);

		// Compute the element-specific data for the current element
		fe->reinit (elem);

		// Zero the element matrices and vectors
		Me.resize (dof_indices.size(), dof_indices.size());
		Ke.resize (dof_indices.size(), dof_indices.size());
		Fe.resize (dof_indices.size());
		Fe_old.resize (dof_indices.size());

		// Compute the RHS and mass and stiffness matrix for this element (Me)
		for (unsigned int qp = 0; qp < qrule.n_points(); qp++){
			for (unsigned int i = 0; i < N.size(); i++){
				
				// This is where the source term will go in future
				Fe(i) += 0; 
				Fe_old(i) += 0; 

				for (unsigned int j = 0; j < N.size(); j++){				
					Me(i,j) += JxW[qp]*(N[i][qp]*rho*cp*N[j][qp]);  // mass matrix
					Ke(i,j) += JxW[qp]*(B[i][qp]*k*B[j][qp]);		// stiffness
				}
			}
		}

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
			u_old(i) = this->old_solution(dof_indices[i]);
		}
		
		// Build K_hat and F_hat (appends existing)
		K_hat.resize(dof_indices.size(), dof_indices.size());
		F_hat.resize(dof_indices.size());		
		build_stiffness_and_rhs(K_hat, F_hat, Me, Ke, Fe_old, Fe, u_old, dt, theta);
		
		// Applies the dirichlet constraints to K_hat and F_hat
		dof_map.heterogenously_constrain_element_matrix_and_vector(K_hat, F_hat, dof_indices);

		// Apply the local components to the global K and F
		this->matrix->add_matrix(K_hat, dof_indices);
		this->rhs->add_vector(F_hat, dof_indices);	

	} // (end) for ( ; el != end_el; ++el)
} // (end) assemble()
	
// A function for applying the general time integration scheme
void HeatEqSystem :: build_stiffness_and_rhs(DenseMatrix<Number>& K_hat, 
			DenseVector<Number>& F_hat, 
			DenseMatrix<Number> Me, 
			DenseMatrix<Number> Ke, 
			DenseVector<Number> Fe_old, 
			DenseVector<Number> Fe, 
			DenseVector<Number> u_old, 
			Real dt, 
			Real theta){

	// Build the stiffness matrix
	K_hat.add(1, Me);
	K_hat.add(dt*theta, Ke);

	// Create temporary storage for the RHS
	DenseVector<Number> a;
	DenseMatrix<Number> B;	
	DenseVector<Number> c;
	int m = K_hat.m();
	int n = K_hat.n();
	a.resize(m);
	B.resize (m, n);
	c.resize(m);
	
	// Generate the a-vector
	a.add(dt*(1-theta), Fe_old); 	// a = dt (1 - \theta) f_t
	a.add(dt*theta, Fe);		 	// a = dt [(1 - \theta) * f_t + \theta f_{t+1}]
		
	// Generate the B-matrix
	B.add(1, Me);					// B = M_e
	B.add((-1)*dt*(1-theta), Ke);   // B = [M_e - dt (1 - \theta) * Ke]
		
	// Build F_hat	
	B.vector_mult(c, u_old);		// c = [M_e - dt (1 - \theta) * Ke] * u_t
	F_hat.add(1, a);				// \hat{f} = a
	F_hat.add(1, c);				// \hat{f} = a + c
}

