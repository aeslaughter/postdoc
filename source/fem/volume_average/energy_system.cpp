 /*! \file energy_system.cpp
  * \brief Source code for the volume averaging EnergySystem class.
  */
  
#include "fem/volume_average/energy_system.h"
using namespace SlaughterFEM;

//! [constructor]
// Energy Equation class constructor
EnergySystem :: EnergySystem(EquationSystems& es, const string& name, const unsigned int number) :
	SystemBase(es, name, number){
	
	// Indidicate that the enthalpy initilization method is NOT used,
	// the default assumes that the initilization function is temperature
	enthalpy_init_ = false;

	// Add the unknown enthalpy, h, and enthalpy rate, h_dot
	this->add_variable("h", FIRST);		
	//this->add_variable("h_dot", order);	
	//this->add_variable("delta_h_dot", order);	

	set_constant<Number>("gamma", 0.5);

	// Add matrices
	this->add_matrix("M");
	this->add_matrix("N");
	this->add_matrix("K");
	this->add_vector("F");
	
	// Add solution matrices
	this->add_vector("h_dot");
	this->add_vector("delta_h_dot");
	
} 
//! [constructor]

// Eq. 12 (used for initialization of EnergySystem)
void EnergySystem :: temperature_to_enthalpy(DenseVector<Number>& output, const Point& p, const Real t){

	// ThermoEq pointer must be defined 
	if (!thermo){
		printf("ERROR: ThermoEq pointer must be defined.\n");
		libmesh_error();
	}

	// Temperature for the current point via the intialization function
	Number T;
	
	// Get T using boost::function
	if(bst_init_func != NULL){
		DenseVector<Number> x(1);
		bst_init_func(x, p, t);
		T = x(0);

	// Get T using a function pointer directly
	} else if (init_func != NULL){
		T = init_func(p, this->get_equation_systems().parameters, this->name(), NULL);
		
	// Get T using EqInitBase class (first component value)
	} else if (init_ptr){
		T = init_ptr->component(0, p, t);
		
	} else {
		printf("ERROR: No initilization function exists\n");
		libmesh_error();
	}
	
	// Compute necessary terms
	Number Tliq = thermo->T_liq(p);
	Number Tsol = thermo->T_sol(p);
	Number hsol = thermo->h_sol(p);	
	
	// Get the necessary constants
	Number Te = thermo->get_constant<Number>("eutectic_temperature");
	Number hf = thermo->get_constant<Number>("latent_heat");
	Number cf = thermo->get_constant<Number>("specific_heat_fluid");
	Number cs = thermo->get_constant<Number>("specific_heat_solid");
	
	// Compute the correct value for enthalpy
	if (T > Tliq){
		Number f = 1;
		output(0) = f * ((cf-cs)*(T-Te) + hf) + cs*T;
	
	} else if (Te < T && T <= Tliq) {
		Number f = thermo->lever_rule(p,T);
		output(0) = f * ((cf-cs)*(T-Te) + hf) + cs*T;
	
	} else if (Tsol < T && T <= Te){
		Number B = ((cf-cs)*(T-Te) + hf);
		output(0) = (B*hsol - hf*cs*T) / (B - hf);
		
	} else if (T <= Tsol){
		output(0) = T*cs;
	}
}
		
// Compute the alpha term of Eq. 69
Number EnergySystem :: alpha(Gradient grad_T, Gradient grad_h, Number f){
	
	Real mag_grad_h = grad_h.size();
		
	if (mag_grad_h > 0){
		return (grad_T * grad_h) / (mag_grad_h * mag_grad_h);
	
	} else {
		const Number cf = thermo->get_constant<Number>("specific_heat_fluid");
		const Number cs = thermo->get_constant<Number>("specific_heat_solid");
		return 1 / (cs + cf*f - cs*f);  
		//return - 1 / ((cf - cs) * (f - 1));
	}
} 

void EnergySystem :: solve(){
printf("here1\n");	
	assemble();
	update_rhs();
	this->matrix->close();
	this->update();
	
	PetscLinearSolver<Number> solver;
printf("here2\n");		
	const double tol = 0.01;
	const unsigned int iter = 100;
printf("here3\n");		
	this->matrix->print(std::cout);
	this->rhs->print(std::cout);
printf("here4\n");	
	
	solver.solve(*this->matrix, this->get_vector("delta_h_dot"), *this->rhs, tol, iter);

	
	
}

//! [initialize]
// Initilize function
void EnergySystem :: initialize(){

	// Initialization with enthalpy function directly
	if(enthalpy_init_){
		SystemBase :: initialize();
	
	// Initialize the enthalpy from	a temperature function (default)
	} else {
		// Create boost::function pointers
		boost::function< void (DenseVector<Number>&, const Point&, const Real)> fptr_shrt;

		// Bind pointer to this class
		fptr_shrt = boost::bind(&EnergySystem::temperature_to_enthalpy, this, _1, _2, _3);

		// Project the solution
		MyAnalyticFunction<Number> fobj(fptr_shrt);	
		this->project_solution(&fobj);
	}
	

	// Set h_dot to zero
	this->get_vector("h_dot").zero();
	
	
	Number dt = get_constant<Number>("dt");
	this->update_solution(0,dt);
	
	// Set initialization flag to true
	_initialized = true;
} 
//! [initialize]
		
void EnergySystem :: update_rhs(){
	
	this->rhs->add(this->get_vector("F"));

		//	R.add(1,Fe);

		//DenseVector<Number> a(dof_indices.size());
		//Me.vector_mult(a, h_dot);
		//R.add(-1, a);

	//	DenseMatrix<Number> B(dof_indices.size(), dof_indices.size());
	//	DenseVector<Number> b(dof_indices.size());
	//	B.add(1,Ne);
	//	B.add(1,Ke);

	//	B.vector_mult(b, h);

	//	R.add(-1,b);

	//	system.rhs->add_vector(R, dof_indices);
	
	
}		
		
		
// Stiffness and RHS assembly
// Equation references are from Samanta and Zabaras, 2005
void EnergySystem :: assemble(){

// GENERAL VARIABLES
	// Get a constant reference to the mesh object.
	const MeshBase& mesh = this->get_mesh();

	// The dimension that we are running
	const unsigned int dim = this->ndim();

// FEM THERMODYNAMIC RELATIONSHIPS (ThermoEq Class) 
	// Determine the FEM type (should be same for all ThermoEq variables)
	FEType fe_type_thermo = thermo->variable_type(0);
	
	// Build FE object; accessed via a pointer
	AutoPtr<FEBase> fe_thermo(FEBase::build(dim, fe_type_thermo));

	// Setup a quadrature rule
	QGauss qrule_thermo(dim, fe_type_thermo.default_quadrature_order());
	
	// Link FE and Quadrature
	fe_thermo->attach_quadrature_rule(&qrule_thermo);

	// References to shape functions and derivatives
	const vector<std::vector<Real> >& N_thermo = fe_thermo->get_phi();
	const vector<std::vector<RealGradient> >& B_thermo = fe_thermo->get_dphi();
		
	// Setup a DOF map
	const DofMap& dof_map_thermo = thermo->get_dof_map(); 

// FEM MOMENTUM EQUATION
	// Determine the FEM type 
	FEType fe_type_momentum = momentum->variable_type(0);
	
	// Build FE object; accessed via a pointer
	AutoPtr<FEBase> fe_momentum(FEBase::build(dim, fe_type_momentum));

	// Setup a quadrature rule
	QGauss qrule_momentum(dim, fe_type_momentum.default_quadrature_order());
	
	// Link FE and Quadrature
	fe_momentum->attach_quadrature_rule(&qrule_momentum);

	// References to shape functions and derivatives
	const vector<std::vector<Real> >& N_momentum = fe_momentum->get_phi();
		
	// Setup a DOF map
	const DofMap& dof_map_momentum = momentum->get_dof_map(); 
	
// FEM ENERGY EQ. RELATIONSHIPS
	// Get a constant reference to the Finite Element type
	// for the first (and only) variable in the system.
	FEType fe_type = this->variable_type(0);

	// Build a Finite Element object of the specified type
	AutoPtr<FEBase> fe      (FEBase::build(dim, fe_type));
	AutoPtr<FEBase> fe_face (FEBase::build(dim, fe_type));

	// A Gauss quadrature rule for numerical integration.
	// Let the \p FEType object decide what order rule is appropriate.
	QGauss qrule (dim,   fe_type.default_quadrature_order());
	QGauss qface (dim-1, fe_type.default_quadrature_order());

	// Tell the finite element object to use our quadrature rule.
	fe->attach_quadrature_rule(&qrule);
	fe_face->attach_quadrature_rule(&qface);

	// Here we define some references to cell-specific data that
	// will be used to assemble the linear system.  We will start
	// with the element Jacobian * quadrature weight at each integration point.   
	const vector<Real>& JxW      = fe->get_JxW();
	const vector<Real>& JxW_face = fe_face->get_JxW();

	// The element shape functions evaluated at the quadrature points.
	const vector<std::vector<Real> >& N = fe->get_phi();
	const vector<std::vector<Real> >& N_face = fe_face->get_phi();

	// Element shape function gradients evaluated at quadrature points
	const vector<std::vector<RealGradient> >& B = fe->get_dphi();

	// The XY locations of the quadrature points used for face integration
	const vector<Point>& qface_points = fe_face->get_xyz();

	// A reference to the \p DofMap objects
	const DofMap& dof_map = this->get_dof_map(); // this system

// DEFINE VECTOR AND MATRIX VARIABLES
	// Define data structures to contain the element matrix
	// and right-hand-side vector contribution (Eq. 107)  
	DenseMatrix<Number> Me;			// [\hat{M} + \hat{M}_{\delta}]
	DenseMatrix<Number> Ne;			// [\hat{N} + \hat{N}_{\delta}]
	DenseMatrix<Number> Ke;			// [\hat{K} + \hat{K}_{\delta}]
	DenseVector<Number> Fe;			// [\hat{F} + \hat{F}_{\delta}]
	
	//DenseVector<Number> Fe_old;		// element force vector (previous time)
	DenseVector<Number> h;			// element enthalpy vector (previous time)
	DenseVector<Number> h_dot;
	//DenseVector<Number> delta_h_dot;
	
	DenseMatrix<Number> Mstar;		// general time integration stiffness matrix (Eq. 125)
	DenseVector<Number> R;			// general time integration force vector (Eq. 126)

	// Storage vectors for the degree of freedom indices 
	std::vector<unsigned int> dof_indices;		// this system (h)
//	std::vector<unsigned int> dof_indices_hdot;		
//	std::vector<unsigned int> dof_indices_deltahdot;		
	std::vector<unsigned int> dof_indices_velocity;	// this system
	std::vector<unsigned int> dof_indices_rho;	// ThermoEq density
	std::vector<unsigned int> dof_indices_tmp;  // ThermoEq temperature
	std::vector<unsigned int> dof_indices_f;  	// ThermoEq liquid fraction
	std::vector<unsigned int> dof_indices_eps;  // ThermoEq epsilon
	
	// Define the necessary constants
	const Number gamma = get_constant<Number>("gamma");
	const Number dt = get_constant<Number>("dt");		// time step
	Real time = this->time;					// current time
	const Number ks = thermo->get_constant<Number>("conductivity_solid");
	const Number kf = thermo->get_constant<Number>("conductivity_fluid");
	const Number cs = thermo->get_constant<Number>("specific_heat_solid");
	const Number cf = thermo->get_constant<Number>("specific_heat_fluid");
	const Number Te = thermo->get_constant<Number>("eutectic_temperature");
	const Number hf = thermo->get_constant<Number>("latent_heat");

	// Index of density variable in ThermoEq system
	const unsigned int rho_idx = thermo->variable_number("density");
	const unsigned int tmp_idx = thermo->variable_number("temperature");
	const unsigned int f_idx   = thermo->variable_number("liquid_mass_fraction");
	const unsigned int eps_idx = thermo->variable_number("epsilon");

	// Loop over all the elements in the mesh that are on local processor
	MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
	const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end(); 

	for ( ; el != end_el; ++el){    

		// Pointer to the element current element
		const Elem* elem = *el;
		
		// Get the degree of freedom indices for the current element
		dof_map.dof_indices(elem, dof_indices, 0);
		//dof_map.dof_indices(elem, dof_indices_hdot, 1);
		//dof_map.dof_indices(elem, dof_indices_deltahdot, 2);
		dof_map_momentum.dof_indices(elem, dof_indices_velocity);		
		dof_map_thermo.dof_indices(elem, dof_indices_rho, rho_idx);
		dof_map_thermo.dof_indices(elem, dof_indices_tmp, tmp_idx);
		dof_map_thermo.dof_indices(elem, dof_indices_f, f_idx);
		dof_map_thermo.dof_indices(elem, dof_indices_eps, eps_idx);

		// Compute the element-specific data for the current element
		fe->reinit (elem);
		fe_thermo->reinit(elem);
		fe_momentum->reinit(elem);

		// Zero the element matrices and vectors
		Me.resize (dof_indices.size(), dof_indices.size());		// [\hat{M} + \hat{M}_{\delta}]
		Ne.resize (dof_indices.size(), dof_indices.size());		// 
		Ke.resize (dof_indices.size(), dof_indices.size());
		Fe.resize (dof_indices.size());
		
		// Extract a vector of quadrature x,y,z coordinates
		const vector<Point> qp_vec = fe->get_xyz(); 

		// Compute the element length, h
		Number elem_length = thermo->element_length(elem);

		// Compute the RHS and mass and stiffness matrix for this element (Me)
		for (unsigned int qp = 0; qp < qrule.n_points(); qp++){
				
			// Get the velocity vector at this point (old value)
			VectorValue<Number> v;
			for (unsigned int i = 0; i < N_momentum.size(); i++){
				for (unsigned int j = 0; j < dim; j++){
					v(j) += N_momentum[i][qp] * momentum->old_solution(dof_indices_velocity[2*i+j]);
				}
			}
		
			// Compute ThermoEq variables; must be mapped from node to quadrature points
			Number T = 0;
			Gradient grad_T;
			Number f = 0;
			Gradient grad_f;			
			Number rho = 0;
			Number rho_old = 0;
			Number eps = 0;

			for (unsigned int i = 0; i < N_thermo.size(); i++){
				T += N_thermo[i][qp] * thermo->current_solution(dof_indices_tmp[i]);
				grad_T.add_scaled(B_thermo[i][qp], thermo->current_solution(dof_indices_tmp[i]));
				f += N_thermo[i][qp] * thermo->current_solution(dof_indices_f[i]);
				grad_f.add_scaled(B_thermo[i][qp], thermo->current_solution(dof_indices_f[i]));
				rho += N_thermo[i][qp] * thermo->current_solution(dof_indices_rho[i]);
				rho_old += N_thermo[i][qp] * thermo->old_solution(dof_indices_rho[i]);
				eps += N_thermo[i][qp] * thermo->current_solution(dof_indices_eps[i]);
			}

			// Compute EnergySystem variables
			Gradient grad_h;
			for (unsigned int i = 0; i < B.size(); i++){
				grad_h.add_scaled(B[i][qp], this->current_solution(dof_indices[i]));
			}

			// Compute T_{,k}^h v_k^h and f_{,k} v_k^h summation terms for F
			Number Tv = 0;
			Number fv = 0;
			for (unsigned int i = 0; i < dim; i++){
				Tv += grad_T(i) * v(i);
				fv += grad_f(i) * v(i);
			}

			// Compute the time derivative of density
			const Number drho_dt = (rho - rho_old)/dt;

			// Compute alpha term of Eq. 69
			const Number alpha = this->alpha(grad_T, grad_h, f);
			
			// Extract tau_1 stabilization term		
			const Number tau_1 = thermo->tau_1(qp_vec[qp], elem_length);
					
			// Loop through the components and construct matrices
			for (unsigned int i = 0; i < N.size(); i++){
				
				// Compute advective stabilization term (Eq. A, p. 1777)
				const Number d = tau_1 * v * B[i][qp] / f - tau_1 * 1/rho * drho_dt * (1-f)/f * N[i][qp];										
				// Force vector, Eq. 77
				Number F1 = JxW[qp] * (N[i][qp] + d) * rho * (1 - f) * (cf - cs) * Tv;   	
				Number F2 = JxW[qp] * (N[i][qp] + d) * rho * fv * ((cf - cs) * (T - Te) + hf);	
				Number F3 = JxW[qp] * (N[i][qp] + d) * drho_dt * (1 - f) * ((cf - cs) * (T - Te) + hf);	
				Fe(i) += F1 + F2 + F3;

				// Build the stiffness matrices
				for (unsigned int j = 0; j < N.size(); j++){
		
					// Mass matrix, Eq. 108
					Me(i,j) += JxW[qp] * rho * ((N[i][qp] + d) * N[j][qp]); 
					
					// Stiffness matrix one, Ne, Eq. 109
					Ne(i,j) += JxW[qp] * rho * ((N[i][qp] + d) * (v * B[j][qp]));
					
					// Stiffness matrix two, Ke, Eq. 110
					Ke(i,j) += JxW[qp]*((eps*kf + (1 - eps)*ks) * alpha * B[i][qp] * B[j][qp]);
				}
			}
		}
		
		printf("Me:\n");
		Me.print(std::cout);
		
		printf("\nNe:\n");
		Ne.print(std::cout);
		
		printf("\nKe:\n");
		Ke.print(std::cout);
		
		printf("\nFe:\n");
		Fe.print(std::cout);		
	
	
		h.resize(dof_indices.size());
		h_dot.resize(dof_indices.size());
	//	delta_h_dot.resize(dof_indices_deltahdot.size());


		for (unsigned int i = 0; i < dof_indices.size(); i++){
			h(i) = this->old_solution(dof_indices[i]);
			h_dot(i) = this->get_vector("h_dot")(dof_indices[i]);
	//		delta_h_dot(i) = this->old_solution(dof_indices_deltahdot[i]);
		}
		
		this->get_matrix("M").add_matrix(Me, dof_indices);
		this->get_matrix("N").add_matrix(Ne, dof_indices);
		this->get_matrix("K").add_matrix(Ke, dof_indices);
		this->get_vector("F").add_vector(Fe, dof_indices);
		
		Mstar.resize(dof_indices.size(), dof_indices.size());	
		R.resize(dof_indices.size());

		// Me + gamma*dt*(Ke + Ne);
		Mstar.add(1,Me);		
		Mstar.add(gamma*dt,Ke);
		Mstar.add(gamma*dt,Ne);

		this->matrix->add_matrix(Mstar, dof_indices);

		R.add(1,Fe);

		DenseVector<Number> a(dof_indices.size());
		Me.vector_mult(a, h_dot);
		R.add(-1, a);

		DenseMatrix<Number> B(dof_indices.size(), dof_indices.size());
		DenseVector<Number> b(dof_indices.size());
		B.add(1,Ne);
		B.add(1,Ke);

		B.vector_mult(b, h);

		R.add(-1,b);

		this->rhs->add_vector(R, dof_indices);

	
		
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
*/		
	} // (end) for ( ; el != end_el; ++el)

	//update_rhs();


} // (end) assemble()

