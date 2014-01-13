/*! \file level_set_system.cpp
  * \brief Source code for LevelSetSystem class.
  */
  
// Include the header file 
#include "fem/volume_average/level_set_system.h"
#include "fem/common/my_dense_matrix.h"
using namespace SlaughterFEM;

// Constructor
LevelSetSystem :: LevelSetSystem(EquationSystems& es, const string& name, const unsigned int number) : ImplicitSystemBase(es, name, number){
	
	
	this->add_variable("phi", FIRST, L2_HIERARCHIC);

	this->time = 0;	

	this->add_vector("_oldest_local_solution");
	//this->add_vector("_old_solution");
	//this->add_matrix("_stiffness_matrix");
	this->add_matrix("_mass_matrix_inverse");	
		
	_count = 0;
	
	this->h = 1;
	
}

LevelSetSystem :: ~LevelSetSystem(){
	this->clear();
	
	//old_local_solution.release();
	//older_local_solution.release();
	//oldest_local_solution.release();
	//temporary_local_solution.release();
	
}

void LevelSetSystem :: initialize(){
	
	if (velocity == NULL){
		printf("ERROR: The pointer to the FrontVelocityEq class must be defined.\n");
		libmesh_error();
	}
	
	ImplicitSystemBase :: initialize();
	
}

/*
Number LevelSetSystem :: oldest_solution(const unsigned int global_dof_number){

	// Check the sizes
	libmesh_assert (global_dof_number < this->get_dof_map().n_dofs());
	libmesh_assert (global_dof_number < oldest_local_solution->size());
  
	// Return the desired value
	return (*oldest_local_solution)(global_dof_number);
}


Number LevelSetSystem :: temporary_solution(const unsigned int global_dof_number){

	// Check the sizes
	libmesh_assert (global_dof_number < this->get_dof_map().n_dofs());
	libmesh_assert (global_dof_number < temporary_local_solution->size());
  
	// Return the desired value
	return (*temporary_local_solution)(global_dof_number);
}
*/

void LevelSetSystem :: update_solution(){
	update_solution(this->time, this->time_step());

}
	
void LevelSetSystem :: update_solution(Real time, Real dt){
	
	 // Apply the new time and time step
	 //this->time = time;
	
	 // Step the system with time
	// *this->oldest_local_solution = *this->older_local_solution;
	 *this->older_local_solution  = *this->old_local_solution;
	 *this->old_local_solution    = *this->current_local_solution;
	 
	 // Update the contraints (dirichlet boundaries)
	// const MeshBase& mesh = this->get_mesh();
	// DofMap& dof_map = this->get_dof_map();
	// dof_map.create_dof_constraints(mesh, time);
}

Number LevelSetSystem :: time_step(){
	Number p = 2; // using a second-order Runge-Kutta method
	Number vmax = this->velocity->system().solution->max();
	Number vmin = std::abs(this->velocity->system().solution->min());
	if (vmin > vmax){
		vmax = vmin;
	}
	printf("vmax = %g; h = %g\n", vmax, h);
	
	
	return 0.5*this->h / (vmax * (2*p + 1));
	//return 0.001;
}

void LevelSetSystem :: solve(){
	
	// Assemble the system
	if (this->assemble_before_solve){
		this->assemble();
	}

	// Get the time-step
	Number dt = this->time_step();
	printf("dt = %g\n", dt);
	this->time += dt;
	
	// Pointer to the mass matix inverse and stiffness matrices
	SparseMatrix<Number>* M = this->request_matrix("_mass_matrix_inverse");
	SparseMatrix<Number>* K = this->matrix;
	
	// Pointers to solution vectors
	AutoPtr<NumericVector<Number> > un = this->solution->clone();
	AutoPtr<NumericVector<Number> > u1 = this->solution->clone();	
	AutoPtr<NumericVector<Number> > Lun = this->solution->clone();
	AutoPtr<NumericVector<Number> > Lu1 = this->solution->clone();
	AutoPtr<NumericVector<Number> > tmp = this->solution->clone();
	//NumericVector<Number>* y = this->request_vector("_temporary_solution");

	// Stage one, compute u1
	tmp->zero();
	u1->zero();
	K->vector_mult(*tmp, *un);
	M->vector_mult(*Lun, *tmp);
	Lun->scale(dt);
	u1->add(*un);
	u1->add(*Lun);
	
	tmp->zero();
	K->vector_mult(*tmp, *u1);
	M->vector_mult(*Lu1, *tmp);
	Lu1->scale(dt);
	
	this->solution->zero();
	this->solution->add(*un);
	this->solution->add(0.5, *Lun);
	this->solution->add(0.5, *Lu1);

	this->update();

	_count++;

}

void LevelSetSystem :: assemble(){
	this->assemble_stiffness();
	this->assemble_mass();
}
			
void LevelSetSystem :: assemble_stiffness(){
	
	const MeshBase& mesh = this->get_mesh();
	const unsigned int dim = mesh.mesh_dimension();

	const DofMap& dof_map = this->get_dof_map();
	
	
	FEType fe_type = this->variable_type(0);
	
	AutoPtr<FEBase> fe(FEBase::build(dim, fe_type));
	AutoPtr<FEBase> fe_elem_face(FEBase::build(dim, fe_type));
	AutoPtr<FEBase> fe_neighbor_face(FEBase::build(dim, fe_type));
	
	QGauss qrule(dim, fe_type.default_quadrature_order());
	QGauss qface(dim-1, fe_type.default_quadrature_order());
	
	fe->attach_quadrature_rule(&qrule);
	fe_elem_face->attach_quadrature_rule(&qface);
	fe_neighbor_face->attach_quadrature_rule(&qface);
	
	const std::vector<Real>& JxW = fe->get_JxW();
    const std::vector<std::vector<Real> >& N = fe->get_phi();	
    const std::vector<std::vector<RealGradient> >& B = fe->get_dphi();
	
	//const std::vector<Real>& JxW_plus = fe_elem_face->get_JxW();
    const std::vector<std::vector<Real> >&  N_plus = fe_elem_face->get_phi();
    const std::vector<std::vector<RealGradient> >& B_plus = fe_elem_face->get_dphi();
    
	//const std::vector<Real>& JxW_minus= fe_neighbor_face->get_JxW();
	const std::vector<std::vector<Real> >&  N_minus = fe_neighbor_face->get_phi();
    const std::vector<std::vector<RealGradient> >& B_minus = fe_neighbor_face->get_dphi();
	
    const std::vector<Point>& qface_normals = fe_elem_face->get_normals();
    const std::vector<Point>& qface_points = fe_elem_face->get_xyz();
     
    DenseMatrix<Number> Ke;  
    DenseMatrix<Number> Kne;
    DenseMatrix<Number> Ken;
    DenseMatrix<Number> Kee;
    DenseMatrix<Number> Knn;
	
	// Extract the velocity vector for the current point
	VectorValue<Number> v(dim);
	
	std::vector<unsigned int> dof_indices;
	std::vector<unsigned int> neighbor_dof_indices;

  
  	// Get the system time
	Real t = this->time;
  
  
    MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end(); 
   
    for ( ; el != end_el; ++el){
		
		const Elem* elem = *el;
		
		dof_map.dof_indices(elem, dof_indices);
		
		const unsigned int n_dofs = dof_indices.size();

		fe->reinit(elem);
		Ke.resize(n_dofs, n_dofs);

				// Extract a vector of quadrature x,y,z coordinates
		const vector<Point> qp_vec = fe->get_xyz(); 
		
		//
		for (unsigned int qp = 0; qp < qrule.n_points(); qp++){
			//v(0) = velocity.point_value(0, qp_vec[qp], t);
			//v(1) = velocity.point_value(1, qp_vec[qp], t);
			v = velocity->system().point_value(qp_vec[qp]);		
		
			for (unsigned int i = 0; i < n_dofs; i++){
				for (unsigned int j = 0; j < n_dofs; j++){
					Ke(i,j) += JxW[qp]*(N[i][qp] * (v*B[j][qp]));
				}
			}
		}
		
		// Add to the global stiffness matrix
		this->matrix->add_matrix(Ke, dof_indices);

		// Loop through the sides
		for (unsigned int k = 0; k < elem->n_sides(); k++){
	
			if (elem->neighbor(k) != NULL){
				
				AutoPtr<Elem> side(elem->build_side(k));
				const Elem* neighbor = elem->neighbor(k);
				
				std::vector<Point> qface_neighbor_point;
  
				std::vector<Point > qface_point;

				fe_elem_face->reinit(elem, k);

				qface_point = fe_elem_face->get_xyz();

				// Refinement dependant
				FEInterface::inverse_map (elem->dim(), fe->get_fe_type(), neighbor, qface_point, qface_neighbor_point);
	
	            fe_neighbor_face->reinit(neighbor, &qface_neighbor_point);

			
				dof_map.dof_indices(elem, dof_indices);	
				dof_map.dof_indices (neighbor, neighbor_dof_indices);
				
				const unsigned int n_dofs = dof_indices.size();
				const unsigned int n_neighbor_dofs = neighbor_dof_indices.size();
  
				Kne.resize (n_neighbor_dofs, n_dofs);
				Ken.resize (n_dofs, n_neighbor_dofs);
				Kee.resize (n_dofs, n_dofs);
				Knn.resize (n_neighbor_dofs, n_neighbor_dofs);
		
				// Perform quadrature
				for (unsigned int qp = 0; qp < qface.n_points(); qp++){
				
				//	v(0) = velocity.point_value(0, qp_vec[qp], t);
				//	v(1) = velocity.point_value(1, qp_vec[qp], t);
				//	v = velocity(qface_point[qp], t);
					v = velocity->system().point_value(qp_vec[qp]);	
					
					Number v_dot_n = v*qface_normals[qp];
					if (v_dot_n >= 0){
					
						for (unsigned int i = 0; i < n_dofs; i++){
							for (unsigned int j = 0; j < n_neighbor_dofs; j++){
								Ken(i,j) += JxW[qp]*(N_plus[i][qp] * N_minus[j][qp]) * v_dot_n;
							}
						}
					
						for (unsigned int i = 0; i < n_neighbor_dofs; i++){
							for (unsigned int j = 0; j < n_neighbor_dofs; j++){
								Knn(i,j) += JxW[qp]*(N_minus[i][qp] * N_minus[j][qp]) * v_dot_n;
							}
						}

					} else {
					
						for (unsigned int i = 0; i < n_dofs; i++){
							for (unsigned int j = 0; j < n_dofs; j++){
								Kee(i,j) += JxW[qp]*(N_plus[i][qp] * N_plus[j][qp]) * v_dot_n;
							}
						}
					
						for (unsigned int i = 0; i < n_neighbor_dofs; i++){
							for (unsigned int j = 0; j < n_dofs; j++){
								Kne(i,j) += JxW[qp]*(N_minus[i][qp] * N_plus[j][qp]) * v_dot_n;
							}
						}
						
					} // v \dot n if
				} // quadrature loop
				
				// Add the contributions of this side to the stiffness matrix 
				this->matrix->add_matrix(Kne, neighbor_dof_indices, dof_indices);
				this->matrix->add_matrix(Ken, dof_indices, neighbor_dof_indices);
				this->matrix->add_matrix(Kee, dof_indices);
				this->matrix->add_matrix(Knn, neighbor_dof_indices);
							
			} // neighbor if
			
		} // side loop
	
	} // element loop
	
	this->matrix->close();
	
} // assemble_stiffness()
			
void LevelSetSystem :: assemble_mass(){

	// Get a constant reference to the mesh object and get the dimension
	const MeshBase& mesh = get_mesh();
	const unsigned int dim = mesh.mesh_dimension();

	// Get a constant reference to the Finite Element type
	// for the first (and only) variable in the system.
	FEType fe_type = variable_type(0);

	// Build a Finite Element object of the specified type.  
	AutoPtr<FEBase> fe      (FEBase::build(dim, fe_type));
	AutoPtr<FEBase> fe_elem_face (FEBase::build(dim, fe_type));
	AutoPtr<FEBase> fe_neighbor_face (FEBase::build(dim, fe_type));

	// A Gauss quadrature rule for numerical integration.
	// Let the \p FEType object decide what order rule is appropriate.
	//QSimpson qrule (dim);
	//QSimpson qface (dim-1);
	QGauss qrule (dim,   fe_type.default_quadrature_order());
	QGauss qface (dim-1, fe_type.default_quadrature_order());

	// Tell the finite element object to use our quadrature rule.
	fe->attach_quadrature_rule      (&qrule);
	fe_elem_face->attach_quadrature_rule (&qface);
	fe_neighbor_face->attach_quadrature_rule (&qface);

	// Here we define some references to cell-specific data that
	// will be used to assemble the linear system.  We will start
	// with the element Jacobian * quadrature weight at each integration point.   
	const std::vector<Real>& JxW      = fe->get_JxW();
	const std::vector<Real>& JxW_face = fe_elem_face->get_JxW();

	// The element shape functions evaluated at the quadrature points.
	const std::vector<std::vector<Real> >& N = fe->get_phi();
	const std::vector<std::vector<Real> >& N_face = fe_elem_face->get_phi();

	// The element shape function gradients evaluated at the quadrature
	// points.
	const std::vector<std::vector<RealGradient> >& B = fe->get_dphi();

	// The XY locations of the quadrature points used for face integration
	//const std::vector<Point>& qface_points = fe_elem_face->get_xyz();

	// A reference to the \p DofMap object for this system.
	const DofMap& dof_map = get_dof_map();

	// Define data structures to contain the element mass matrix
	// and right-hand-side vector contribution.  
	MyDenseMatrix<Number> Me;			// element mass matrix

	// Define a vector for storing the velocity
	//VectorValue<Number> v(dim);
	
	// Storage for the degree of freedom indices
	std::vector<unsigned int> dof_indices;
	
	// Get the system time
	Real t = this->time;

	// Loop over all the elements in the mesh that are on local processor
	MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
	const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end(); 
		
	for ( ; el != end_el; ++el){    
	
		// Pointer to the element current element
		const Elem* elem = *el;
		
		// Update element length
		if (h > elem->hmin()){
			h = elem->hmin();
		}

		// Get the degree of freedom indices for the current element
		dof_map.dof_indices(elem, dof_indices);

		// Compute the element-specific data for the current element
		fe->reinit (elem);

		// Extract a vector of quadrature x,y,z coordinates
		const vector<Point> qp_vec = fe->get_xyz(); 

		// Zero the element matrices and vectors
		Me.resize (dof_indices.size(), dof_indices.size());

		// Compute the RHS and mass and stiffness matrix for this element (Me)
		for (unsigned int qp = 0; qp < qrule.n_points(); qp++){
			
			// Extract the velocity vector for the current point
			//v(0) = velocity.point_value(0, qp_vec[qp], t);
			//v(1) = velocity.point_value(1, qp_vec[qp], t);
			//v = velocity(qp_vec[qp], t);
			//v = velocity->system().point_value(qp_vec[qp]);	
			
			for (unsigned int i = 0; i < N.size(); i++){
				for (unsigned int j = 0; j < N.size(); j++){				
					Me(i,j) += JxW[qp]*(N[i][qp] * N[j][qp]);  // mass matrix 
				}
			}
		}
	
		// Invert the matrix
		Me.inverse();

		// Apply the local components to the global mass matrix
		this->request_matrix("_mass_matrix_inverse")->add_matrix(Me, dof_indices);
		//this->mass_matrix_inverse->add_matrix(Me, dof_indices);

	} // (end) for ( ; el != end_el; ++el)
	
	// Indicate that the mass matrix is complete
	//this->mass_matrix_inverse->close();
	this->request_matrix("_mass_matrix_inverse")->close();
	
} // (end) assemble_stiffness()

