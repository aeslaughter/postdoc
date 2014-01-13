/*! \file thermo_system.cpp
 * \brief Source code for a class that stores volume average thermodynamic nodal data
 */

// Standard C++ includes
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>

// BOOST includes
#include <boost/current_function.hpp>

// Include the header for this source  
#include "fem/volume_average/thermo_system.h"
using namespace SlaughterFEM; 
  
ThermoSystem :: ThermoSystem(EquationSystems& es, const string& name, const unsigned int number):
	ExplicitSystemBase(es, name, number){};
	
/*
ThermoSystem :: ThermoSystem(EquationSystems es, base_ptr m, base_ptr h, base_ptr c) : 
	ExplicitSystemBase(es, "thermo_system"){) : 
	ExplicitSystemBase(es, "thermo_system"){

// A private class constructor used when the system is cloned
ThermoSystem :: ThermoSystem(EquationSystems es, base_ptr m, base_ptr h, base_ptr c) : 
	ExplicitSystemBase(es, "thermo_system"){
		
	// Call the general constructor function		
	constructor();	
	
	// Set the EqVariableLinker pointers
	this->momentum = m;
	this->energy = h;
	this->concentration = c;
}
*/
//A priviate class for constructing
void ThermoSystem :: constructor(){
				
	// Setup constants
	set_default_constants();	// sets the default values of constants
	
	// Add Variables
	// - All must have the same order (e.g., FIRST) 
	add_variable("temperature", FIRST);			// temperature
	add_variable("density", FIRST);				// vol. avg. density	
	add_variable("epsilon", FIRST);				// fluid volume fraction
	add_variable("fluid_concentration", FIRST);	// fluid concentration
	add_variable("liquid_mass_fraction", FIRST);// liquid mass fraction, Eq. 22
	
	// Attach initilization
	this->attach_init_object(*this);
	
	// Set initilization flag to false (EqCore)
	SystemBase::_initialized = false;
}		
				
// This is called when the solution is projected
Number ThermoSystem :: component(unsigned int index, const Point& p, Real t){

	// Determine the variable name
	string var_name = this->variable_name(index);
	
	// Compute the temperature
	if (var_name.compare("temperature") == 0){
		return temperature(p);
		
	// Compute the fluid volume fraction
	} else if (var_name.compare("epsilon") == 0){
		return epsilon(p);
	
	// Compute the concentration of the fluid
	} else if (var_name.compare("fluid_concentration") == 0){
		return fluid_concentration(p);

	// The volume average density	
	} else if (var_name.compare("density") == 0) {
		return density(p);
		
	// The liquid mass fraction
	} else if (var_name.compare("liquid_mass_fraction") == 0){
		return liquid_mass_fraction(p);
	} 
}		
	
		
// Creates a copy of this class		
AutoPtr<FunctionBase<Number> > ThermoSystem :: clone() const{
	return AutoPtr<FunctionBase<Number> >(new ThermoSystem(this->get_equation_systems(), momentum, energy, concentration));
}
	
// Initialize function
void ThermoSystem :: initialize(){

	// Check EqVariableLinker pointers, exist
	if (!momentum){
		printf("ERROR: A valid pointer for the velocity must be set.\n");
		libmesh_error();
	} else if (!energy){
		printf("ERROR: A valid pointer for the enthalpy must be set.\n");
		libmesh_error();
	} else if (!concentration){
		printf("ERROR: A valid pointer for the concentration must be set.\n");
		libmesh_error();
	}
	
	// Check that the EqVariableLinker pointers have been initilized
	if (!momentum->initialized()){
		printf("ERROR: The velocity must be initilized.\n");
		libmesh_error();
	} else if (!energy->initialized()){
		printf("ERROR: The enthalpy must be initilized.\n");
		libmesh_error();
	} else if (!concentration->initialized()){
		printf("ERROR: The concentration must be initilized.\n");
		libmesh_error();
	}	

	// Project the solution (uses the component member)
	this->project_solution(this);

	// Initially the old and current solutions are the same, this is
	// needed to compute time derivatives 
	*this->old_local_solution = *this->current_local_solution;

	// This ThermoSystem system is now initialized
	SystemBase::_initialized = true;
}	
	
// Sets the defaults for the various constants 	
void ThermoSystem :: set_default_constants(){
	
	// Physical constants (Table I of Samanta and Zabaras, 2005)
	set_constant<Number>("conductivity_solid", 3.97e-2); 		// conductivity of solid, k_s [kW/m/C]
	set_constant<Number>("conductivity_fluid", 2.29e-2);		// conductivity of fluid, k_l [kW/m/C]	
	set_constant<Number>("specific_heat_solid", 0.1779);		// specific heat of solid, c_s [kJ/kg/C]
	set_constant<Number>("specific_heat_fluid", 0.1547);		// specific heat of fluid, c_l[kJ/kg/C]
	set_constant<Number>("latent_heat", 30.162);				// latent heat, h_f [kJ/kg]
	set_constant<Number>("partition_coefficient", 0.31);		// equilibrium partition ratio,\kappa_p []
	set_constant<Number>("thermal_expansion", 1.09e-4);			// thermal expansion coef., \Beta_T [1/C]
	set_constant<Number>("solute_expansion", 0.354);			// solute expansion coef.. \Beta_s []
	set_constant<Number>("density_solid", 10800);				// density of solid, \rho_s [kg/m^3]
	set_constant<Number>("density_fluid", 10000);				// desnity of fluid, \rho_l [kg/m^3]
	set_constant<Number>("viscosity", 0.0023);					// viscosity of fluid, \mu [kg/m/s]
	set_constant<Number>("eutectic_temperature", 183);			// eutectic temperature, T_{eut} [C]
	set_constant<Number>("melting_temperature", 327);			// melting temperature, T_m [C]
	//set_constant<Number>("initial_temperature", 287);			// initial temperature, T_i [C]
	set_constant<Number>("ambient_temperature", 20);			// ambient temperature, T_{amb} [C]
	//set_constant<Number>("initial_concentration", 19.2);		// initial concentration, C_{l,0} [wt%]
	set_constant<Number>("gravity", 9.81);						// gravity constant, g [m/s^2]
	set_constant<Number>("liquidus_slope", -232.63); 			// dimensionless slope of liquidus, m_liq [C]
	set_constant<Number>("diffusion", 1.05e-9);					// diffusion coefficient, D_l [m/s]
	set_constant<Number>("dentrite_arm_spacing", 0.001);		// dentritic arm spacing, d [m]
	
	// Numerical parameters
	set_constant<Number>("dt", 0.01);							// time step [s]
	
	// Iteration parameters for solving Temperature
	set_constant<unsigned int>("temp_max_iter", 100);			// maximum number of iterations
	set_constant<Number>("temp_min_error", 0.001);				// minium acceptable error	
}

// Prints equation results for the given element
void ThermoSystem :: test(Elem* elem){
	
	printf("\nEQS. (17) to (21):\n");
	for (unsigned int i = 0; i < elem->n_nodes(); i++){
		Point& p = elem->point(i);
		printf("\t Point %d: (%g, %g)\n", i, p(0), p(1));
		printf("\t\tT_liq = %g\n", T_liq(p));
		printf("\t\tT_sol = %g\n", T_sol(p));
		printf("\t\th_liq = %g\n", h_liq(p));
		printf("\t\th_sol = %g\n", h_sol(p));
		printf("\t\th_e   = %g\n", h_e(p));			
	}		
}						

// Validation function for pointers
void ThermoSystem :: valid(boost::shared_ptr<SystemBase<TransientNonlinearImplicitSystem> > x, string name, string func){
		
	// Check that system is defined
	if (!x){
		printf("ERROR (%s): The %s class must exist.\n", func.c_str(), name.c_str());
		libmesh_error();
		
	// Check that the associated system is initialized	
	} else if (!x->initialized()){
		printf("ERROR (%s): The %s class must be initialized.\n", func.c_str(), name.c_str());
		libmesh_error();
	}
	
}

// Volume average specific heat
Number ThermoSystem :: specific_heat(const Point& p, const Real t){
	
	const Number cs = get_constant<Number>("specific_heat_solid");
	const Number cf = get_constant<Number>("specific_heat_fluid");
	const Number eps = epsilon(p);
	
	return eps * cf + (1 - eps) * cs;
}
	
// Volume average thermal conductivity heat
Number ThermoSystem :: conductivity(const Point& p, const Real t){
	
	const Number ks = get_constant<Number>("conductivity_solid");
	const Number kf = get_constant<Number>("conductivity_fluid");
	const Number eps = epsilon(p);
	
	return eps * kf + (1 - eps) * ks;
}
	
// Returns the element length
Number ThermoSystem :: element_length(const Elem* elem){
	
	// Set the FEType to first order lagrange
	FEType fe_type(FIRST, LAGRANGE); 

	// Build a Finite Element object of the specified type (ndim() from EqCore)
	AutoPtr<FEBase> fe (FEBase::build(ndim(), fe_type));

	// The element shape functions evaluated at the quadrature points.
	const vector<vector<RealGradient> >& B = fe->get_dphi();
		
	// Get the node locations associated with this element
	vector<Point> pvec;
	for (unsigned int i = 0; i < elem->n_nodes(); i++){
		pvec.push_back(elem->point(i));
	}
		
	// Re-initialize the fe system for this element
	fe->reinit(elem, &pvec);
		
	// Intialize the element length
	Number h = 0;
	 
	// Loop over Gauss points (nodal summation of Eq. 69)
	for (unsigned int i = 0; i < pvec.size(); i++){	

		// Extract velocity for current point and compute the norm
		VectorValue<Number> s(momentum->point_value(pvec[i]));
		Number norm = s.size();
		
		// Continue only if the norm is non-zero
		if (norm != 0){
			s = s*(1/norm); 			// normalize the velocity vector
			h += std::abs(s * B[i][i]);	// s \cdot \nabla N_a	
		}
	}

	// Return the value of the element length
	return 2/h;	
}

// Returns the \tau_1 value for the advective stabilization term (Eq. 40)
Number ThermoSystem :: tau_1(const Point& p, const Number h){
	
	// Return the minimum of the two values computed
	return std::min(tau_supg(p, h), tau_K(p));
}

// Returns the SUPG portion of the \tau_1 equation
Number ThermoSystem :: tau_supg(const Point& p, const Number h){

	// Test the momentum
	valid(momentum, "MomentumEq", BOOST_CURRENT_FUNCTION);
	
	// Get the necessary constants
	const Number mu = get_constant<Number>("viscosity");
	
	// Compute the norm of the velocity
	const VectorValue<Number> v(momentum->point_value(p));
	const Number v_norm = v.size();				// ||v^h||
	
	// Get the current temperature
	const Number T = this->point_value("temperature", p);
	
	// Get the density
	const Number rho = this->point_value("density", p);
	
	// Get the value for f
	const Number f = this->point_value("liquid_mass_fraction", p);
	
	// Compute Reynolds number from velocity
	const Number Re = v_norm * h / (2 * mu/rho); // Re_v, Eq. 67

	// z(Re) of Eq. 70
	Number zRe;
	if (Re >= 0 && Re <= 3){
		zRe = Re/3;
	} else {
		zRe = 1;
	}			

	// The \tau_{SUPG} value
	return f * h / (2 * v_norm) * zRe ;
}

// Returns the K portion of the \tau_1 equation
Number ThermoSystem :: tau_K(const Point& p){
	
	// Gather relevant constants
	const Number mu = get_constant<Number>("viscosity");
	const Number pf = get_constant<Number>("density_fluid");
	
	// Get the value of epsilon
	unsigned int idx = this->variable_number("epsilon");
	const Number eps = System::point_value(idx, p);
	
	// Get the value for Kozeny-Carmen (Eq. 6);
	const Number K = kozeny_carman(p, eps);
	
	// Return desired value
	return (K * pf) / (eps * mu);
}

// Eq. 6
Number ThermoSystem :: kozeny_carman(const Point& p, const Number e){
	
	// Get arm spacing
	const Number d = get_constant<Number>("dentrite_arm_spacing");

	// Return the desired value
	return (std::pow(d,2)/180 * std::pow(e,3)) / std::pow((1-e),2);
}

// Eq. 11
Number ThermoSystem :: reference_enthalpy(){
	
	// Collect the constants needed for computation
	const Number cs = get_constant<Number>("specific_heat_solid");	
	const Number cf = get_constant<Number>("specific_heat_fluid");	
	const Number hf	= get_constant<Number>("latent_heat");			
	const Number Te = get_constant<Number>("eutectic_temperature");
	
	// Compute the reference enthalpy, Eq. 11.
	return (cs - cf)*Te + hf;
}

// Eq. 17
Number ThermoSystem :: T_liq(const Point& p){
	
	// Test the concentration
	valid(concentration, "ConcentrationEq", BOOST_CURRENT_FUNCTION);
	
	// Gather necessary terms
	const Number Tm = get_constant<Number>("melting_temperature");
	const Number m = get_constant<Number>("liquidus_slope");
	const Number C = concentration->point_value(0, p);
	
	// Compute T_liq
	return Tm + m * C;
}

// Eq. 18
Number ThermoSystem :: T_sol(const Point& p){
	
	// Test the concentration
	valid(concentration, "ConcentrationEq", BOOST_CURRENT_FUNCTION);
	
	// Gather necessary terms
	Number Tm = get_constant<Number>("melting_temperature");
	Number Te = get_constant<Number>("eutectic_temperature");
	Number m  = get_constant<Number>("liquidus_slope");
	Number kp = get_constant<Number>("partition_coefficient");
	Number C = concentration->point_value(0, p);
	
	// Return the maximum
	return std::max(Tm + m/kp*C, Te);
}

// Eq. 19
Number ThermoSystem :: h_liq(const Point& p){
	
	// Gather constants
	const Number cf = get_constant<Number>("specific_heat_fluid");
	const Number h0 = reference_enthalpy();

	// Return the desired value
	return cf * T_liq(p) + h0;
} 

// Eq. 20
Number ThermoSystem :: h_sol(const Point& p){
	
	// Gather constants
	const Number cs = get_constant<Number>("specific_heat_solid");

	// Return the desired value
	return cs* T_sol(p);
}		

// Eq. 21
Number ThermoSystem :: h_e(const Point& p){
	
	// Gather the required constants
	const Number Te = get_constant<Number>("eutectic_temperature");
	const Number hf = get_constant<Number>("latent_heat");
	const Number cs = get_constant<Number>("specific_heat_solid");

	// Gather liquid mass fraction at Te
	Number fe = lever_rule(p,Te);

	// Return the desired value
	return fe * hf + cs * Te; 
}

// Eq. 22
Number ThermoSystem :: lever_rule(const Point& p, const Number T){
	
	// Gather the required constants
	const Number kp = get_constant<Number>("partition_coefficient");
	const Number Tm = get_constant<Number>("melting_temperature");

	// Return the value
	return 1 - 1 / (1-kp) * (T - T_liq(p)) / (T - Tm);
}

// Temperature
Number ThermoSystem :: temperature(const Point& p){
	
	// Test the energy
	valid(energy, "EnergyEq", BOOST_CURRENT_FUNCTION);
	
	// Get the necessary constants
	const Number cf = get_constant<Number>("specific_heat_fluid");
	const Number cs = get_constant<Number>("specific_heat_solid");
	const Number h0 = reference_enthalpy();
	
	// Get the various enthalpy parameters
	const Number hliq = h_liq(p);
	const Number hsol = h_sol(p);
	const Number he = h_e(p);	
	
	// Get the current enthalpy
	const Number h = energy->point_value(0, p);

	// Compute the correct value for theta
	if (h > hliq){
		return (h - h0) / cf;
	
	} else if (he < h && h <= hliq) {
		return temperature_iterative(p);
	
	} else if (hsol < h && h <= he){
		return get_constant<Number>("eutectic_temperature");
		
	} else if (h <= hsol){
		return h / cs;
	}
	
}

// Iterative solver for temperature
Number ThermoSystem :: temperature_iterative(const Point& p){
	
	// Check that energy system is initialized
	if (!energy){
		printf("ERROR: The EnergyEq class must be initialized.\n");
		libmesh_error();
	}
		
	// Gather/set the iteration related constants
	unsigned int n_max = get_constant<unsigned int>("temp_max_iter");	// maximum number of iterations
	Number err_min = get_constant<Number>("temp_min_error");	// minium acceptable error	
	unsigned int cnt = 0; // iteration counter
	
	// Get the necessary constants
	const Number cf = get_constant<Number>("specific_heat_fluid");
	const Number cs = get_constant<Number>("specific_heat_solid");
	const Number h0 = reference_enthalpy();
	
	// Get the various enthalpy parameters
	const Number hliq = h_liq(p);
	const Number hsol = h_sol(p);
	const Number he = h_e(p);	
	
	// Get the current enthalpy and temperature
	const Number h = energy->point_value(0, p);
	Number T = point_value("temperature", p);

	// Perform iterative solution
	Number err = 1;
	while (err > err_min){
		
		// Compute the value of f (lever-rule)
		Number f = lever_rule(p, T);

		// Compute the new temperature
		Number T_new = (h - f*h0) / (f*cf + (1-f)*cs);
		
		// Compute the error estimate
		err = std::abs(T_new - T) / T;
		
		// Update the epsilon term
		T = T_new;
		
		// Increment count; and exit if the max count is reached
		cnt++;
		if (cnt > n_max){ 
			printf("WARNING: Maximum iterations of %i reached when computing temperature.\n", cnt);
			break;
		}			
	}
	
	// Return the temperature
	return T;
} 

// Density
Number ThermoSystem :: density(const Point& p){

	// Test the energy
	valid(energy, "EnergyEq", BOOST_CURRENT_FUNCTION);
		
	// Get the necessary constants
	const Number hf = get_constant<Number>("latent_heat");
	const Number pf = get_constant<Number>("density_fluid");
	const Number ps = get_constant<Number>("density_fluid");
	
	// Get the various enthalpy parameters
	const Number hliq = h_liq(p);
	const Number hsol = h_sol(p);
	const Number he = h_e(p);	
	
	// Get the current enthalpy
	const Number h = energy->point_value(0, p);
	
	// Get the current temperature
	const Number T = temperature(p);
//	const Number T = point_value("temperature", p);

	// Compute the value of density (see p. 1774 of Samanta and Zabaras, 2005)
	if (h > hliq){
		return pf;
	
	// Compute the density via Eq. 16
	} else if (he < h && h <= hliq) {
		Number f = lever_rule(p, T);
		return 	1 / (f/pf + (1-f)/ps);
	
	// Melt in control volume
	} else if (hsol < h && h <= he){
		Number f = (h - hsol) / hf;
		return 1 / (f/pf + (1-f)/ps);
		
	} else if (h <= hsol){
		return ps;
		
	} else {
		printf("ERROR: Enthalpy, %f, is in an undefined range.\n", h);
		libmesh_error();
	}	
}

// Volume fraction of fluid
Number ThermoSystem :: epsilon(const Point& p){

	// Test the energy
	valid(energy, "EnergyEq", BOOST_CURRENT_FUNCTION);
		
	// Get the necessary constants
	const Number ps = get_constant<Number>("density_solid");
	const Number pf = get_constant<Number>("density_fluid");
	const Number hf = get_constant<Number>("latent_heat");
	
	// Get the various enthalpy parameters
	const Number hliq = h_liq(p);
	const Number hsol = h_sol(p);
	const Number he = h_e(p);	
	
	// Get the current temperature, density, and enthalpy
	const Number T = temperature(p);
	const Number rho = density(p);
	const Number h = energy->point_value(0, p);
	
	// Compute the correct value for epsilon 
	if (h > hliq){
		return 1;
	
	} else if (he < h && h <= hliq) {
		Number f = lever_rule(p, T);
		return rho * f / pf;
	
	} else if (hsol < h && h <= he){
		Number f = (h - hsol) / hf;
		return (rho - ps) / (pf - ps);	
				
	} else if (h <= hsol){
		return 0;
	}
}	

// Concentration of the fluid
Number ThermoSystem :: fluid_concentration(const Point& p){
	
	// Test the energy
	valid(energy, "EnergyEq", BOOST_CURRENT_FUNCTION);
		
	// Get the necessary constants
	Number m = get_constant<Number>("liquidus_slope");
	Number Tm = get_constant<Number>("melting_temperature");	
	Number Te = get_constant<Number>("eutectic_temperature");

	// Get the various enthalpy parameters
	Number hliq = h_liq(p);
	Number hsol = h_sol(p);
	Number he = h_e(p);	
	
	// Get the current temp., density, and enthalpy
	//! \todo Using point_value is causing a segmentation fault, why?
	//const Number T = point_value("temperature", p);
	//const Number rho = point_value("density", p);
	const Number T = temperature(p);
	const Number rho = density(p);
	const Number h = energy->point_value(0, p);
	
	// Compute the correct value for epsilon 
	if (h > hliq){
		return concentration->point_value(0, p);
	
	} else if (he < h && h <= hliq) {
		return (T - Tm) / m;
	
	} else if (hsol < h && h <= he){
		return (Te - Tm) / m;
				
	} else if (h <= hsol){
		return 0;
	}
}

// Liquid mass fraction
Number ThermoSystem :: liquid_mass_fraction(const Point& p){

	// Gather necessary variables
	const Number hf = get_constant<Number>("latent_heat");
	const Number h = energy->point_value(0, p);
	const Number T = this->point_value("temperature", p);
	const Number hliq = h_liq(p);
	const Number hsol = h_sol(p);
	const Number he = h_e(p);
	
	// Compute the correct value for epsilon 
	if (h > hliq){
		return 1;
	
	} else if (he < h && h <= hliq) {
		return lever_rule(p, T);
	
	} else if (hsol < h && h <= he){
		return (h - h_sol(p))/hf;
				
	} else if (h <= hsol){
		return 0;
	} else {
		return 0;
	}
	
}
