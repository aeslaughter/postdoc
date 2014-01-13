 /*! \file thermo_system.h
  * \brief Header file for the ThermoSystem class.
  */
  
// Avoid multiple includes
#ifndef volume_average_thermo_system_h
#define volume_average_thermo_system_h  
   
// C++ standard library includes  
#include <iostream>
#include <string>
#include <vector>
using std::string;
using std::vector;

// BOOST includes
#include <boost/shared_ptr.hpp>

// libMesh includes
#include <libmesh.h>
#include <libmesh_common.h>
#include <point.h>
#include <elem.h>
#include <mesh.h>
#include <mesh_generation.h>
#include <system.h>
#include <equation_systems.h>
#include <explicit_system.h>
#include <nonlinear_implicit_system.h>
#include <transient_system.h>
#include <analytic_function.h>
#include <function_base.h>
#include <fe.h>
#include <fe_type.h>
#include <dof_map.h>
#include <quadrature_gauss.h>
#include <point_locator_tree.h>
#include <auto_ptr.h>
#include <vector_value.h>
#include <dense_vector.h>

// My includes
#include "fem/common/explicit_system_base.h"

// Add to the VolumeAverage namespace
namespace SlaughterFEM{

//! Boost shared pointer to SystemBase for linking equations
typedef boost::shared_ptr<SystemBase<TransientNonlinearImplicitSystem> > base_ptr;

/*! \class ThermoSystem thermo_system.h "fem/volume_average/thermo_system.h"
 * \brief A class for solving the volume averaged momentum equation with libMesh
 * 
 * This class contains the necessary components for computing the various
 * element and nodal parameters for the volume averaging finite element 
 * technique. All references to equations in the docmentation of this class,
 * unless noted otherwise, it taken from Samanta and Zabaras (2005),
 * "Modelling convection in solidification process using stabilized 
 * finite element techniques."
 * 
 * \ingroup FEMVolumeAverage
 */  	
class ThermoSystem : public ExplicitSystemBase{
	friend class EnergySystem;
	
	public:
		// EnergyEq is a friend b/c it relies on this class to initialize
		// enthalpy from temperature
	
	
		// Calls to initialized() are from EqCore class
		using SystemBase::initialized;
		
		//! Class constructor 
		/*!
		 * This is the general class constructor that should be called.
		 * There is also a private version used when cloning this class
		 * via the clone() function.
		 * 
		 * In order for this class to behave correctly, three pointers
		 * must be set: velocity_ptr, enthalpy_ptr, and concentration_ptr.
		 * 
		 * \param sys The over-riding EquationSystems object that contains
		 * the various Systems, including the data in this class.
		 * 
		 * \see clone() 
		 * \see MomentumEq
		 * \see EnergyEq
		 * \see ConcentrationEq
		 */
		ThermoSystem(EquationSystems& es, const string& name, const unsigned int number);
		
		//! Called when the variables are projected
		/*!
		 * The definition of the pure virtual function derived from the 
		 * FunctionBase class of libmesh.
		 * 
		 * \param index The variable to compute
		 * \param p The Point at which to compute the value
		 * \param t The time, which should not be used because libMesh
		 * 		uses t = 0 when projecting. If time is needed it should
		 * 		be extracted from the System.
		 * \return A libmesh:Number for the variable at the specified point
		 * \see operator()
		 */ 
		Number component(unsigned int index, const Point& p, Real t = 0);	
			
		//! Returns a vector of all the system variables
		/*!
		 * The definition of the pure virtual function derived from the 
		 * FunctionBase class of libmesh.
		 * 
		 * Computes the value of a variables associated with the System and
		 * returns them in the output variable.
		 * 
		 * \param p The Point at which to compute the value
		 * \param t The time, which should not be used because libMesh
		 * 		uses t = 0 when projecting. If time is needed it should
		 * 		be extracted from the System.
		 * \param output The vector in which the data is returned
		 * \see componenent
		 */ 	
		void operator() (const Point& p, const Real t, DenseVector<Number>& output); 			
			
		//! This function clones the class, required for libmesh
		/*!
		 * The definition of the pure virtual function derived from the 
		 * FunctionBase class of libmesh.
		 * 
		 * Creates a copy of this class, it automatically copies the
		 * velocity_ptr, concentration_ptr, and enthalpy_ptr pointers.
		 * 
		 * \return A libmesh::AutoPtr to a new instance of this class.
		 */ 
		AutoPtr<FunctionBase<Number> > clone () const;

		//! Computes the temperature at a point
		/*! 
		 * \f[
		 * T = \begin{cases} 
				\frac{h-h_0}{c_s} &\textrm{if } h > h_{liq} \\[1em] 
				\frac{h-h_0}{fc_l + (1-f)c_s} &\textrm{if } h_e < h \leq h_{liq} \\[1em] 
				T_e &\textrm{if } h_{sol} < h \leq h_e \\[1em] 
				\frac{h}{c_s} &\textrm{if } h \leq h_{sol} \\[1em] 
			\end{cases}
		 * \f]	
		 * \param p The current location (libMesh::Point)
		 * \return The scalar value of temperature
		 */ 
		Number temperature(const Point& p);
			
		//! Compute the volume average density
		/*! 
		 * \f[
		 * \rho = \begin{cases} 
				\rho_f &\textrm{if } h > h_{liq} \\[1em] 
				\frac{1}{\frac{f}{\rho_f} + \frac{1-f}{\rho_s}}\,(f \textrm{ via Eq. 22})
					&\textrm{if } h_e < h \leq h_{liq} \\[1em] 
				\frac{1}{\frac{f}{\rho_f} + \frac{1-f}{\rho_s}}\,(f \textrm{ via Eq. 12}) &\textrm{if } h_{sol} < h \leq h_e \\[1em] 
				\rho_s &\textrm{if } h \leq h_{sol} \\[1em] 
			\end{cases}
		 * \f]	
		 * \param p The current location (libMesh::Point)
		 * \return The scalar value of temperature
		 */ 
		Number density(const Point& p);
			
		//! A function for returning the volume fraction
		/*! 
		 * \f[
		 * \epsilon = \begin{cases} 
				1 &\textrm{if } h > h_{liq} \\[1em] 
				\frac{\rho f}{\rho_f}\,(f \textrm{ via Eq. 22}) &\textrm{if } h_e < h \leq h_{liq} \\[1em] 
				\frac{\rho - \rho_s}{\rho_s - \rho_f}&\textrm{if } h_{sol} < h \leq h_e \\[1em] 
				0 &\textrm{if } h \leq h_{sol} \\[1em] 
			\end{cases}
		 * \f]	
		 * \param p The current location (libMesh::Point)
		 * \return The scalar value of fluid volume fraction
		 */ 
		Number epsilon(const Point& p);		

		//! Returns the fluid concentration, \$f C_l\$f.
		/*!
		 * Refer to Section 4 of Zabaras and Samanta (2004) for details.
		 * 
		 * \f[
		 * C_l = \begin{cases}
		 * 		C & \textrm{if } h > h_{liq} \\[1em] 
		 *		\frac{T - T_m}{m} & \textrm{if } h_e < h \leq h_{Liquidus} \\[1em] 
		 *		C_e & \textrm{if } h_{Solidus} < h \leq h_e \\[1em]
		 * 		0 & \textrm{if } h_{Solidus} < h \leq h_e \\
		 * \end{cases}
		 * \f]
		 * \param p The current location (libMesh::Point)
		 * \return The scalar value of fluid concentration
		 */
		Number fluid_concentration(const Point& p);
		
		//! Liquid mass fraction, \$f f \$f.
		/*!
		 * Refer to Section 4 of Zabaras and Samanta (2004) for details.
		 * 
		 * \f[
		 * f = \begin{cases}
		 * 		1 &\textrm{if}\, h > h_{liq} \\[1em] 
		 *		1 - \frac{1}{1 - \kappa_p} \left(\frac{T - T_{liq}}{T - T_m}\right) 
					&\textrm{if}\, h_e < h \leq h_{Liquidus} \\[1em] 
		 *		\frac{h-h_{sol}}{h_f} &\textrm{if}\, h_{Solidus} < h \leq h_e \\[1em]
		 * 		0 &\textrm{if}\, h_{Solidus} < h \leq h_e \\
		 * \end{cases}
		 * \f]
		 * \param p The current location (libMesh::Point)
		 * \return The scalar value of liquid mass fraction
		 */
		Number liquid_mass_fraction(const Point& p);
		
		//! Volume average specific heat
		/*!
		 * \param p The Point at which to compute the value
		 * \param t The time, which should not be used because libMesh
		 * 		uses t = 0 when projecting. If time is needed it should
		 * 		be extracted from the System.
		 * \return Scalar of volume average specific heat
		 */ 
		Number specific_heat(const Point& p, const Real t);
		
		//! Volume average thermal conductivity
		/*!
		 * \param p The Point at which to compute the value
		 * \param t The time, which should not be used because libMesh
		 * 		uses t = 0 when projecting. If time is needed it should
		 * 		be extracted from the System.
		 * \return Scalar of volume average themal conductivity
		 */ 
		Number conductivity(const Point& p, const Real t);

		//! Computes the element length, Eq. 69 (Zabaras & Samanta, 2004)
		/*!
		 * Eq. 69 is misleading, the equation that should be examined
		 * is in Tezduyar (1992), Eq. 4.11. This is the relationship
		 * shown here.
		 * \f[
		 * h = 2 \left( \sum_{a=1}^{n_{en}} \left | \hat{v} \cdot \nabla N_{a} \right |   \right)^{-1}
		 * \f]
		 * \param elem A pointer to the libMesh::Elem for which the length is to be computed
		 * \return Scalar of the element length
		 */ 
		Number element_length(const Elem* elem);
			
		//! Returns the \f$\tau_1^e \f$ value for the advective stabilization term (Eq. 40)
		/*!
		 * \f[
		 * \tau_1^e = \min_{\vec{x}\in \Omega^e}\left [ \tau_{SUPG}, \frac{K(\epsilon)\rho_f}{\epsilon \mu} \right ]
		 * \f]
		 * 
		 * \param p The Point at which to compute the value
		 * \param h The element length for the element containing the 
		 * point specified.
		 * \return The scalar value of the advective stabilization term
		 * 
		 * \see tau_supg
		 * \see tau_epsilon
		 * \see element_length
		 */ 
		Number tau_1(const Point& p, const Number h);		
		
		//! A test function for this class
		/*!
		 * This function will print the results of the various equations
		 * in this class. These results should be compared to the MATLAB
		 * program: \ref matlab_volavg_thermo .
		 * \param elem A pointer for the element to test
		 */ 
		void test(Elem* elem);
		
		//! Pointer to the velocity variables
		base_ptr momentum;
		
		//! Pointer the enthalpy variable
		base_ptr energy;
		
		//! Pointer to the concentration variable
		base_ptr concentration;
		
	private:
	
		//! Initialization function
		void initialize();
				
		//! Set the default values for all of the required constants
		/*! 
		 * This functions sets the default values for the various constants
		 * required. The values listed in Tables I of Samanta and Zabaras (2005)
		 * are used by default.
		 */ 
		void set_default_constants();
						
		//! Computes the SUPG term for \f$ \tau_1^e \f$ (i.e., Eq. 42)
		/*!
		 * \f[
		 * \tau_{SUPG} = \frac{f h}{2 || \vec{v}^h ||} z(Re_{\vec{v}})
		 * \f]
		 * 
		 * \param p The Point at which to compute the value
		 * \param h The element length for the element containing the 
		 * point specified.
		 *
		 * \see tau_1
		 * \see element_length
		 */
		Number tau_supg(const Point& p, const Number h);
		
		//! Computes the \f$K\f$-based term for \f$ \tau_1^e \f$ in Eq. 40
		/*!
		 * \f[
		 * \tau_K = \frac{K(\epsilon)\rho_f}{\epsilon \mu}
		 * \f]
		 * 
		 * \param p The Point at which to compute the value
		 * \return The scalar value of \f$\tau_K\f$ for Eq. 40
		 * \see tau_1
		 * \see element_length
		 */
		Number tau_K(const Point& p);	
		
		//! Equation 6
		/*!
		 * \f[
		 * K(\epsilon) = \frac{K_0 \epsilon^3}{(1-\epsilon)^2}
		 * \f]
		 * \param p The Point at which to compute the value
		 * \param e The value of \f$\epsilon\f$ for the current point
		 * \return The scalar value of the Kozeny-Carman relationship
		 */
		Number kozeny_carman(const Point& p, const Number e);
		
		//! Equation 11
		/*!
		 * \f[
		 * h_0 = (c_s - c_f)T_e + h_f;
		 * \f]
		 * \return Returns the constant reference enthalpy
		 */ 
		Number reference_enthalpy();
		
		//! Equation 17
		/*!
		 * \f[
		 * T_{liq} = T_m + m_{liq}C
		 * \f]
		 * 
		 * \param p The current point on the element
		 * \return \f$T_{liq}\f$ [C] 
		 */
		Number T_liq(const Point& p);
		
		//! Equation 18
		/*!
		 * \f[
		 * T_{sol} = \textrm{max}\left\{T_m + \frac{m}{\kappa_p}C, T_e\right\}
		 * \f]
		 * \param p The current point on the element
		 * \return \f$T_{sol}\f$ [C]
		 */
		Number T_sol(const Point& p); 

		//! Equation 19
		/*!
		 * \f[
		 * h_{liq} = c_f T_{liq} + h_{l,0}
		 * \f]
		 * \param p The current point on the element
		 * \return \f$h_{liq}\f$ [kJ/kg]
		 */
		Number h_liq(const Point& p); 

		//! Eqution 20
		/*!
		 * \f[
		 * h_{sol} = c_s T_{sol}
		 * \f]
		 * \param p The current point on the element
		 * \return \f$h_{sol}\f$ [kJ/kg]
		 */
		Number h_sol(const Point& p);		
		
		//! Equation 21
		/*!
		 * \f[
		 * h_e = f_e h_f + c_s T_e
		 * \f]
		 * \param p The current point on the element
		 * \return \f$h_e\f$ [kJ/kg]
		 */
		Number h_e(const Point& p);
		
		//! Equation 22 
		/*!
		 * \f[
		 * f = 1 - \frac{1}{1 - \kappa_p} \left(\frac{T - T_{liq}}{T - T_m}\right)
		 * \f]
		 * \param p The current point on the element
		 * \param T The temperature at which to compute the value of f
		 * \return The liquid mass fraction, \f$f\f$ 
		 */
		Number lever_rule(const Point& p, const Number T);		
			
		//! Iterative solver for temperature when \f$ h_e < h \leq h_{liq} \f$
		/*!
		 * Solves f (Eq. 22) and T iterativily, as detailed on p. 1774
		 * of Samanta and Zabaras (2005).
		 * \param p The current point on the element
		 * \return The scalar value of f of Eq. 22
		 */ 
		Number temperature_iterative(const Point& p);		

		//! A private constructor function.
		/*!
		 * Performs the general constructor operations that is used
		 * by both constructors available for this class.
		 */ 
		void constructor();
	
		//! A protected class constructor used when the system is cloned
		/*!
		 * \param sys The over-riding EquationSystems object that contains
		 * the various Systems, including the data in this class.
		 * \param m Pointer to the MomentumEq containing velocity as the link variable
		 * \param h Pointer to the EnergyEq containing enthalpy as the link variable
		 * \param c Pointer to the ConcentrationEq containing concentration as the link variable
		 */ 
		ThermoSystem(EquationSystems es, base_ptr m, base_ptr h, base_ptr c);

		//! Validates that a pointer exists and is initialized
		/*!
		 * \param x A pointer to the EqBaseLink class being tested
		 * \param name the name of the class being tested (for error message)
		 * \param func the function from which the valid function is called
		 */
		void valid(boost::shared_ptr<SystemBase<TransientNonlinearImplicitSystem> > x, 
			string name, string func);  

}; // class
} // namespace

#endif
