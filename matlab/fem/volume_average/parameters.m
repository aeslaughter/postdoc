%% Volume Average constants (test values)
% A function for that return an arbitrary value for testing the various
% calculations performed for the volume averaging finite element methods.

%% Function definition
function [T, C, V] = parameters
% PARAMETERS
% Returns various paramters for testing volume average calculations
%
% Syntax:
%   [T, C, V] = parameters;
%
% Description:


%% Define the constants
	% Physical constants (Table I of Samanta and Zabaras, 2005)
	setpref('thermo','conductivity_solid', 3.97e-2); 		% conductivity of solid, k_s [kW/m/C]
	setpref('thermo','conductivity_fluid', 2.29e-2);		% conductivity of fluid, k_l [kW/m/C]	
	setpref('thermo','specific_heat_solid', 0.1779);		% specific heat of solid, c_s [kJ/kg/C]
	setpref('thermo','specific_heat_fluid', 0.1547);		% specific heat of fluid, c_l[kJ/kg/C]
	setpref('thermo','latent_heat', 30.162);				% latent heat, h_f [kJ/kg]
	setpref('thermo','partition_coefficient', 0.31);		% equilibrium partition ratio,\kappa_p []
	setpref('thermo','thermal_expansion', 1.09e-4);		% thermal expansion coef., \Beta_T [1/C]
	setpref('thermo','solute_expansion', 0.354);			% solute expansion coef.. \Beta_s []
	setpref('thermo','density_solid', 10800);				% density of solid, \rho_s [kg/m^3]
	setpref('thermo','density_fluid', 10000);				% desnity of fluid, \rho_l [kg/m^3]
	setpref('thermo','viscosity', 0.0023);				% viscosity of fluid, \mu [kg/m/s]
	setpref('thermo','eutectic_temperature', 183);		% eutectic temperature, T_{eut} [C]
	setpref('thermo','melting_temperature', 327);			% melting temperature, T_m [C]
%	setpref('thermo','initial_temperature', 287);			% initial temperature, T_i [C]
	setpref('thermo','ambient_temperature', 20);			% ambient temperature, T_{amb} [C]
%	setpref('thermo','initial_concentration', 19.2);		% initial concentration, C_{l,0} [wt%]
	setpref('thermo','gravity', 9.81);					% gravity constant, g [m/s^2]
	setpref('thermo','liquidus_slope', -232.63); 			% dimensionless slope of liquidus, m_liq [C]
	setpref('thermo','diffusion', 1.05e-9);				% diffusion coefficient, D_l [m/s]
	setpref('thermo','dentrite_arm_spacing', 0.001);		% dentritic arm spacing, d [m]
	
	% Iteration parameters for solving Temperature
	setpref('thermo','temp_max_iter', 100);               % maximum number of iterations
	setpref('thermo','temp_min_error', 0.001);			% minium acceptable error

%% Define arbitrary nodal velocities
V(1,:) = [1,2];
V(2,:) = [1,0]; 
V(3,:) = [1,1];
V(4,:) = [0,0]; 

%% Define arbitrary intial concen tration
C = [0.192, 0.192, 0.192, 0.192];

%% Define arbitrary initial temperature
T = [287, 287, 287, 287];

