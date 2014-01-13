%% Thermodynamic Relationships Test Function
% Section 4 of Zabaras and Samanta (2004) details the thermodynamic
% relationships that relate temperature (\(\theta\)), volume fraction
% (\(\epsilon\)), and liquid concentration (\(C_l\)). This function is
% designed to test that these relationships are working properly. All
% equation references refer to Samanta and Zabaras (2005).

%% Function definition
function varargout = thermodynamics(varargin)
% THERMODYNAMICS
% Test for thermodynamic relationships of Samanta and Zabaras 2005
%
% Syntax:
%   h = thermodynamics(T,C);

%% Gather the constants and element data
% C = liquid concentration (non-dimensional)
% T = nodal temperatuer (non-dimensional)
% [T, C, V] = parameters;
% p = quad4;

%% Return the value of a single function
    if nargin == 1 && nargout == 1;
        [~,C] = parameters;
        for i = 1:length(C);
            varargout{1}(i) = feval(varargin{1}, C(i));
        end
        return;
    end

%% Compute the enthalpy from temperature
    if nargin == 2 && nargout == 1;
        T = varargin{1};
        C = varargin{2};
        for i = 1:length(C);
            h(i) = temperature_to_enthalpy(T(i),C(i));
        end
        varargout{1} = h;
        return;
    end

%% Compute temperature(T), epsilon(E), density(D), fluid concentration(C), and mass fraction(F)
    if nargout == 5 && nargin == 3;

        % Initial values of enthalpy, concentration, and velocity
        T_int = varargin{1};
        C_int = varargin{2};
        V_int = varargin{3};

        % Loop through each of the nodes
        for i = 1:lenght(h);

            % enthalpy, concentration, and velocity for current node
            h = temperature_to_enthalpy(T_int(i),C_int(i));
            c = C_int(i);
            v = V_int(i);

            % Collect necessary constants
            h0 = reference_enthalpy;
            hf = getpref('thermo', 'latent_heat');
            cf = getpref('thermo', 'specific_heat_fluid');
            Tm = getpref('thermo', 'melting_temperature');
            Te = getpref('thermo', 'eutectic_temperature');
            m  = getpref('thermo','liquidus_slope');
            rho_f = getpref('thermo', 'density_fluid');
            rho_s = getpref('thermo', 'density_solid');

            if h > h_liq(c);
                T(i) = (h - h0) / cf;
                P(i) = getpref('thermo','density_fluid');
                E(i) = 1;
                C(i) = c;
                F(i) = 1;

            elseif h_e(c) < h && h <= h_liq(c);
                T(i) = temperature_iterative(h,T_int(i),c);
                F(i) = lever_rule(T(i),c);
                P(i) = 1 / (F(i)/rho_f + (1-F(i))/rho_s);
                E(i) = P(i)*F(i)/rho_f;
                C(i) = (T(i) - Tm)/m;

            elseif h_sol(c) < h && h <= h_e(c);
                T(i) = Te;
                C(i) = (Te - Tm)/m;
                F(i) = (h - h_sol(c))/hf;
                P(i) = 1 / (F(i)/rho_f + (1-F(i))/rho_s);
                E(i) = (P(i) - rho_s)/(rho_f - rho_s);

            elseif h <= h_sol(c);
                T(i) = h/cs;
                C(i) = 0;
                P(i) = rhos_s;
                F(i) = 0;
                E(i) = 0;
            end % if statement
        end % loop
        return;
    end
end % main function

%% Compute and dislpay the nodal values
% for i = 1:length(C);
%     
% %     % Display the current point
% %     str = sprintf('ELEMENT 0; Point %i (%g, %g)', i, p(i,1), p(i,2)); disp(str);
% %     out('reference_enthalpy', reference_enthalpy);
% %     out('T_liq', T_liq(C(i)));
% %     out('T_sol', T_sol(C(i)));
% %     out('h_liq', h_liq(C(i)));
% %     out('h_sol', h_sol(C(i)));
% %     out('h_e', h_e(C(i)));
%     
%     
%     
%     
%     
%     
% end 


%% Iterative temperature solver
function T = temperature_iterative(h,T,C)
    
    % Iteration control
    nmax = getpref('thermo','temp_max_iter');
    min_err  = getpref('thermo','temp_min_error');	

    % Gather constants
    h0 = reference_enthalpy;
    cf = getpref('thermo', 'specific_heat_fluid');
    cs = getpref('thermo', 'specific_heat_solid');

    % Iteration
    for i = 1:nmax;
        f = lever_rule(C,T);
        Tnew = (h - f*h0) / (f*cf + (1-f)*cs);
        err = abs(Tnew - T)/T;
        T = Tnew;
        if err <= min_err; break; end
    end % loop

end % function

%% Eq. 6

%% Eq. 11
function x = reference_enthalpy
    Te = getpref('thermo', 'eutectic_temperature');
    hf = getpref('thermo', 'latent_heat');
    cs = getpref('thermo', 'specific_heat_solid');
    cf  = getpref('thermo', 'specific_heat_fluid');
    x = (cs - cf)*Te + hf;
end

%% Eq. 17
function x = T_liq(C)
    Tm = getpref('thermo', 'melting_temperature');
    m  = getpref('thermo', 'liquidus_slope');
    x = Tm + m * C;
end

%% Eq. 18
function x = T_sol(C)
    Tm = getpref('thermo', 'melting_temperature');
    Te = getpref('thermo', 'eutectic_temperature');
    m  = getpref('thermo', 'liquidus_slope');
    kp = getpref('thermo', 'partition_coefficient');
    x = max(Tm + m/kp*C, Te);
end

%% Eq. 19
function x = h_liq(C)
    cf = getpref('thermo', 'specific_heat_fluid');
    h0 = reference_enthalpy;
    x = cf*T_liq(C) + h0;
end

%% Eq. 20
function x = h_sol(C)
    cs = getpref('thermo', 'specific_heat_solid');
    x = cs*T_sol(C);
end

%% Eq. 21
function x = h_e(C)
    hf = getpref('thermo', 'latent_heat');
    cs = getpref('thermo', 'specific_heat_solid');
    Te = getpref('thermo', 'eutectic_temperature');
    fe = lever_rule(C,Te);
    x = fe*hf + cs*Te;
end

%% Eq. 22
function f = lever_rule(C,T)
    kp = getpref('thermo', 'partition_coefficient');
    Tm = getpref('thermo', 'melting_temperature');
    f = 1 - 1/(1-kp)*(T-T_liq(C))/(T-Tm);
end

%% Convert initial temp. to initial enthalpy
function h = temperature_to_enthalpy(T,C)
    h0 = reference_enthalpy;
    hf = getpref('thermo', 'latent_heat');
    cf = getpref('thermo', 'specific_heat_fluid');
    cs = getpref('thermo', 'specific_heat_solid');
    Te = getpref('thermo', 'eutectic_temperature');
    f = lever_rule(C,T);   
    
    if T > T_liq(C);
       h = T*cf + h0;
    
    elseif T > Te && T <= T_liq(C);
        h = T*(f*cf + (1-f)*cs) + f*h0;
        
    elseif T > T_sol(C) && T <= Te;
        h = f*hf + h_sol(C);
        
    elseif T <= T_sol(C);
        h = T*cs;
    end  
end

%% Display function
function out(str,num)
    disp(['  ', str, ' = ', num2str(num)]);
end

%% References
% <html>
% <dd>Samanta, D. and Zabaras, N.</dd>
% <dd><i>Modelling convection in solidification process using stabilized finite element techniques.</i></dd>
% <dd>International Journal for Numerical Methods in Engineering, <b>2005</b>, Vol. 64, pp. 1769&ndash;1799.</dd>
% </html>