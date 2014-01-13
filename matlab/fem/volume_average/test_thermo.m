%% Test volume average thermodynamic relationships
% The C++ class ThermoEq may be tested using the test_thermo C++ program,
% the results should be comparable to the results of this program.

%% Function definintion
function test_thermo

%% Initial Temp. to Enthalpy Converstoin
    [T,C] = parameters;
    p = quad4;
    h = thermodynamics(T,C);

    fprintf('\n\nTEMP. TO ENTHALPY CONVERSION:\n');
    for i = 1:length(h);
        fprintf('\th = %g (%g,%g)\n', h(i), p(i,1), p(i,2));
    end

%% Test Eqs. (17) to (21)
    T_liq = thermodynamics('T_liq');
    T_sol = thermodynamics('T_sol');
    h_liq = thermodynamics('h_liq');
    h_sol = thermodynamics('h_sol');    
    h_e   = thermodynamics('h_e'); 
    
    fprintf('\n\nEQS. (17) to (21):\n');
    for i = 1:length(p);
        fprintf('\tPoint %i (%g, %g)\n', i, p(i,1), p(i,2));
        fprintf('\t\tT_liq = %g\n', T_liq(i));
        fprintf('\t\tT_sol = %g\n', T_sol(i));
        fprintf('\t\th_liq = %g\n', h_liq(i));
        fprintf('\t\th_sol = %g\n', h_sol(i));
        fprintf('\t\th_e = %g\n', h_e(i));
    end
