function test_thermo

Te = 183;
Tm = 327;
cl = 0.1547;
cs = 0.1779;
kp = 0.31;
hf = 30.162;
h0 = (cs-cl)*Te + hf;
m  = -232.63;

Ti = 185;
Ci = 0.192;

T_liq = @(C) Tm + m*C;
T_sol = @(C) max([Tm+m/kp*C, Te]);

f     = @(T,C) 1 - 1/(1-kp)*(T-T_liq(C))/(T-Tm);

h_liq = @(C) cl*T_liq(C) + h0;
h_sol = @(C) cs*T_sol(C);
h_e   = @(C) f(Te,C)*hf + cs*Te;

T     = @(h,f) (h-f*h0)/(f*cl+(1-f)*cs);

B = @(T) (cl-cs)*(T-Te) + hf;
h1 = @(f,T) f * B(T) + cs*T;
h2 = @(f,T) T*(f*(cl-cs) + cs) + f*h0;
h3 = @(C,T) (h_sol(C)*B(T) - hf*cs*T)/(B(T)-hf);

disp('Temperature levels:');
disp(['  T_liq = ', num2str(T_liq(Ci))]);
disp(['  T_sol = ', num2str(T_sol(Ci))]);
disp(['  T_e   = ', num2str(Te)]);

disp('Enthalpy levels:');
disp(['  h_liq = ', num2str(h_liq(Ci))]);
disp(['  h_sol = ', num2str(h_sol(Ci))]);
disp(['  h_e   = ', num2str(h_e(Ci))]);

disp(['Initial temperature: ', num2str(Ti)]);
if Ti > T_liq(Ci);
    disp('Case 1');
    h = h1(1,Ti);
   
elseif Te < Ti && Ti <= T_liq(Ci);
    disp('Case 2');
    h = h2(f(Ti,Ci),Ti);

elseif T_sol(Ci) < Ti && Ti <= Te;
    disp('Case 3');
    h = h3(Ci,Ti);

elseif Ti < T_sol(Ci);
   disp('Case 4');
    h = Ti*cs;
    
end
disp(['Enthalpy: ', num2str(h)]);


if h > h_liq(Ci);
    disp('Case 1');
    Ti = (h - h0)/cl;
    
elseif h_e(Ci) < h && h <= h_liq(Ci);
    disp('Case 2');
    Ti = iterative(f,T,h,Ti,Ci);

elseif h_sol(Ci) < h && h <= h_e(Ci);
    disp('Case 3');
    Ti = Te;

elseif h < h_sol(Ci);
    disp('Case 4');
    Ti = h/cs;
end
disp(['Computed temperature: ', num2str(Ti)]);

function T = iterative(lever_rule,temperature,h,Ti,Ci)

nmax = getpref('thermo','temp_max_iter');             
min_err  = getpref('thermo','temp_min_error');

cnt = 0;
err = 1;
T   = Ti;
while err > min_err;
    f    = lever_rule(T,Ci);
    Tnew = temperature(h,f);
    err  = abs(Tnew - T)/T;
    cnt = cnt + 1;
    if cnt > nmax; break; end
end
