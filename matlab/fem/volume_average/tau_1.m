%% Compute \(\tau_1^e\)
% Computes the convective stabilization term, \(\tau_1^e\), from Eq. 63 of
% Zabaras and Samanta (2004).

%% Function definition
function tau = tau_1(p, varargin)
% Computes the tau_1 term from Eq. 63.

% Computes the convective stabilization term, \(\tau_1^e\), from Eq. 63 of
% Zabaras and Samanta (2004).
%
% Syntax:
%   tau = tau_1(idx), where idx is the node index that ranges from 1 to 4
%   tau = tau_1(p) computes the value at the perscribed point, where p
%       gives the locatoin (e.g., p = [-1,-1]).
%   tau = tau_1(..., false) disables the printing of the results


%% Set the flag for displaying output
show = true;
if nargin == 2 && islogical(varargin{1});
    show = varargin{1};
elseif nargin > 1;
    error('Input Error');
end

%% Collect the parameters for testings
[C, V] = parameters;

%% Get the velocity vector at the desired location
if isscalar(p);
    v = V(p,:);     % value at node
else
    v = quad4(p,V); % value at supplied point
end

%% Compute the \(\tau_{SUPG}\) term

% Element length
h = elem_length;

% Eq. 67
n = norm(v);
Re_v = n*h / (2 *C.Pr);

% Eq. 70
if Re_v <= 3;
    zRe = Re_v/3;
else
    zRe = 1;
end

% Eq. 65
if n == 0;
    tau_supg = 0;
else
    tau_supg = C.epsilon*h / (2*n) * zRe;
end

%% Return the minimum of the two terms

% Case when tau_1b is infinite
if C.epsilon == 1;
    tau = tau_supg;
    if show;
        disp(['tau_1a = inf; tau_1b = ', num2str(tau)]);
    end
    
% Standard case    
else
    tau_epsilon = C.epsilon^2 / (1 - C.epsilon)^2 * C.Da / C.Pr;
    tau = min(tau_epsilon, tau_supg);
    if show;
        disp(['tau_1a = ', num2str(tau_epsilon),'; tau_1b = ', num2str(tau_supg)]);
    end
end

