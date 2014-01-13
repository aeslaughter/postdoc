%% 4-Node Quadrilateral Element
% Returns the shape functions, shapve function derivatives, and velocity
% vector for a four node quadrilateral element.

%% Function definition
function [N, B, varargout] = quad4(p,varargin)
% QUAD4
% 4-Node Quadrilateral
%
% Returns the shape functions, shapve function derivatives, and velocity
% vector for a four node quadrilateral element.
%
% Syntax:
%  [p,qp] = quad4;
%  [N, B] = quad4(p);
%  [N, B, v] = quad4(p, V);
%
% Description:
%  [p,qp] = quad4 returns the node location and Gauss quadrature points
%       for the element.
%  [N, B] = quad4(p) returns the shape functions and gradients evaluated
%       at the specified point.
%  [N, B, v] = quad4(p, V) same as above but also returns the velocity 
%       vector projected at the specified point given the nodal velocities
%       in V (see parameters.m).

%% Return the nodes and quadrature points
if nargin == 0;
    
    % Define local node locations (\(\xi, \eta\))
    p(1,:) = [-1, -1];
    p(2,:) = [1, -1];
    p(3,:) = [1, 1];
    p(4,:) = [-1, 1];

    % Define the Guass quadrature points
    qp(1,:) = [-0.57735,-0.57735];
    qp(2,:) = [0.57735, -0.57735];
    qp(3,:) = [-0.57735, 0.57735];
    qp(4,:) = [0.57735,  0.57735];

    % Return these values via N and B; halt the program
    N = p; B = qp;
    return;
end

%% Get the shape function and derivate matrices
% Uses the two sub-functions to compute the values of the shape functions
% and thier gradient at the supplied point. The functions names mimic those
% used by libMesh (see help on libMesh::FEBase).

N = get_phi(p(1),p(2));
B = get_dphi(p(1),p(2));

%% Compute the velocity at the given location
% If the user supplies the nodal velocities, then return the velocity
% vector at the point as computed using the shape functions.

if nargin == 2 && ~isscalar(varargin{1});
    NN = [N', zeros(size(N')); zeros(size(N')), N'];
    d = reshape(varargin{1}, numel(varargin{1}),1);
    varargout{1} = (NN*d)';
end

end % quad4

%% Subfunction: Element shape functions
function N = get_phi(xi,eta)
% GET_PHI
% Retruns the element shape functions evalauted at eta and xi

% Fish & Belytschko (2007), Eq. 7.27, p. 165
xiI = [-1,1,1,-1];
etaI = [-1,-1,1,1];
N = 1/4 * (1 + xiI * xi)' .* (1 + etaI * eta)';

end % get_phi

%% Subfunction: Element shape function gradients
function B = get_dphi(xi,eta)
% GET_DPHI
% Returns the shape function gradients evaluated at eta and xi

% Fish & Belytschko (2007), Ex. 8.2, p. 198
B(:,1) = 1/4 * [eta - 1;    xi - 1];
B(:,2) = 1/4 * [1 - eta;    -xi - 1];
B(:,3) = 1/4 * [1 + eta;    1 + xi];
B(:,4) = 1/4 * [-eta - 1;   1 - xi];

end % get_dphi

%% References
% <html>
% <dd>Fish, J. & Belytschko, T. </dd>
% <dd><i>A first course in finite elements. </i></dd>
% <dd>John Wiley and Sons, Ltd., <b>2007</b></dd>
% </html>
