%% Element Length Test Function
%%
%% Element length
% Solution to Eq. 4.11 of Tezduyar (1992, p. 16) and Eq. 69 of 
% Zabaras and Samanta (2004):
%
% \[ h = 2 \left( \sum_{a=1}^{n_{en}} | \vec{s} \cdot \nabla N_a | \right)^{-1} \].

%% Function definition
function h = elem_length
% ELEM_LENGTH
% Compute the element length ot a test element
%
% Syntax:
%   h = elem_length;
%
% Description:
%   h = elem_length return the value of the element length, h, for the test
%       4-node quadraterial element

%% Gather constant and test parameters
[P,~] = quad4;      % gets the node positions
[~,V] = parameters; % gets the nodal velocities

%% Compute summation components of the element length
% Creates a vector of length \(n_{en}\) that
% contains the value to be summed: \(\left|\vec{s} \cdot \nabla N_a\right|\). If the
% velocity is zero then so is the associated value.

% Initialize h
h = 0;

for i = 1:length(P);

    % Compute the norm of the velocity vector, if it is zero set the value
    % to zero and progress to the next node
    n = norm(V(i,:));
    if n == 0;
       continue;
    end
    
    % Compute the shape function gradient
    [~, B] = quad4(P(i,:));
        
    % Compute the normal vector from the velocity
    s = V(i,:)/ n;

    % Compute the summation component
    h = h + abs(dot(s,B(:,i)));
end

%% Compute the scalar element length, this is the value returned
h = 2/h;

%% References
% <html>
% <dd>Tezduyar, T.</dd>
% <dd><i>Stabilized finite element formulations for incompressible flow computations</i></dd>
% <dd>Advances in applied mechanics, <b>1992</b>, Vol. 28(1), pp. 1-44.</dd>
% </html>
