function example2
rho = 1600;
cp = 800;
k = 3;
h = 200;
Tinf = 50;
q = 0;
Tbc = 300;

% Bhatti numbering
X(:,1) = [0, 2, 2] / 100;
Y(:,1) = [0, 0, 4] / 100;
%X(:,2) = [0, 2, 0] / 100; % Bhatti
%Y(:,2) = [0, 4, 2] / 100; % Bhatti;
X(:,2) = [0, 2, 0] / 100;   % Box
Y(:,2) = [0, 4, 4] / 100;

% Convection and flux BC's
map(:,1) = [1, 3, 4]; bc{1} = [1, 2; 2, 3]; alpha{1} = [-h, 0]; beta{1} = [Tinf*h, -q];
map(:,2) = [1, 4, 2]; bc{2} = [1, 2; 2, 3]; alpha{2} = [0, 0]; beta{2} = [-q, 0];

n = 4;
theta = 0.5;
dt = 1;
tmax = 300;

    K = sparse(n,n);
    M = sparse(n,n);
    F = zeros(n,1);

    for i = 1:2;
        [Me, Kke, fe] =  compute_matrices(X(:,i), Y(:,i),k, rho, cp);
    fe   
        [Kae, fbe] = compute_boundary_components(X(:,i), Y(:,i), alpha{i}, beta{i}, bc{i});
    fbe 
        [M, K, F] = append_global(M, K, F, Me, Kke + Kae, fe + fbe, map(:,i));
    end
 
    Tall = [50;50;50;50];
    Mglobal = [M + theta*dt*K];
    Fglobal = dt*((1-theta)*F + theta*F) + (M - dt*(1-theta)*K)*Tall;
    
    Mglobal\Fglobal
    
    full(Mglobal)
    full(Fglobal)
    
    
    
    m = M([1,3], [1,3]);
    k = K([1,3], [1,3]); 
    f = F([1,3]) - K([1,3],[2,4])*[Tbc;Tbc];
    
    K_hat = m + theta*dt*k;
    F_hat_a = dt * ( (1-theta)*f + theta*f );
    F_hat_b = (m - dt*(1-theta)*k);
    
    
    
    
T = [50;50];
%Tout = zeros(tmax/dt + 1, 3);
Tout(1,2:3) = T;
j = 2;
for t = 1:dt:tmax
    F_hat = F_hat_a + F_hat_b * T;
    T = K_hat\F_hat;
    Tout(j,1) = t;
    Tout(j,2:3) = T;
    j = j + 1;
end

Tout
%Tout(1:30:end,:,:)
%plot(Tout(:,1), Tout(:,[2,3]));



function [Ka_out, fb_out] = compute_boundary_components(x, y, a, b, BC)

Ka_out = zeros(3,3);
fb_out = zeros(3,1);

for i = 1:length(BC);
    fb = zeros(3,1);
    Ka = zeros(3,3);
    
    
    bc = BC(i,:);
    L = sqrt( (x(bc(2)) - x(bc(1)))^2 + (y(bc(2)) - y(bc(1)))^2 );
    Ka(bc(1),bc(1)) = 2;
    Ka(bc(1),bc(2)) = 1;
    Ka(bc(2),bc(1)) = 1;
    Ka(bc(2),bc(2)) = 2;

    Ka_out = -a(i) * L / 6 * Ka + Ka_out;
    
    fb(bc(1)) = 1;
    fb(bc(2)) = 1;
    
    fb_out = b(i) * L / 2 * fb + fb_out;
    
end



function [M, K, f] = append_global(M, K, f, Me, Ke, fe, map)

for i = 1:3;
   f(map(i)) = f(map(i)) + fe(i);
    
    for j = 1:3;
        K(map(i),map(j)) = K(map(i),map(j)) +  Ke(i,j);
        M(map(i),map(j)) = M(map(i),map(j)) +  Me(i,j);
    end
end

function [M, K, f] = compute_matrices(x, y ,k, rho, cp)

b = [y(2) - y(3);
     y(3) - y(1);
     y(1) - y(2)];
     
c = [x(3) - x(2);
     x(1) - x(3);
     x(2) - x(1)];
 
f = [x(2)*y(3) - x(3)*y(2);
     x(3)*y(1) - x(1)*y(3);
     x(1)*y(2) - x(2)*y(1)];
 
A = 1/2 * sum(f);

for i = 1:3;
    N(i) = 1 / (2*A) * (x(i)*b(i) + y(i)*c(i) + f(i));
end

B = 1/(2*A)*[b, c];

M = rho * cp * A/12 * [2, 1, 1; 1, 2, 1; 1, 1, 2];

C = [k, 0; 0, k];

K = A*B*C*B';

f = [0,0,0]';