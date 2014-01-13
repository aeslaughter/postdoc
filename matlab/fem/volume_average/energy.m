%% Energy Equation Assembly

%%
% The energy equation
function energy

[~,qp] = quad4;
[Ti, Ci, Vi] = parameters;

Me = zeros(length(qp),length(qp));
Ne = zeros(length(qp),length(qp));

for i = 1:length(qp);
    [N,B,v] = quad4(qp(i,:), Vi);
%     tau = tau_1(qp(i,:),false);
% 
%     for j = 1:length(B);
%         P(j,1) = 1/C.epsilon * tau * v * B(:,j);
%     end
% 
% %     disp(['qp[',num2str(i),'] = ', num2str(qp(i,1)),', ',num2str(qp(i,2))]);
% %     disp(['  tau_1 = ', num2str(tau)]);
% %     disp(['  B = ', num2str(B(1,i)),', ', num2str(B(2,i))]);
% %     disp(['  v = ', num2str(v(1)),', ', num2str(v(2))]);
% %     disp(['  P = ', num2str(P)]);
% %     disp(' ');
%     
%     Me = Me + (N + P) * N';
%     
%     Ne = Ne + (N + P) * (v * B);
    
    
end
% 
% Me
% Ne

   