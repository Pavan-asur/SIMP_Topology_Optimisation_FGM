function [Hn,epsilon_gauss] = Hn_evaluation(nelx,nely,u,Hn,Bu,edofMat,lambda,mu,Psi_c)
%% history functional evaluation
epsilon_gauss = reshape(Bu*u(edofMat)',3,4*nely*nelx);
eps_tr_plus = operator_trace_plus(epsilon_gauss);
eps_plus = zeros(3,4*nely*nelx);
for i = 1:4*nely*nelx % parfor
  eps_plus(:,i) = operator_plus(epsilon_gauss(:,i));
end
Psi_plus = (lambda/2).*eps_tr_plus.^2 + mu.*(eps_plus(1,:).^2+eps_plus(2,:).^2+0.5*eps_plus(3,:).^2); % 
Psi = (abs(Psi_plus-Psi_c)+(Psi_plus-Psi_c))/2; % Psi_c
Hn = max(reshape(Psi,4,nelx*nely),Hn);

% evaluation of the trace of the positive epsilon
function eps_tr_plus = operator_trace_plus(eps)
eps_tr_plus = ((eps(1,:) + eps(2,:)) + abs(eps(1,:)+eps(2,:)))/2;

% evaluation of the positive epsilon
function eps_plus = operator_plus(eps)
[V,D] = eig([eps(1),eps(3)/2; eps(3)/2,eps(2)]);
eps_plus = (D(1,1)+abs(D(1,1)))/2*kron(V(:,1),V(:,1)')...
  + (D(2,2)+abs(D(2,2)))/2*kron(V(:,2),V(:,2)');
eps_plus = [eps_plus(1,1); eps_plus(2,2); 2*eps_plus(1,2)];