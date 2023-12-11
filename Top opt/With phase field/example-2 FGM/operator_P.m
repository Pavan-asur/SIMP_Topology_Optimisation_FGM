function [P_plus,P_minus]=operator_P(eps)

[V,D] = eig([eps(1),eps(3)/2; eps(3)/2,eps(2)]);
M1 = V(:,1)*V(:,1)'; M2 = V(:,2)*V(:,2)';
M12 = V(:,1)*V(:,2)'+ V(:,2)*V(:,1)';
lambda = diag(D);
g = (lambda+abs(lambda))/2;
beta = 0.5*(g(1)-g(2))/(lambda(1)-lambda(2));

d1 = 0; d2 = 0; 
if lambda(1)>0; d1=1; end
if lambda(2)>0; d2=1; end
if 0.98<lambda(1)/lambda(2)<1.02; beta = d1/2; end

T = d1*kron2(M1,M1)+d2*kron2(M2,M2)+beta*kron2(M12,M12);
P_plus = [T(1,1);T(3,3);2*T(1,3);T(2,2);T(4,4);2*T(2,4);T(1,2);T(3,4);2*T(1,4)];
P_minus = [1;0;0;0;1;0;0;0;1] - P_plus;

% 2nd order tensor Kron multiplication
function C = kron2(A,B)
C = [A(1,1)*B A(1,2)*B;
  A(2,1)*B A(2,2)*B];
  