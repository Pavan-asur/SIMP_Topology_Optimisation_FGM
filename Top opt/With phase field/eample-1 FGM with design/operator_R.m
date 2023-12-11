function [R_plus,R_minus] = operator_R(eps)
R_plus = (sign(eps(1,:)+eps(2,:))+1)./2; % R+
R_minus = (sign(-eps(1,:)-eps(2,:))+1)./2; % R-