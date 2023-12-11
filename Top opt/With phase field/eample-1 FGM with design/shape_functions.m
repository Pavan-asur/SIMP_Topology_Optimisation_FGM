function [Nd,Ndp,Bdp,Bu,J] = shape_functions(elemsize)

J = elemsize*elemsize;
% gp4 gp3 [0,1]
% gp1 gp2 [0,1]
gauss = [0.2113 0.2113; 0.7887 0.2113;...
  0.7887 0.7887; 0.2113 0.7887];

Nd = zeros(4,4); Bd = zeros(8,4);
Nu = zeros(8,8); Bu = zeros(12,8);

for i = 1:4
  s = gauss(i,1); t = gauss(i,2);
  
  N1 = (1-s)*(1-t); N2 = s*(1-t);
  N3 = s*t; N4 = (1-s)*t;
  Nd(i,:) = [N1 N2 N3 N4];
  
  Bd1 = [t-1; s-1]; Bd2 = [1-t; -s];
  Bd3 = [t; s]; Bd4 = [-t; 1-s];
  Bd(2*i-1:2*i,:) = [Bd1 Bd2 Bd3 Bd4]./elemsize;
  
  Nu1 = [N1 0; 0 N1]; Nu2 = [N2 0; 0 N2];
  Nu3 = [N3 0; 0 N3]; Nu4 = [N4 0; 0 N4];
  Nu(2*i-1:2*i,:) = [Nu1 Nu2 Nu3 Nu4];
  
  Bu1 = [t-1 0; 0 s-1; s-1 t-1]; Bu2 = [1-t 0; 0 -s; -s 1-t];
  Bu3 = [t 0; 0 s; s t]; Bu4 = [-t 0; 0 1-s; 1-s -t];
  Bu(3*i-2:3*i,:) = [Bu1 Bu2 Bu3 Bu4]./elemsize;
end

Ndp1 = Nd(1,:)'*Nd(1,:); Ndp2 = Nd(2,:)'*Nd(2,:);
Ndp3 = Nd(3,:)'*Nd(3,:); Ndp4 = Nd(4,:)'*Nd(4,:);
Ndp = [Ndp1(:),Ndp2(:),Ndp3(:),Ndp4(:)];
Bdp1 = Bd(1:2,:)'*Bd(1:2,:); Bdp2 = Bd(3:4,:)'*Bd(3:4,:);
Bdp3 = Bd(5:6,:)'*Bd(5:6,:); Bdp4 = Bd(7:8,:)'*Bd(7:8,:);
Bdp = [Bdp1(:),Bdp2(:),Bdp3(:),Bdp4(:)];