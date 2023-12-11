function [d,dfixed,dfree,u,uload,uload1,uload2,ufree,f] = boundary_conditions(nelx,nely)

num_nodes = (nely+1)*(nelx+1);    % number of nodes
dist=48;
%% boundary conditions of displacement field for tension test
ualldofs = (1:2*num_nodes);       % dofs of displacement field
uload1 = (2:2*(nely+1):2*num_nodes);
uload2 = (2*(nely+1):2*(nely+1):2*num_nodes);
uload=union(uload1,uload2);

ufixed1 = (1:2:2*(nely+1)-1);
ufixed2 = (2*(nely+1)*nelx+1:2:2*num_nodes-1);
ufixed=union(ufixed1,ufixed2);


ufree = setdiff(ualldofs,union(ufixed,uload));
u = zeros(2*num_nodes,1); 
f = zeros(2*num_nodes,1);

%% boundary conditions of crack phase field
dalldofs = (1:num_nodes);
dfixed =  50:81:698;
dfixed=[dfixed,4080:81:4890 ];
%dfree = dalldofs
setdiff(dalldofs,dfixed);
d = zeros(num_nodes,1);
d(dfixed) = 1;


passive(50, 1:5)=1;
passive(30, 55:60)=1;


% dalldofs = (1:num_nodes);
% dfixed1=(nely+1)*80+1: (nely+1)*80+9;
% dfixed2=(nely+1)*80+(nely-8):(nely+1)*80+(nely) ;
% dfixed=union(dfixed1,dfixed2);
%dfixed = (nely/2+1):(nely+1):((nely+1)*(nelx/4)+nely/2+1);
dfree = setdiff(dalldofs,dfixed);
d = zeros(num_nodes,1);
d(dfixed) = 1;