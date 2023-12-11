function [d,dfixed,dfree,u,uload,uload1,uload2,ufree,f] = boundary_conditions(nelx,nely)

num_nodes = (nely+1)*(nelx+1);    % number of nodes

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
dfixed = (nely/2+1):(nely+1):((nely+1)*(nelx/4)+nely/2+1);
dfree = setdiff(dalldofs,dfixed);
d = zeros(num_nodes,1);
d(dfixed) = 1;