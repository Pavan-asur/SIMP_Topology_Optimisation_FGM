% function eso_crack(nelx,nely,X0,rmin)
% nelx=60; nely=120; volfrac=0.05; er=0.04; rmin=1.667;
close all;
clear all;
clc

nelx=60; nely=80; X0=0.0512; 
%% initialization for topology optimization
x = zeros(nely,nelx); % 1-inclusion, 0-matrix
x(1:nely,1:nelx) = X0;
% vol = mean(x(:));
maxiter=200;
elsize = 100./nelx; % millimeter (mm)
l = elsize; % mm
rmin=elsize;
loadN=150;

passive = zeros(nely,nelx);
passive(49:51, 1:9)=1;
passive(29:31, 51:60)=1;
% for i=1:nelx
%     for j=1:nely
%         if sqrt((j-40)^2+(i-30)^2 )< nely/10
%             passive(j,i)=1;
%         end
% %         if sqrt((j-60)^2+(i-30)^2 )< nely/10
% %             passive(j,i)=1;
% %         end
%         
%     end
% end
%passive(70, 1:20)=1;
pass_val=0.01;
x(passive==1) = pass_val;
pore=find(passive==1);

X=x(:);
X(pore) = [];
vol = 0.05;
% vol2=mean(x(:));
%% MMA prepare
m=1;
n=nely*nelx-length(pore);
xold1=X;
xold2=xold1;
% move=0.5;
%% topology optimization loop
obj = zeros(maxiter,1);
xx = zeros(nely*nelx,maxiter);
dd = zeros((nely+1)*(nelx+1),maxiter);
 change = 1;
voliter=2;
Fload=zeros(loadN+6,maxiter); Uload=zeros(loadN+6,maxiter);

for iter = 1:maxiter
    [obj(iter),dc,d,d0,Fload(:,iter),Uload(:,iter)] = crack_propagation(x,iter, elsize, loadN, l);
    %dc = (dc+flipud(dc))/2;
%     dc = sqrt(abs(dc));
    %% filtering and modifying sensitivity numbers
    [H,Hs,rr] = indices_filter(nelx,nely,rmin,voliter,iter); % filtering indices
    dc(:) = H*dc(:)./Hs;
    %         dc=dc(:);
    %         dc(pore) = [];
    
    post_field_plot(nelx,nely,x,d0);
    saveas(gcf,['design_iter_',num2str(iter),'_top.png']);
    post_field_plot(nelx,nely,x,d);
    saveas(gcf,['design_iter_',num2str(iter),'_crack.png']);
    
    close all; pause(1e-6);
    if iter >= voliter; change = abs((obj(iter)-obj(iter-1))/obj(iter-1)); end
    fprintf('it.:%3i rmin.:%4.3f obj.:%9.6f vol.:%4.3f ch.:%7.5f \n',...
        iter,rmin,obj(iter),mean(x(:)),change);
    
    X=x(:);
    X(pore) = [];
    xx(:,iter) = x(:); dd(:,iter) = d(:); % record x-d into xx-dd
    if change <= 1e-6; break; end
     %  [x,olddc] = optimizer_beso(x,dc,vol,passive);
    [x,olddc] = optimizer_simp(x,dc,vol,passive,pass_val);
    %% MMA
%          xold1=reshape(xold1,n,1);
%          xold2=reshape(xold2,n,1);
%          fval = sum(X)-(vol*n);
%          dfdx = ones(1,n);
%          xmin(1:n,1)=0.01;xmax(1:n,1)=1;
%          low(1:n,1)=0;
%          upp=xmax;
%          a0=1;a=0;
%          c0=10000;ddd=0;
%     %
%     %
%          [xmma,y,z,lam,xsi,eta,mu,zet,s,low,upp]=mmasub(m,n,iter,X,xmin,xmax,xold1,xold2, ...
%              obj(iter),dc,fval,dfdx,low,upp,a0,a,c0,ddd);
%          olddc=dc;
%          xold2 = xold1;
%          xold1 = X;
%     %     %
%          xxx=zeros(1,length(xmma)+length(pore));
%          xxx(pore)=1;
%          xxx(~xxx)=xmma;
%          xxx(pore)=0.001;
%     %     %
%          x = reshape(xxx,nely,nelx);
end
 save('results_data.mat', 'obj', 'nelx', 'nely', 'xx', 'dd', 'Uload', 'Fload', 'd', 'elsize', 'rmin', 'l')
