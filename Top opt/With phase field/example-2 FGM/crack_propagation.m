function [obj,dc,d,d0,Fload,Uload] = crack_propagation(x,iter, elsize, loadN, l)
%% sensitivity number
[nely,nelx] = size(x); x = x(:);
dc = zeros(nely*nelx,1); obj = 0;
%% phase field parameters
p = 2; k = 1e-9; % small positive value d = 1
q=3;

E_1=10.4;
E_2=5.2;
E_inc=52;

nu_1=0.3;
nu_2=0.3;
nu_inc=0.3; 

sigmac_1=0.01;
sigmac_2=0.005;
sigmac_inc=0.03;

Vf=repmat(0, nely, nelx); 
power=1;

for i=1:nelx
    Vf(:, i)=i; % get the coords
end

Vf(:,:)=Vf(:,:)./max(max(Vf));
Vf=Vf(:);
%% material properties (Nguyen_CMAME_2015,Cement paste and Sand)

% Young's modulus (GPa - KN/mm2)

E1=repmat(E_1, nely, nelx); E1=E1(:);
E2=repmat(E_2, nely, nelx); E2=E2(:);

E_mat=E1+(E2-E1).*(Vf).^power; 
Einc=repmat(E_inc, nely, nelx); Einc=Einc(:);

%Poisons ratio
nu1=repmat(nu_1, nely, nelx); nu1=nu1(:);
nu2=repmat(nu_2, nely, nelx); nu2=nu2(:);

nu_mat=nu1+(nu2-nu1).*(Vf).^power;
nuinc=repmat(nu_inc, nely, nelx); nuinc=nuinc(:);

%Sigma_c for matrix and inclusion GPa 0.05/0.01

sigma_c1=repmat(sigmac_1, nely, nelx); sigma_c1=sigma_c1(:);
sigma_c2=repmat(sigmac_2, nely, nelx); sigma_c2=sigma_c2(:);

sigma_c_mat=sigma_c1+(sigma_c2-sigma_c1).*(Vf).^power; 
sigma_cinc=repmat(sigmac_inc, nely, nelx); sigma_cinc=sigma_cinc(:);

%Compute Lame constants
lambda_mat=(E_mat.*nu_mat)./((1+nu_mat).*(1-2*nu_mat)); 
lambda_inc=(Einc.*nuinc)./((1+nuinc).*(1-2*nuinc)); 
mu_mat = E_mat./(2*(1+nu_mat)); 
mu_inc = Einc./(2*(1+nuinc));

Psi_c_mat = (sigma_c_mat.^2)./E_mat/2;  % specific fracture energy density (GPa)
Psi_cinc = (sigma_cinc.^2)./Einc/2;  % specific fracture energy density (GPa)

%Compute the properties at guass point based on the SIMP
lambda_gauss = kron(lambda_inc.*x.^q+lambda_mat.*(1-x.^q),ones(4,1))';
mu_gauss = kron(mu_inc.*x.^q+mu_mat.*(1-x.^q),ones(4,1))';
Psi_c_gauss = kron(Psi_cinc.*x.^q+Psi_c_mat.*(1-x.^q),ones(4,1))';
%% boundary conditions of both displacement and crack phase field
[d,dfixed,dfree,u,uload,uload1,uload2,ufree,f] = boundary_conditions(nelx,nely);
u_inc = [2e-3*ones(5,1);0.5e-3*ones(loadN,1)];
%% preprocessing, shape funcitons and indices for the three fields
[Nd,Ndp,Bdp,Bu,J] = shape_functions(elsize);
[iK,jK,iKd,jKd,iFd,jFd,edofMat,dedofMat] = indices_fields(nelx,nely);
%% Staggered solution procedure
Fload = zeros(length(u_inc)+1,1); Uload = zeros(length(u_inc)+1,1);
Hn = zeros(4,nely*nelx); P_plus = zeros(9,4*nely*nelx); P_minus = P_plus;
sK_sen = zeros(64,nely*nelx);
for ij = 1:length(u_inc)
  
    u_old = u; sK_sen_old = sK_sen; f_old = f;
    %% evaluate historical tensorial energy
    [Hn,epsilon_gauss] = Hn_evaluation(nelx,nely,u,Hn,Bu,edofMat,lambda_gauss,mu_gauss,Psi_c_gauss);
    %% crack phase filed solution
    Bterm = reshape(2*Psi_c_gauss*l^2,4,nely*nelx); % B terms of Miehe_CMAME_2015
    Nterm = reshape(2*Psi_c_gauss,4,nely*nelx)+2*Hn; % N terms of Miehe_CMAME_2015
    sKd_B = reshape(Bdp*Bterm,16*nelx*nely,1); % stiffness of B terms
    sKd_D = reshape(Ndp*Nterm,16*nelx*nely,1); % stiffness of D terms
    sKd = 0.25*J*(sKd_D+sKd_B); % ensemble of B and D terms
    sFd = 0.25*J*reshape(2*Nd'*Hn,4*nelx*nely,1); % force terms
    Kd = sparse(iKd,jKd,sKd); %Kd= (Kd+Kd')/2; % global sparse matrix assembly
    Fd = sparse(iFd,jFd,sFd); % global sparce force vector

    d(dfree) = Kd(dfree,dfree)\(Fd(dfree)-Kd(dfree,dfixed)*d(dfixed));
    if ij == 1; d0 = d; end
%                    figure(10)
% %       %  imagesc(x)
%       hold on 
%       [X,Y] = meshgrid(0.5:1:nelx+0.5,0.5:1:nely+0.5); 
%           surf(X,Y,1*reshape(d,nely+1,nelx+1))
%       view(0,90); shading interp;
%       colormap jet;  caxis([0.4,1]);caxis([0,1]);
% %     post_field_plot(nelx,nely,reshape(x,nely,nelx),d);
% %     saveas(gcf,['step',num2str(ij),'_crack.png']);
%     close all;
    %% displacement field solution
    d_gauss = reshape(Nd*d(dedofMat)',1,4*nely*nelx);
    [R_plus,R_minus] = operator_R(epsilon_gauss); % operator R+ R-
    for i = 1:4*nely*nelx % parfor
        [P_plus(:,i),P_minus(:,i)] = operator_P(epsilon_gauss(:,i)); % operator P+ P-
    end
    I = [1;1;0]; Ip = reshape(I*I',9,1);
    D = kron(((1-d_gauss).^p+k),ones(9,1)).*(kron(lambda_gauss.*R_plus,Ip)+2*kron(mu_gauss,repmat([1;1;0.5],3,1)).*P_plus)...
        + kron(lambda_gauss.*R_minus,Ip) + 2*kron(mu_gauss,repmat([1;1;0.5],3,1)).*P_minus;
    sK = zeros(64,nely*nelx);
    for i = 1:nely*nelx % parfor
        sK(:,i) = reshape((0.25*J*Bu'*blkdiag(reshape(D(:,4*(i-1)+1),3,3),reshape(D(:,4*(i-1)+2),3,3),...
            reshape(D(:,4*(i-1)+3),3,3),reshape(D(:,4*(i-1)+4),3,3))*Bu),64,1);
    end
    K = sparse(iK,jK,sK(:));% K = (K+K')/2; % global sparse matrix assembly
    u(uload1) = sum(u_inc(1:ij)); u(uload2) = -sum(u_inc(1:ij));
    u(ufree) = -K(ufree,ufree)\(K(ufree,uload)*u(uload));
    f = K*u;
    Fload(ij+1) = sum(f(uload1));
    Uload(ij+1) = sum(u_inc(1:ij));
    %% sensitivity evaluation
    du=u - u_old;

    lambda_sen = kron((lambda_inc-lambda_mat)*q.*x.^(q-1),ones(4,1))';
    mu_sen = kron((mu_inc-mu_mat)*q.*x.^(q-1),ones(4,1))';
%     D_sen = kron(((1-d_gauss).^p+k),ones(9,1)).*(kron(lambda_sen.*R_plus,Ip)+2*kron(mu_gauss,repmat([1;1;0.5],3,1)).*P_plus)...
%         + kron(lambda_sen.*R_minus,Ip) + 2*kron(mu_gauss,repmat([1;1;0.5],3,1)).*P_minus;
    D_sen = kron(((1-d_gauss).^p+k),ones(9,1)).*(kron(lambda_sen.*R_plus,Ip)+2*kron(mu_sen,repmat([1;1;0.5],3,1)).*P_plus)...
        + kron(lambda_sen.*R_minus,Ip) + 2*kron(mu_sen,repmat([1;1;0.5],3,1)).*P_minus;
    sK_sen = zeros(64,nely*nelx);
    for i = 1:nely*nelx % parfor
        sK_sen(:,i) = reshape((0.25*J*Bu'*blkdiag(reshape(D_sen(:,4*(i-1)+1),3,3),reshape(D_sen(:,4*(i-1)+2),3,3),...
            reshape(D_sen(:,4*(i-1)+3),3,3),reshape(D_sen(:,4*(i-1)+4),3,3))*Bu),64,1);
        dc(i) = dc(i)-1/2*(u_old(edofMat(i,:))'*reshape(sK_sen_old(:,i),8,8)*du(edofMat(i,:))+...
            u(edofMat(i,:))'*reshape(sK_sen(:,i),8,8)*du(edofMat(i,:)));
    end
    obj = obj + 1/2*(f_old+f)'*(u - u_old);
    if ij>10
    if abs(Fload(ij+1)) <= 0.05; break; end
    end
end
%% field plothold off
figure(2)
dc = reshape(dc,nely,nelx);
plot(Uload(1:ij+1),Fload(1:ij+1),'-r','LineWidth',2); axis tight;
tx = xlabel('displacement [mm]'); ty = ylabel('load [KN]');
tt = title('load-deflection curve'); pause(1e-6);
set(tt,'Fontsize',24); set(tx,'Fontsize',24); set(ty,'Fontsize',24);
set(gca,'FontName','Times','FontSize',20) ;
saveas(gcf,['design_iter_',num2str(iter),'_curve.png']);
