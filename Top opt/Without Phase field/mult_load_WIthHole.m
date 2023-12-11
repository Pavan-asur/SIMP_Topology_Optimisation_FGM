%%%% AN 88 LINE TOPOLOGY OPTIMIZATION CODE Nov, 2010 %%%%
nelx=60; nely=60; volfrac=0.4; penal=3;  rmin=6; ft=2; 
%% MATERIAL PROPERTIES FGM
nu = 0.3;
E0=zeros(nely, nelx); % for the FGM function
Emin=repmat(0.001, nely, nelx); % as a dummy Stiffness
E1=1; % Define the Youngs modulus of material 1
E2=1; % define the Youngs modulus of Material2
VF_alp=1; % Power for the volume fraction
%Define the properties of the FGM material via defining volume fraction
for elx=1:nelx
    for ely=1:nely 
        %E1 is towards the left and E2 is assigned to the right
        %Uncomment for the above variation
        E0(ely, elx)=E1+(E2-E1)*(((elx-1)/nelx)^VF_alp); % define the volume fraction[(ely-1)/nely)         
        %E1 is assignes to the top and E2 is assigned to the bottom 
        %Uncomment for the above variation
       % E0(ely, elx)=E1+(E2-E1)*(((ely-1)/nely)^VF_alp); % define the volume fraction[(ely-1)/nely)]
    end   
end
E0(:, nelx)=E2; % for the first condition
%E0(nely, :)=E2; % for the second condition 
% adjusting a bit for the volume fraction (otherwise, the last element wil be little less than E2)

%% PREPARE FINITE ELEMENT ANALYSIS
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]); % note than nu here is kept constant
%%%%%%This nu value is constant only under the assumption that E is much
%%%%%%higher than the nu value. If E value is close to 10, nu will have the
%%%%%%effects and that needs to be added as a matrix as well. 
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
% DEFINE LOADS AND SUPPORTS (cantilever beam)
F = sparse([2*(nely+1)*nelx+2, 2*(nely+1)*(nelx+1)], [1,2], [1, -1], 2*(nely+1)*(nelx+1), 2); 
U = zeros(2*(nely+1)*(nelx+1),2);
fixeddofs =[1:2*nely+1]; 
alldofs = [1:2*(nely+1)*(nelx+1)];
freedofs = setdiff(alldofs,fixeddofs);
%% PREPARE FILTER
iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for i1 = 1:nelx
  for j1 = 1:nely
    e1 = (i1-1)*nely+j1;
    for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
      for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
        e2 = (i2-1)*nely+j2;
        k = k+1;
        iH(k) = e1;
        jH(k) = e2;
        sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
      end
    end
  end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);
%% INITIALIZE ITERATION
x = repmat(volfrac,nely,nelx);
xPhys = x;
loop = 0;
change = 1;
passive=zeros(nely, nelx); 
for i=1:nelx
    for j=1:nely
        if sqrt((j-nely/2)^2+(i-nelx/3)^2 )< nely/3
            passive(j,i)=1;
        end
    end
end
%% START ITERATION
while change > 0.01
  loop = loop + 1;
  %% FE-ANALYSIS
  %sK = reshape(KE(:)*(xPhys(:)'.^penal*E0+(1-xPhys(:)'.^penal)*Emin),64*nelx*nely,1);
  sK = reshape(KE(:)*(xPhys(:)'.^penal.*E0(:)'+(1-xPhys(:)'.^penal).*Emin(:)'),64*nelx*nely,1);
  K = sparse(iK,jK,sK); K = (K+K')/2;
  U(freedofs, :) = K(freedofs,freedofs)\F(freedofs,:);  
  %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
  c=0; %chnaged here
  dc=0; %changed here
  for i=1:size(F,2)
      Ui=U(:,i);
  ce = reshape(sum((Ui(edofMat)*KE).*Ui(edofMat),2),nely,nelx);
  c = c+ sum(sum(( Emin+xPhys.^penal.*(E0-Emin)).*ce));
  dc = dc-penal.*(E0-Emin).*xPhys.^(penal-1).*ce;
  end
  dv = ones(nely,nelx);  
  %% FILTERING/MODIFICATION OF SENSITIVITIES
  if ft == 1
    dc(:) = H*(x(:).*dc(:))./Hs./max(1e-3,x(:));
  elseif ft == 2
    dc(:) = H*(dc(:)./Hs);
    dv(:) = H*(dv(:)./Hs);
  end
  %% OPTIMALITY CRITERIA UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES
  l1 = 0; l2 = 1e9; move = 0.2;
  while (l2-l1)/(l1+l2) > 1e-3
    lmid = 0.5*(l2+l1);
    xnew = max(0,max(x-move,min(1,min(x+move,x.*sqrt(-dc./dv/lmid)))));
    if ft == 1
      xPhys = xnew;
    elseif ft == 2
      xPhys(:) = (H*xnew(:))./Hs;
    end
    xPhys(passive==1)=0;
    xPhys(passive==2)=1;
    if sum(xPhys(:)) > volfrac*nelx*nely, l1 = lmid; else l2 = lmid; end
  end
  change = max(abs(xnew(:)-x(:)));
  x = xnew;
  %% PRINT RESULTS
  fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,c, ...
    mean(xPhys(:)),change);
  %% PLOT DENSITIES
  colormap(winter); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;
end

disp('The Optimisation is done')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code is intended for educational purposes and theoretical details    %
% are discussed in the paper                                               %
% "SIMP Phase field topology Optimisation framework to maximise fracture   %
%resistance in FGMs", Pavan Kumar Asur, Pengfie Li, Jose Reinoso, Qi Chang %
%He, Julien Yvonnet, Marco Paggi, In Composites Structure.                 %
%                                                                          %
%Please cite the above paper if you use the codes for research papers      %
%                                                                          % 
%                                                                          %
%Contact the first author Pavan Kumar Asur, for Comments, and suggestion   %
% at pavan.kumar@ilsb.tuwien.ac.at                                         %
%                                                                          %
%                                                                          % 
%                                                                          %
% This Matlab code was based on the codes from                             %
% "Efficient topology optimization in MATLAB using 88 lines of code,       %
% E. Andreassen, A. Clausen, M. Schevenels,                                %
% B. S. Lazarov and O. Sigmund, Struct Multidisc Optim, 2010               %
%                                                                          %
%                                                                          %
%                                                                          %
% Disclaimer:                                                              %
% The authors reserves all rights but do not guaranty that the code is     %
% free from errors. Furthermore, we shall not be liable in any event       %
% caused by the use of the program.                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%