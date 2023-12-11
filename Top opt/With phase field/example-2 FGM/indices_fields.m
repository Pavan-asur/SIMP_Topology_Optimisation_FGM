function [iK,jK,iKd,jKd,iFd,jFd,edofMat,dedofMat] = indices_fields(nelx,nely)

% nodes number
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);

% indices for the displacement field
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);

% indices for the crack phase field
dedofVec = reshape(nodenrs(2:end,1:end-1),nelx*nely,1);
dedofMat = repmat(dedofVec,1,4)+repmat([0 nely+[1 0] -1],nelx*nely,1);
iKd = reshape(kron(dedofMat,ones(4,1))',16*nelx*nely,1);
jKd = reshape(kron(dedofMat,ones(1,4))',16*nelx*nely,1);
iFd = reshape(dedofMat',4*nelx*nely,1);
jFd = ones(4*nelx*nely,1);

