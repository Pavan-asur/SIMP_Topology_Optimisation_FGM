function [xnew,dc] = optimizer_simp(x,dc,volfrac,passive,pass_val)
[nely,nelx] = size(x);
l1 = 0; l2 = 100000; move = 0.2;
while (l2-l1 > 1e-13)
  lmid = 0.5*(l2+l1);
  xnew = max(0.001,max(x-move,min(1.,min(x+move,x.*sqrt(-dc./lmid)))));
  xnew(passive==1) = pass_val; % passive region
  if sum(sum(xnew)) - volfrac*nelx*nely > 0
    l1 = lmid;
  else
    l2 = lmid;
  end
end
end