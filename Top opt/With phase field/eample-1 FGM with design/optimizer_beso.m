function [x,dc] = optimizer_beso(x,dc,vol,passive) 
l1 = min(min(dc)); l2 = max(max(dc));
xold = x; move = 0.02;
% determine alpha_th
while ((l2-l1)/l2 > 1.0e-5)
   th = (l1+l2)/2.0;
   x = max(0,sign(dc-th));
   x(passive==1) = 0; % passive region
   if mean(x(:)) - vol > 0
      l1 = th;
   else
      l2 = th;
   end
end
% determine alpha_th^add and alpha_th^del
if sum(x(xold==0))/sum(xold(:)) > move
   l1 = min(min(dc)); l2 = max(max(dc));
   while ((l2-l1)/l2 > 1.0e-5)
      thadd = (l1+l2)/2.0;
      x(xold==0) = max(0,sign(dc(xold==0)-thadd));
      x(passive==1) = 0; % passive region
      if sum(x(xold==0))/sum(xold(:)) > move
         l1 = thadd;
      else
         l2 = thadd;
      end
   end
   l1 = min(min(dc)); l2 = max(max(dc));
   while ((l2-l1)/l2 > 1.0e-5)
      thdel = (l1+l2)/2.0;
      x(xold==1) = max(0,sign(dc(xold==1)-thdel));
      x(passive==1) = 0; % passive region
      if mean(x(:)) - vol > 0
         l1 = thdel;
      else
         l2 = thdel;
      end
   end
end