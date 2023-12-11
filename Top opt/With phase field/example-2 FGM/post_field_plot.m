function post_field_plot(nelx,nely,x,d)
figure();
imagesc(-x); hold on;
d(d<=0.4) = NaN;
[X,Y] = meshgrid(0.5:1:nelx+0.5,0.5:1:nely+0.5); 
surf(X,Y,1.0*reshape(d,nely+1,nelx+1));
view(0,90); shading interp;
colormap jet; caxis([-1,1.0]);%
axis equal; axis tight; axis off;