function [] = task6plot(Nt,Nx,Ny,zeta,x,y)
for i = 1:Nt
    zetatemp = zeros(Nx-2,Ny-2);
    zetatemp(1:Nx-2,1:Ny-2) = zeta(i,2:Nx-1,2:Ny-1);
    pcolor(x(2:Nx-1),y(2:Ny-1),zetatemp')
    shading flat
    axis off
    %colorbar    
    drawnow
end
end