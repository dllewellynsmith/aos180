%% Vorticity
T = 30*24*60*60;
dt = 12*60;
Gamma = 8*10^-5;
beta = 0; %2*10^-11;
lx = 2000*1000; % domain size
ly = 2000*1000;
Nx = 201; % number of grid points
Ny = 201;
Nt = T/dt+1;
dx = lx/(Nx-1);
dy = ly/(Ny-1);
a = 180*1000;
b = 300*1000;
nu = 0.01;
kappa = nu;
x01 = lx/2+b/2;
y01 = ly/2;
x02 = lx/2-b/2;
y02 = ly/2;
%x03 = lx/2;
%y03 = 2*ly/3;
x = linspace (0,lx,Nx);
y = linspace (0,ly,Ny);
t = linspace (0,T,Nt);
vorticity = zeros(Nx,Ny);
tol = 10^-10;
sigma = 1/(1+dx^2/dy^2) * (cos(pi/(Nx-2))+dx^2/dy^2 * cos(pi/(Ny-2)));
alpha = 2/(1+sqrt(1-sigma^2));
for i = 2:Nx-1
    for j = 2:Ny-1
        vorticity(i,j) = Gamma * (exp(-1/a^2 * ((x(i)-x01).^2+(y(j)-y01).^2))-exp(-1/a^2 * ((x(i)-x02).^2+(y(j)-y02).^2)));%+exp(-1/a^2 * ((x(i)-x03).^2+(y(j)-y03).^2)));
    end
end
vorticity = vorticity - mean(mean(vorticity));
S1 = zeros(Nx,Ny);
for i = 2:Nx-1
    for j = 2:11
            S1(i,j) = 1;
    end
end
%%
phi = poisson(Nx,Ny,dx,dy,vorticity,tol,alpha);
vorticity(1,:) = (6*phi(2,:)-3/2*phi(3,:)+2/9*phi(4,:))/(dx^2);
vorticity(:,1) = (6*phi(:,2)-3/2*phi(:,3)+2/9*phi(:,4))/(dx^2);
vorticity(Nx,:) = (6*phi(Nx-1,:)-3/2*phi(Nx-2,:)+2/9*phi(Nx-3,:))/(dx^2);
vorticity(:,Ny) = (6*phi(:,Ny-1)-3/2*phi(:,Ny-2)+2/9*phi(:,Ny-3))/(dx^2);
u = zeros(Nx,Ny);
v = zeros(Nx,Ny);
for i = 2:Nx-1
    for j = 2:Ny-1
        u(i,j) = -(phi(i,j+1) - phi(i,j-1))/(2*dy);
        v(i,j) = (phi(i+1,j) - phi(i-1,j))/(2*dx);
    end
end
K = sum(sum(1/2*(u.^2 + v.^2)));
omega = sum(sum(vorticity.^2));
cmax = max(max(abs(u)*dt/dx + abs(v)*dt/dy));
zeta = zeros(Nt,Nx,Ny);
count = zeros(1,Nt);
vnew = vorticity;
%%
n=1;
J_a = Arakawa(phi,vorticity,dx,dy,Nx,Ny);
for i = 2:Nx-1
    for j = 2:Ny-1
        vnew(i,j) = vorticity(i,j) - dt*J_a(i,j) - beta*2*dt/dx*(phi(i+1,j) - phi(i-1,j));
        vnew(i,j) = vnew(i,j) + nu * (vorticity(i+1,j) - 2*vorticity(i,j) + vorticity(i-1,j)) + nu * (vorticity(i,j+1) - 2*vorticity(i,j) + vorticity(i,j-1));
    end
end
vold = vorticity;
vorticity = vnew;
zeta(n,:,:) = vorticity;

J_S1 = Arakawa(phi,S1,dx,dy,Nx,Ny);
S1new = S1;
S1old = S1;
for i = 2:Nx-1
    for j = 2:Ny-1
        S1new(i,j) = S1(i,j) - dt*J_S1(i,j);
        S1new(i,j) = S1new(i,j) + 2*kappa * (S1old(i+1,j) - 2*S1old(i,j) + S1old(i-1,j)) + 2*kappa * (S1old(i,j+1) - 2*S1old(i,j) + S1old(i,j-1));
    end
end
S1old = S1;
S1 = S1new;
scalar1 = zeros(Nt,Nx,Ny);
scalar1(n,:,:) = S1;

[phi,count(n)] = poisson(Nx,Ny,dx,dy,vorticity,tol,alpha);
vorticity(1,:) = (6*phi(2,:)-3/2*phi(3,:)+2/9*phi(4,:))/(dx^2);
vorticity(:,1) = (6*phi(:,2)-3/2*phi(:,3)+2/9*phi(:,4))/(dx^2);
vorticity(Nx,:) = (6*phi(Nx-1,:)-3/2*phi(Nx-2,:)+2/9*phi(Nx-3,:))/(dx^2);
vorticity(:,Ny) = (6*phi(:,Ny-1)-3/2*phi(:,Ny-2)+2/9*phi(:,Ny-3))/(dx^2);
for i = 2:Nx-1
    for j = 2:Ny-1
        u(i,j) = -(phi(i,j+1) - phi(i,j-1))/(2*dy);
        v(i,j) = (phi(i+1,j) - phi(i-1,j))/(2*dx);
    end
end
K(n) = sum(sum(1/2*(u.^2 + v.^2)));
omega(n) = sum(sum(vorticity.^2));
cmax(n) = max(max(abs(u)*dt/dx + abs(v)*dt/dy));
%%
for n = 2:Nt
    J_a = Arakawa(phi,vorticity,dx,dy,Nx,Ny);
    for i = 2:Nx-1
        for j = 2:Ny-1
            vnew(i,j) = vold(i,j) - 2*dt*J_a(i,j) - beta*2*dt/dx*(phi(i+1,j) - phi(i-1,j));
            vnew(i,j) = vnew(i,j) + 2*nu * (vold(i+1,j) - 2*vold(i,j) + vold(i-1,j)) + 2*nu * (vold(i,j+1) - 2*vold(i,j) + vold(i,j-1));
        end
    end
    vold = vorticity;
    vorticity = vnew;
    zeta(n,:,:) = vorticity;

    J_S1 = Arakawa(phi,S1,dx,dy,Nx,Ny);
    S1new = S1;
    for i = 2:Nx-1
        for j = 2:Ny-1
            S1new(i,j) = S1old(i,j) - 2*dt*J_S1(i,j);
            S1new(i,j) = S1new(i,j) + 2*kappa * (S1old(i+1,j) - 2*S1old(i,j) + S1old(i-1,j)) + 2*kappa * (S1old(i,j+1) - 2*S1old(i,j) + S1old(i,j-1));
        end
    end
    S1old = S1;
    S1 = S1new;
    scalar1(n,:,:) = S1;

    [phi,count(n)] = poisson(Nx,Ny,dx,dy,vorticity,tol,alpha);
    vorticity(1,:) = (6*phi(2,:)-3/2*phi(3,:)+2/9*phi(4,:))/(dx^2);
    vorticity(:,1) = (6*phi(:,2)-3/2*phi(:,3)+2/9*phi(:,4))/(dx^2);
    vorticity(Nx,:) = (6*phi(Nx-1,:)-3/2*phi(Nx-2,:)+2/9*phi(Nx-3,:))/(dx^2);
    vorticity(:,Ny) = (6*phi(:,Ny-1)-3/2*phi(:,Ny-2)+2/9*phi(:,Ny-3))/(dx^2);
    for i = 2:Nx-1
        for j = 2:Ny-1
            u(i,j) = -(phi(i,j+1) - phi(i,j-1))/(2*dy);
            v(i,j) = (phi(i+1,j) - phi(i-1,j))/(2*dx);
        end
    end
    K(n) = sum(sum(1/2*(u.^2 + v.^2)));
    omega(n) = sum(sum(vorticity.^2));
    cmax(n) = max(max(abs(u)*dt/dx + abs(v)*dt/dy));
end
%%

function [J] = Jacobian(a,b,dx,dy,Nx,Ny)
J = zeros(Nx,Ny);
for i = 2:Nx-1
    for j = 2:Ny-1
        J(i,j) = 1/(2*dx) * (a(i+1,j) - a(i-1,j)) * 1/(2*dy) * (b(i,j+1)-b(i,j-1)) - 1/(2*dy) * (a(i,j+1) - ...
            a(i,j-1)) * 1/(2*dx) * (b(i+1,j)-b(i-1,j));
    end
end
end

function [Jhat] = JacobianHat(a,b,dx,dy,Nx,Ny)
Jhat = zeros(Nx,Ny);
for i = 2:Nx-1
    for j = 2:Ny-1
        Jhat(i,j) = 1/(4*dx*dy) * (a(i+1,j)*(b(i+1,j+1) - b(i+1,j-1)) - a(i-1,j)*(b(i-1,j+1) - b(i-1,j-1)) - a(i,j+1)*(b(i+1,j+1) ...
            - b(i-1,j+1)) + a(i,j-1)*(b(i+1,j-1) - b(i-1,j-1)) );
    end
end

end

function [Jtilda] = JacobianTilda(a,b,dx,dy,Nx,Ny)
Jtilda = zeros(Nx,Ny);
for i = 2:Nx-1
    for j = 2:Ny-1
        Jtilda(i,j) = 1/(4*dx*dy) * (b(i,j+1)*(a(i+1,j+1) - a(i-1,j+1)) - b(i,j-1)*(a(i+1,j-1) - a(i-1,j-1)) - b(i+1,j) ...
            *(a(i+1,j+1) - a(i+1,j-1)) + b(i-1,j)*(a(i-1,j+1) - a(i-1,j-1)) );
    end
end
end

function[J_a] = Arakawa(a,b,dx,dy,Nx,Ny)
J_a = 1/3*(Jacobian(a,b,dx,dy,Nx,Ny) + JacobianHat(a,b,dx,dy,Nx,Ny) + JacobianTilda(a,b,dx,dy,Nx,Ny));
end