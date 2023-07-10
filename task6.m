%% Vorticity
T = 20*24*60*60;
dt = 24*60;
Gamma = 8*10^-5;
beta = 0; %2*10^-11;
lx = 2000*1000; % domain size
ly = 2000*1000;
Nx = 101; % number of grid points
Ny = 101;
Nt = T/dt+1;
dx = lx/(Nx-1);
dy = ly/(Ny-1);
a = 180*1000;
b = 600*1000;
x01 = lx/2+b/2;
y01 = ly/2;
x02 = lx/2-b/2;
y02 = ly/2;
x = linspace (0,lx,Nx);
y = linspace (0,ly,Ny);
t = linspace (0,T,Nt);
vorticity = zeros(Nx,Ny);
tol = 10^-7;
sigma = 1/(1+dx^2/dy^2) * (cos(pi/(Nx-2))+dx^2/dy^2 * cos(pi/(Ny-2)));
alpha = 2/(1+sqrt(1-sigma^2));
for i = 1:Nx
    for j = 1:Ny
        vorticity(i,j) = Gamma * (exp(-1/a^2 * ((x(i)-x01).^2+(y(j)-y01).^2))+exp(-1/a^2 * ((x(i)-x02).^2+(y(j)-y02).^2)));
    end
end
vorticity(1,:) = 0;
vorticity(:,1) = 0;
vorticity(Nx,:) = 0;
vorticity(:,Ny) = 0;
vorticity = vorticity - mean(mean(vorticity));
%%
phi = poisson(Nx,Ny,dx,dy,vorticity,tol,alpha);
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
    end
end
vold = vorticity;
vorticity = vnew;
zeta(n,:,:) = vorticity;
[phi,count(n)] = poisson(Nx,Ny,dx,dy,vorticity,tol,alpha);
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
        end
    end
    vold = vorticity;
    vorticity = vnew;
    zeta(n,:,:) = vorticity;
    % normally would need to deal with boundary conditions here but not
    % necessary in this case
    [phi,count(n)] = poisson(Nx,Ny,dx,dy,vorticity,tol,alpha);
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