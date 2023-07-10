%% IC
Ttot = 620;
g = 9.81;
dt = 0.5;
lx = 2000; % domain size
lz = 2000;
Nx = 403; % number of grid points
Nz = 401;
Nt = Ttot/dt+1;
dx = lx/(Nx-3); 
dz = lz/(Nz-1);
r0 = 250;
theta0 = 300;
dtheta = 0.5;
x0 = lx/2;
z0 = 260;
x = linspace (-dx,lx+dx,Nx);
z = linspace (0,lz,Nz);
t = linspace (0,Ttot,Nt);
theta = zeros(Nx,Nz);
u = zeros(Nx,Nz);
w = zeros(Nx,Nz);
v = zeros(Nx,Nz);
tol = 10^-8;
sigma = 1/(1+dx^2/dz^2) * (cos(pi/(Nx-2))+dx^2/dz^2 * cos(pi/Nz));
alpha = 2/(1+sqrt(1-sigma^2));
U = sqrt(2*r0*g*dtheta/theta0);
T = r0/U;
nu = 0.001;
vn = nu*dt/dx^2;
for i = 1:Nx
    for j = 1:Nz
        if ((x(i) - x0)^2 + (z(j) - z0)^2)^(1/2) <= r0
            theta(i,j) = theta0 + dtheta;
        else
            theta(i,j) = theta0;
        end
    end
end
%%
phi = SOR(Nx,Nz,dx,dz,w,tol,alpha);
for i = 2:Nx-1
    for j = 2:Nz-1
        u(i,j) = -(phi(i,j+1) - phi(i,j-1))/(2*dz);
        v(i,j) = (phi(i+1,j) - phi(i-1,j))/(2*dx);
    end
end
u(1,:) = u(Nx-2,:);
u(Nx,:) = u(3,:);
v(1,:) = v(Nx-2,:);
v(Nx,:) = v(3,:);

cmax = max(max(abs(u)*dt/dx + abs(v)*dt/dz));
zeta = zeros(Nt,Nx,Nz);
count = zeros(1,Nt);
wnew = w;
%%
n=1;
J_w = Arakawa(phi,w,dx,dz,Nx,Nz);
J_theta = Arakawa(phi,theta,dx,dz,Nx,Nz);
for i = 2:Nx-1
    for j = 2:Nz-1
        wnew(i,j) = w(i,j) - dt*J_w(i,j) - g/theta0*2*dt/dx*(theta(i+1,j) - theta(i-1,j)); % add diffusion
        wnew(i,j) = wnew(i,j) + nu * (w(i+1,j) - 2*w(i,j) + w(i-1,j)) + nu * (w(i,j+1) - 2*w(i,j) + w(i,j-1));
    end
end
wold = w;
w = wnew;
zeta(n,:,:) = w;
w(1,:) = w(Nx-2,:);
w(Nx,:) = w(3,:);

thetanew = theta;
for i = 2:Nx-1
    for j = 2:Nz-1
        thetanew(i,j) = theta(i,j) - dt*J_theta(i,j);
        thetanew(i,j) = thetanew(i,j) + nu * (theta(i+1,j) - 2*theta(i,j) + theta(i-1,j)) + nu * (theta(i,j+1) - 2*theta(i,j) + theta(i,j-1));
    end
end
thetaold = theta;
theta = thetanew;
temp = zeros(Nt,Nx,Nz);
temp(n,:,:) = theta;
theta(1,:) = theta(Nx-2,:);
theta(Nx,:) = theta(3,:);

[phi,count(n)] = SOR(Nx,Nz,dx,dz,w,tol,alpha);
for i = 2:Nx-1
    for j = 2:Nz-1
        u(i,j) = -(phi(i,j+1) - phi(i,j-1))/(2*dz);
        v(i,j) = (phi(i+1,j) - phi(i-1,j))/(2*dx);
    end
end
u(1,:) = u(Nx-2,:);
u(Nx,:) = u(3,:);
v(1,:) = v(Nx-2,:);
v(Nx,:) = v(3,:);

cmax(n) = max(max(abs(u)*dt/dx + abs(v)*dt/dz));
vmax(n) = max(max(v));
thetamax(n) = max(max(theta));

%%
for n = 2:Nt
    J_w = Arakawa(phi,w,dx,dz,Nx,Nz);
    J_theta = Arakawa(phi,theta,dx,dz,Nx,Nz);
    for i = 2:Nx-1
        for j = 2:Nz-1
            wnew(i,j) = wold(i,j) - 2*dt*J_w(i,j) - g/theta0*2*dt/dx*(theta(i+1,j) - theta(i-1,j));
            wnew(i,j) = wnew(i,j) + 2*nu * (wold(i+1,j) - 2*wold(i,j) + wold(i-1,j)) + 2*nu * (wold(i,j+1) - 2*wold(i,j) + wold(i,j-1));
        end
    end
    wold = w;
    w = wnew;
    zeta(n,:,:) = w;
    w(1,:) = w(Nx-2,:);
    w(Nx,:) = w(3,:);
    
    for i = 2:Nx-1
        for j = 2:Nz-1
            thetanew(i,j) = thetaold(i,j) - 2*dt*J_theta(i,j);
            thetanew(i,j) = thetanew(i,j) + 2*nu * (thetaold(i+1,j) - 2*thetaold(i,j) + thetaold(i-1,j)) + 2*nu * (thetaold(i,j+1) - 2*thetaold(i,j) + thetaold(i,j-1));

        end
    end
    thetaold = theta;
    theta = thetanew;
    temp(n,:,:) = theta;
    theta(1,:) = theta(Nx-2,:);
    theta(Nx,:) = theta(3,:);

    [phi,count(n)] = SOR(Nx,Nz,dx,dz,w,tol,alpha);
    for i = 2:Nx-1
        for j = 2:Nz-1
            u(i,j) = -(phi(i,j+1) - phi(i,j-1))/(2*dz);
            v(i,j) = (phi(i+1,j) - phi(i-1,j))/(2*dx);
        end
    end
    u(1,:) = u(Nx-2,:);
    u(Nx,:) = u(3,:);
    v(1,:) = v(Nx-2,:);
    v(Nx,:) = v(3,:);

    cmax(n) = max(max(abs(u)*dt/dx + abs(v)*dt/dz));
    vmax(n) = max(max(v));
    thetamax(n) = max(max(theta));
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