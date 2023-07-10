%% Vorticity
epsilon = 0;
dx = 2;
dy = 2;
A = 1; % initial amplitude
l = 100; % domain size
Nx = l/dx+3; % number of grid points (including two ghost nodes)
Ny = l/dy+3;
r = 25;
x01 = 50;
y01 = 25;
x02 = 50;
y02 = 75;
x = linspace (-2,l+2,Nx);
y = linspace (-2,l+2,Ny);
vorticity = zeros(Nx,Ny);
tol = 10^-7;
beta = dx/dy;
sigma = 1/(1+beta^2) * (cos(2*pi/(Nx-2))+beta^2 * cos(2*pi/(Ny-2)));
%alpha = 2/(1+sqrt(1-sigma^2));
alpha = 1.99;
for i = 2:Nx-1
    for j = 2:Ny-1
        d1 = min(1,1/r * sqrt( (x(i) - x01).^2 + (y(j) - y01).^2 ) );
        d2 = min(1,1/r * sqrt( (x(i) - x02).^2 + (y(j) - y02).^2 ) );
        vorticity(i,j) = 1/2 * A * (cos(d1*pi)+1) + 1/2 * A * (cos(d2*pi)+1);  % initial condition
    end
end
vorticity(1,:) = vorticity(Nx-2,:); vorticity(:,1) = vorticity(:,Ny-2);
vorticity(Nx-1,:) = vorticity(2,:); vorticity(:,Ny-1) = vorticity(:,2);
vorticity(Nx,:) = vorticity(3,:); vorticity(:,Ny) = vorticity(:,3);

while (mean(mean(vorticity(2:52,2:52)))>tol)
    vorticity = vorticity - mean(mean(vorticity(2:52,2:52)));
end

num2 = max(max(abs(vorticity(2:52,2:52))));
num1 = 2/dx^2 + 2/dy^2;
%% Phi = 0
phi = zeros(Nx,Ny);
R = zeros(Nx,Ny);
for i = 2:Nx-1
    for j = 2:Ny-1
        R(i,j) = 1/dx^2 * (phi(i-1,j) - 2*phi(i,j)+phi(i+1,j)) + 1/dy^2 * (phi(i,j-1) - 2*phi(i,j)+phi(i,j+1)) - vorticity(i,j);
    end
end

R(1,:) = R(Nx-2,:);
R(:,1) = R(:,Ny-2);
R(Nx,:) = R(3,:);
R(:,Nx) = R(:,3);

e = max(max(abs(R(2:52,2:52)))) / ( num1*sum(abs(phi(2:52,2:52)),'all') + num2);
epsilon(1) = e;
%% Main Loop - Run A
count = 0;
while e>tol
    count = count+1;
    for i = 2:Nx-1
        if i == Nx-1
            phi(Nx,:) = phi(3,:);
        end
        for j = 2:Ny-1
            if j == Ny-1
                phi(i,Ny) = phi(i,3);
            end
            R(i,j) = 1/dx^2 * (phi(i-1,j) - 2*phi(i,j)+phi(i+1,j)) + 1/dy^2 * (phi(i,j-1) - 2*phi(i,j)+phi(i,j+1)) - vorticity(i,j);
            phi(i,j) = phi(i,j) + alpha/num1 * R(i,j);
        end
        phi(i,1) = phi(i,Ny-2);
    end
    phi(1,:) = phi(Nx-2,:);
    e = max(max(abs(R(2:52,2:52)))) / ( num1*sum(abs(phi(2:52,2:52)),'all') + num2);
    epsilon(count+1) = e;
end
phi = phi-mean(mean(phi));
%% u and v
u = zeros(Nx,Ny);
v = zeros(Nx,Ny);
for i = 2:Nx-1
    for j = 2:Ny-1
        u(i,j) = -(phi(i,j+1) - phi(i,j-1))/(2*dy);
        v(i,j) = (phi(i+1,j) - phi(i-1,j))/(2*dx);
    end
end
u(1,:) = u(Nx-2,:); u(:,1) = u(:,Ny-2); v(1,:) = v(Nx-2,:); v(:,1) = v(:,Ny-2);
u(Nx,:) = u(3,:); u(:,Ny) = u(:,3); v(Nx,:) = v(3,:); v(:,Ny) = v(:,3);
%%
vort = zeros(Nx,Ny);
for i = 2:Nx-1
    for j = 2:Ny-1
        vort(i,j) = -(u(i,j+1) - u(i,j-1))/(2*dy) + (v(i+1,j) - v(i-1,j))/(2*dx);
    end
end
vort(1,:) = vort(Nx-2,:); vort(:,1) = vort(:,Ny-2);
vort(Nx,:) = vort(3,:); vort(:,Ny) = vort(:,3);