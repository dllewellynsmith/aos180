T = 400; % final time in s
dt = 0.4;
dx = 1;
dy = 1;
vx = 1; % propagation velocity
vy = 1;
gammax = 0.1;
gammay = 0.1;
vnx = gammax*dt/dx^2;
vny = gammay*dt/dy^2;
cx = vx*dt/dx;
cy = vy*dt/dy;
A = 1; % initial amplitude
l = 50; % domain size
Nx = l/dx+3; % number of grid points (including two ghost nodes)
Ny = l/dy+3;
Nt = floor(T/dt+1); % number of time steps
r = 12.5;
save = 10;
x0 = 25;
y0 = 25;
x = linspace (-1,l+1,Nx);
y = linspace (-1,l+1,Ny);

uc = zeros(Nx,Ny);
for i = 1:Nx
    for j = 1:Ny
        d = min(1,1/r * sqrt( (x(i) - x0).^2 + (y(j)-y0).^2 ) );
        uc(i,j) = 1/2 * A * (cos(d*pi)+1);  % initial condition
    end
end
%%
% allocate arrays for speed
uf = zeros(Nx,Ny);
error = zeros(1,Nt);
mass = zeros(1,Nt);
streamwise = zeros(Nt,Nx);
crosswise = zeros(Nt,Nx);
u = zeros(floor((Nt-1)/save+1),Nx,Ny);

% n=1 requires FTBS

for j = 2:Ny-2
    for i = 2:Nx-2
        uf(i,j) = uc(i,j) - cx * (uc(i,j) - uc(i-1,j)) - cy * (uc(i,j) - uc(i,j-1));
        uf(i,j) = uf(i,j) + vnx * (uc(i+1,j) - 2*uc(i,j) + uc(i-1,j)) + vny * (uc(i,j+1) - 2*uc(i,j) + uc(i,j-1));
    end
end

uf(1,:) = uf(Nx-2,:);
uf(:,1) = uf(:,Ny-2);
uf(Nx-1,:) = uf(2,:);
uf(:,Nx-1) = uf(:,2);
uf(Nx,:) = uf(3,:);
uf(:,Nx) = uf(:,3);
up = uc;
uc = uf;

umax = zeros(1,Nt);
umax(1)=max(max(up));
u(1,:,:) = up;
ureal = up(2:51,2:51);
mass(1) = sum(sum(ureal));
for i = 1:53
    streamwise(1,i) = up(i,i);
    crosswise(1,i) = up(Ny-i+1,i);
end

% main loop
for n = 2:Nt

    for j = 2:Ny-2
        for i = 2:Nx-2
            uf(i,j) = up(i,j) - cx * (uc(i+1,j) - uc(i-1,j)) - cy * (uc(i,j+1) - uc(i,j-1));
            uf(i,j) = uf(i,j) + 2 * vnx * (up(i+1,j) - 2*up(i,j) + up(i-1,j)) + 2 * vny * (up(i,j+1) - 2*up(i,j) + up(i,j-1));
        end
    end
    uf(1,:) = uf(Nx-2,:);
    uf(:,1) = uf(:,Ny-2);
    uf(Nx-1,:) = uf(2,:);
    uf(:,Nx-1) = uf(:,2);
    uf(Nx,:) = uf(3,:);
    uf(:,Nx) = uf(:,3);
    up = uc;
    uc = uf;

    umax(n) = max(max(up));
    ureal = up(2:51,2:51);
    mass(n) = sum(sum(ureal));

    if mod(n-1,save) == 0
        u((n-1)/save+1,:,:) = up;
    end
    for i = 1:53
        streamwise(n,i) = up(i,i);
        crosswise(n,i) = up(Ny-i+1,i);
    end

end

