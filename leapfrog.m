function [t,error,amp,u,uexact] = leapfrog(dt,LX,save,square)
%   Leapfrog scheme
%   Allows running Leapfrog scheme for any wave at specified delta t and
%   wavelength, saving at specified intervals, and allows use of square
%   wave if desired

T = 2000; % final time in s
dx = 1;
v = 0.5; % propagation velocity
courant = v*dt/dx;
A = 1; % initial amplitude
l = 50; % domain size
Nx = l/dx+1; % number of grid points
Nt = T/dt+1; % number of time steps
t = linspace(0,T,Nt); % time array
x = linspace(0,l,Nx); % distance array
if ~exist('save','var') % default save interval
    save = 10;
end
if ~exist('square','var') % default is non-square wave
    square = 0;
end

up(1:Nx) = A * sin(2*pi/LX*x); % initial condition
if square == 1 % for square wave
    up = sign(up);
end

% allocate arrays for speed
uf = zeros(1,Nx);
error = zeros(1,Nt);
amp = zeros(1,Nt);
u = zeros(floor((Nt-1)/save+1),Nx);
uexact = zeros(floor((Nt-1)/save+1),Nx);
uc = zeros(1,Nx);

% n=1 requires FTBS
uc(1) = up(1)-courant*(up(1)-up(Nx-1));
for i = 2:Nx
    uc(i) = up(i)-courant*(up(i)-up(i-1));
end
exact = A*sin(2*pi/LX * (x-v*t(1)));
if square==1
    exact=sign(exact);
end
e = abs(exact-up);
error(1)=max(e);
amp(1)=max(up);
u(1,:) = up;
uexact(1,:) = exact;

% main loop
for n = 2:Nt
    uf(1) = up(1)-courant*(uc(2)-uc(Nx-1));
    for i = 2:(Nx-1)
        uf(i) = up(i)-courant*(uc(i+1)-uc(i-1));
    end
    uf(Nx) = up(Nx)-courant*(uc(2)-uc(Nx-1));
    exact = A*sin(2*pi/LX * (x-v*t(n))); % exact solution at time t(n)
    if(square == 1)
        exact = sign(exact);
    end
    e = abs(exact-uc);
    error(n) = max(e);
    amp(n) = max(uc);
    if mod(n-1,save) == 0
        u((n-1)/save+1,:) = uc;
        uexact((n-1)/save+1,:) = exact;
    end
    up = uc;
    uc = uf;
end
end