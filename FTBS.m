function [t,error,amp,u,uexact] = FTBS(dt,LX,save,square)
%   FTBS scheme
%   Allows running FTBS scheme for any wave at specified delta t and
%   wavelength, saving at specified intervals, and allows use of square
%   wave if desired

T = 200; % final time in s
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

uc(1:Nx) = A * sin(2*pi/LX*x); % initial condition
if square == 1 % for square wave
    uc = sign(uc);
end

% allocate arrays for speed
uf = zeros(1,Nx);
error = zeros(1,Nt);
amp = zeros(1,Nt);
u = zeros(floor((Nt-1)/save+1),Nx);
uexact = zeros(floor((Nt-1)/save+1),Nx);

% main loop
for n = 1:Nt
    uf(1) = uc(1)-courant*(uc(1)-uc(Nx-1));
    for i = 2:Nx
        uf(i) = uc(i)-courant*(uc(i)-uc(i-1));
    end
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
    uc = uf;
end
end