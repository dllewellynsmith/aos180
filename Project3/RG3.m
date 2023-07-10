function [t,error,amp,u,uexact] = RG3(dt,LX,save,square)
%   Runge-Kutta order 3 scheme
%   Allows running Runge-Kutta scheme for any wave at specified delta t and
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
temp1 = zeros(1,Nx);
temp2 = zeros(1,Nx);

% main loop
for n = 1:Nt

    temp1(1) = uc(1)-courant/6*(uc(2)-uc(Nx-1));
    temp1(Nx) = temp1(1);
    for i = 2:(Nx-1)
        temp1(i) = uc(i)-courant/6*(uc(i+1)-uc(i-1));
    end

    temp2(1) = uc(1)-courant/4*(temp1(2)-temp1(Nx-1));
    temp2(Nx) = temp2(1);
    for i = 2:(Nx-1)
        temp2(i) = uc(i)-courant/4*(temp1(i+1)-temp1(i-1));
    end

    uf(1) = uc(1)-courant/2*(temp2(2)-temp2(Nx-1));
    uf(Nx) = uf(1);
    for i = 2:(Nx-1)
        uf(i) = uc(i)-courant/2*(temp2(i+1)-temp2(i-1));
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