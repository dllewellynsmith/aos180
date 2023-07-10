%% Declare Constants
g = 9.81; % acceleration due to gravity in m/s^2
theta = 290; % temperature in K
dtheta_dz = 5/1000; % change in temperature as height increases in K/m
z1 = 100; % initial height in m
w1 = 0; % initial velocity in m/s
T = 2000; % final time in s
N = sqrt((g/theta)*dtheta_dz); % Brunt-Väisäla frequency in Hz
%% EF delta_t = 10
delta_t = 10; % time step in s

t10 = linspace(0,T,T/delta_t+1); % time array with dt = 10
z_exact_10 = z1*cos(N*t10); % exact solution for height
w_exact_10 = -z1*N*sin(N*t10); % exact solution for velocity
n_total = T/delta_t; % total points

w_approx_10(1) = w1; % set I.C.
z_approx_10(1) = z1;
for i = 1:n_total
    w_approx_10(i+1) = w_approx_10(i) + delta_t*(-N^2)*z_approx_10(i);
    z_approx_10(i+1) = z_approx_10(i) + delta_t*w_approx_10(i);
end

epsilon_z_10 = abs(z_exact_10 - z_approx_10); % error
epsilon_w_10 = abs(w_exact_10 - w_approx_10);

%% EF delta_t = 1
delta_t = 1; % time step in s

t1 = linspace(0,T,T/delta_t+1);
z_exact_1 = z1*cos(N*t1); % exact solution for height
w_exact_1 = -z1*N*sin(N*t1);
n_total = T/delta_t;

w_approx_1(1) = w1;
z_approx_1(1) = z1;
for i = 1:n_total
    w_approx_1(i+1) = w_approx_1(i) + delta_t*(-N^2)*z_approx_1(i);
    z_approx_1(i+1) = z_approx_1(i) + delta_t*w_approx_1(i);
end

epsilon_z_1 = abs(z_exact_1 - z_approx_1);
epsilon_w_1 = abs(w_exact_1 - w_approx_1);
%% Plot 1 for EF
figure(1)
clf
hold on
plot(t10,z_exact_10, 'LineWidth', 2.0)
plot(t10,z_approx_10, 'LineWidth', 2.0)
plot(t1,z_approx_1, 'LineWidth', 2.0)
legend('Exact Solution', 'Euler Forward for dt = 10', 'Euler Forward for dt = 1', 'Location', 'northwest')
title('Vertical position using Euler Forward scheme')
xlabel('Time (s)')
ylabel('Height (m)')
%% Plot 2 for EF
figure(2);
clf
hold on
plot(z_exact_10,w_exact_10, 'LineWidth', 2.0)
plot(z_approx_10,w_approx_10)
plot(z_approx_1,w_approx_1)
legend('Exact Solution', 'Euler Forward for dt = 10', 'Euler Forward for dt = 1')
title('Phase Space using Euler Forward scheme')
xlabel('Height (m)')
ylabel('Velocity (m/s)')
%% Plot 3 for EF
figure(3);
clf
loglog(t10,epsilon_z_10, 'LineWidth', 2.0)
hold on
loglog(t1,epsilon_z_1, 'LineWidth', 2.0)
legend('Euler Forward for dt = 10', 'Euler Forward for dt = 1', 'Location', 'northwest')
title('Error in vertical position using Euler Forward scheme')
xlabel('Time (s)')
ylabel('Error (m)')
%% Leapfrog delta = 10
delta_t = 10;
n_total = T/delta_t;

w_approx_10(1) = w1;
z_approx_10(1) = z1;
w_approx_10(2) = w1 + delta_t*(-N^2)*z_approx_10(1);
z_approx_10(2) = z1 + delta_t*w_approx_10(1);
for i = 2:n_total
    w_approx_10(i+1) = w_approx_10(i-1) + 2*delta_t*(-N^2)*z_approx_10(i);
    z_approx_10(i+1) = z_approx_10(i-1) + 2*delta_t*w_approx_10(i);
end

epsilon_z_10 = abs(z_exact_10 - z_approx_10);
epsilon_w_10 = abs(w_exact_10 - w_approx_10);
%% Leapfrog delta = 1
delta_t = 1;
n_total = T/delta_t;

w_approx_1(1) = w1;
z_approx_1(1) = z1;
w_approx_1(2) = w1 + delta_t*(-N^2)*z_approx_1(1);
z_approx_1(2) = z1 + delta_t*w_approx_1(1);
for i = 2:n_total
    w_approx_1(i+1) = w_approx_1(i-1) + 2*delta_t*(-N^2)*z_approx_1(i);
    z_approx_1(i+1) = z_approx_1(i-1) + 2*delta_t*w_approx_1(i);
end

epsilon_z_1 = abs(z_exact_1 - z_approx_1);
epsilon_w_1 = abs(w_exact_1 - w_approx_1);
%% Plot 1 for Leapfrog
figure(4)
clf
hold on
plot(t10,z_exact_10)
plot(t10,z_approx_10)
plot(t1,z_approx_1)
legend('Exact Solution', 'Leapfrog scheme for dt = 10', 'Leapfrog scheme for dt = 1')
title('Vertical position using Leapfrog scheme')
xlabel('Time (s)')
ylabel('Height (m)')
%% Plot 2 for EF
figure(5);
clf
hold on
plot(z_exact_10,w_exact_10)
plot(z_approx_10,w_approx_10)
plot(z_approx_1,w_approx_1)
legend('Exact Solution', 'Leapfrog scheme for dt = 10', 'Leapfrog scheme for dt = 1')
title('Phase Space using Leapfrog scheme')
xlabel('Height (m)')
ylabel('Velocity (m/s)')
%% Plot 3 for EF
figure(6);
clf
loglog(t10,epsilon_z_10)
hold on
loglog(t1,epsilon_z_1)
legend('Leapfrog scheme for dt = 10', 'Leapfrog scheme for dt = 1', 'Location', 'northwest')
title('Error in vertical position using Leapfrog scheme')
xlabel('Time (s)')
ylabel('Error (m)')