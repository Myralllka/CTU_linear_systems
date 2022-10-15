%% For a linearized model determine the input u(t) that brings the pendulum
% from initial angle 5◦ and zero initial angular velocity to the origin in 0.5
% and 1 second with minimum energy and plot x1(t), x2(t) and u(t) as the
% functions of t as Matlab figures (i.e. 3 figures each with 2 courses, do not
% forget units).

% u(t)        [N]       the torque applied to the pendulum [N]
% x1(t)       [rad]     deviation of the pendulum from the lower position
% x2(t)       [rad/s]   angular velocity

clear all;
close all;
clc;

n_iter = 400;

m = 0.1;    % [kg]      mass
l = 0.3;    % [m]       length
b = 0.05;   % [Nms]     fiction
g = 10;     % [m/s^2]   gravitation const

A = [0,     1; 
    -g/l,   -b/(m*(l^2))];

B = [0; 1/(m*(l^2))];

a = 5;      % [deg]     initial angle
T = 0.5;

system = ss(A, B, [0, 0], [0]);

opt = gramOptions('TimeInterval', [0, T]);
Wc = gram(system, 'c', opt);

x_from = [a*pi/180; 0];
x_to = [0; 0];

u = zeros(1, n_iter);

t = linspace(0, T, n_iter);
for k=1:n_iter
    u(k) = transpose(B) * expm(transpose(A)*(T - t(k))) * inv(Wc) * (x_to - expm(A * T)* x_from);
end

[y,tOut,x] = lsim(system, u, t, x_from);
figure(1)
plot(tOut, x)

grid on
