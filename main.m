clear, clc

%%
% u = [theta; d];
a = 0.1; % OA = 0.1 m
b = 0.2; % AB = 0.2 m 
omega = 1; % in rad/s
number_time_steps = 101;

theta_array = zeros(number_time_steps, 1);
x_array = zeros(number_time_steps, 1);
theta_derivativ_array = zeros(number_time_steps, 1);
x_derivativ_array = zeros(number_time_steps, 1);

% set a reasonable starting point
u0 = [pi ; 0];

%%
i = 1;
u = u0;

for t = linspace(0, 1, number_time_steps)
    
    phi = pi / 6 + omega * t;
    
    % create function handles
    F = @(u) constraint(u, a, b, phi);
    J = @(u) jacobian(u, b);
       
    eps = 1e-9;
    
    [u, iteration_counter] = NR_method(F, J, u, eps); 
    
    df_dt = [- a * omega * sin(omega * t)
            a * omega * cos(omega * t)];
    df_dq = [- b * sin(u(1)), -1
            - b * cos(u(1)), 0];

    u_derivative = - df_dt \ df_dq;
    
    theta_array(i) = u(1);
    x_array(i) = u(2);
    theta_derivativ_array(i) = u_derivative(1);
    x_derivativ_array(i) = u_derivative(2);
    
    i = i + 1;
end

%%
t = linspace(0, 1, number_time_steps);

homework_4_results = figure
tiledlayout(2,2)

nexttile
plot(t, theta_array);
title('Theta');

nexttile
plot(t, x_array);
title('X');

nexttile
plot(t, theta_derivativ_array);
title('Theta derivative');

nexttile
plot(t, x_derivativ_array);
title('X derivative');

saveas(homework_4_results, 'homework_4_results.png');
% fprintf('\nMechanism valid position is for d = %.3g m and theta = %g deg\n\n', ...
%         u(t, 2), rad2deg(u(1)));

%%
function P = constraint(u, a, b, phi)
theta = u(1);
d = u(2);

P = [a * cos(phi) + b * cos(theta) - d
    a * sin(phi) - b * sin(theta)];
end

%%
function P = jacobian(u, b)
theta = u(1);
P = [- b * sin(theta), -1
    - b * cos(theta), 0];
end
