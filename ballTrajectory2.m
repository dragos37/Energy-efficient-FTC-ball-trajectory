clc, clearvars

xt = 4.16;
yt = 1;

v1 = 5;
v2 = 15;
theta0_min = 0.5;
theta0_max = 1.4;
theta0_step = 0.01;
tol = 0.01;

solution_found = false;



figure
hold on
xlabel('x (m)')
ylabel('y (m)')
yline(yt, 'r--', 'Target height')
xline(xt, 'g--', 'Target x')

min_energy = inf;
for theta0 = theta0_min : theta0_step : theta0_max
    v = bisect(theta0, xt, yt, v1, v2, tol);
    if ~isnan(v)
        energy = 0.5 * 0.07 * v^2; % E = mv^2 / 2
        if energy < min_energy
            min_energy = energy;
            optimal_theta0 = theta0;
            optimal_v = v;
            solution_found = true;
        end
    end
end

if ~solution_found
    disp('No solution found in the given angle range');
end

if solution_found
    fprintf('Optimal v0: %.4f m/s\n', optimal_v);
    fprintf('Optimal theta: %.4f rad (%.2f deg)\n', optimal_theta0, rad2deg(optimal_theta0));

    vx0 = optimal_v * cos(optimal_theta0);
    vy0 = optimal_v * sin(optimal_theta0);
    q0 = [0; 0; vx0; vy0];
    options_plot = odeset('Events', @(t,q) myEvent(t, q, yt));
    [~, q_full] = ode45(@ballODE, [0,10], q0, options_plot);
    plot(q_full(:,1), q_full(:,2));
    xlabel('x (m)')
    ylabel('y (m)')
    hold on
    plot(xt, yt, 'ro', 'MarkerSize', 10)  

end

function dqdt = ballODE(t, q)
% Unpack the state vector
x  = q(1);
y  = q(2);
vx = q(3);
vy = q(4);

m = 0.07;    
g = 9.81;   
b = 3.8*10^-5; %= 0.5*rho*Cd*A    

s = sqrt(vx^2 + vy^2); %speed

% Equations of motion
ax = -(b/m) * s * vx;
ay = -g - (b/m) * s * vy;

% Pack the derivatives into output vector
dqdt = [vx; vy; ax; ay];
end

function [value, isterminal, direction] = myEvent(t, q, yt)

value      = [q(2) - yt;  q(4)];  % two events: y=yt and vy=0
    isterminal = [1;           0];     % stop on y=yt, log apex
    direction  = [-1;          0];     % y descending, vy any
end

function error = shoot(v0, theta, xt, yt)
tspan = [0, 10];
vx0 = v0 * cos(theta);
vy0 = v0 * sin(theta);
x0 = 0;
y0 = 0;
q0 = [x0 y0 vx0 vy0];

options = odeset('Events', @(t,q) myEvent(t, q, yt));
[~, ~, te, qe, ie] = ode45(@ballODE, tspan, q0, options);
if isempty(qe)
    error = NaN;
    return
end
% ie tells you which event fired last
% if the stopping event (ie=1) fired without an apex (ie=2) first, reject
apex_times = te(ie == 2);
land_time  = te(ie == 1);
if isempty(apex_times) || isempty(land_time)
    error = NaN;
    return
end
if apex_times(end) > land_time
    error = NaN;  % apex after landing — ball never went above yt properly
    return
end
error = qe(end, 1) - xt;
end

function optimal_v = bisect(theta0, xt, yt, v1, v2, tol)
e_a = shoot(v1, theta0, xt, yt);
e_b = shoot(v2, theta0, xt, yt);

if isnan(e_a) || isnan(e_b) || sign(e_a) == sign(e_b)
    optimal_v = NaN;
    return
end

v_mid = (v1 + v2) / 2;

while v2 - v1 > tol
    v_mid = (v2 + v1)/2;

    e_mid = shoot(v_mid, theta0, xt, yt);
    if e_mid < 0
        v1 = v_mid;
    else
        v2 = v_mid;
    end



end  
optimal_v = v_mid;
end




