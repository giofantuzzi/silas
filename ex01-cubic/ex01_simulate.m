% Simulate the dynamics of a two-dimensional cubic ODE with 5 equilibria,
% an attracting limit cycle, and a compact absorbing set described by a
% quartic Lyapunov function
%
% Authors: Albert Alcalde & Giovanni Fantuzzi

% Clean up
clear, clc
close all

% Parameters
rng(0)
savedata = true;
nPoints = 40;      % points per unit time, so dt = 1/(nPoints-1)
T = 10;            % final simulation time (change me)

% Dynamics
f1 = @(x) x(1)-x(2) - x(1).*x(2).^2 + x(2).^3; 
f2 = @(x) x(1) + x(2) - x(2).*x(1).^2 - x(1).^3;
dynamics = @(x) [f1(x); f2(x)];

% Simulate for a long time, until we essetially reach the limit cycle
% Thenm collect data as requested
tot_points = round(T*nPoints);
opts = odeset('RelTol',1e-10,'AbsTol',1e-12);
[~, xsim] = ode45(@(t,x)dynamics(x), [0, 1e2], [0.5, 1], opts);
[tsim, xsim] = ode45(@(t,x)dynamics(x), linspace(0,T,tot_points), xsim(end,:), opts);

% ODE vector field, this time working on columns of matrices
% We use this to save the exact vector field at the data points
f1 = @(x) x(:,1)-x(:,2) - x(:,1).*x(:,2).^2 + x(:,2).^3; 
f2 = @(x) x(:,1) + x(:,2) - x(:,2).*x(:,1).^2 - x(:,1).^3;
f = [f1(xsim), f2(xsim)];

% Save data
% NOTE: WE save the data in row-major ordering consistent with python,
% which is used to generate data for other examples.
dt = T / (tot_points-1);
fname = sprintf('data/data_n%i.hdf5', tot_points);
if isfile(fname); delete(fname); end % overwrite previous runs
h5create(fname,'/dt',[1 1]);
h5create(fname,'/x',[2, tot_points]);
h5create(fname,'/f',[2, tot_points]);
h5write(fname, '/dt', dt)
h5write(fname, '/x', xsim'); % consistent with python row-major ordering
h5write(fname, '/f', f');    % consistent with python row-major ordering

% Plot the dataset (clean)
plot(xsim(:,1), xsim(:,2), '.')