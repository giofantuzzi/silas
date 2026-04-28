% Data-driven discovery of a stable model for a 2D cubic ODE without
% absorbing ellipsoid, but with an absorbing set described by a quartic
% polynomial Lyapunov function
%
% Authors: Albert Alcalde & Giovanni Fantuzzi

% Clean up
clear, clc
yalmip clear
close all
addpath('../utils/')

% Parameters for ODE discovery
df_vals = 2:9;
param = make_parameters();
param.dv = 4;               % degree of Lyapunov dictionary *MUST BE EVEN and >=2*
param.maxIter = 5;          % how many alternating minimizations
param.beta = 1e2;           % force b=0 if c=0
param.eps(1) = 1e-1;        % penalty for c==1
param.noise = 0.01;         % noise parameter: E(|error|) = this value * exact measurement
param.Lambda = diag( 2./[5 5] ); % rescaling matrix
param.mu = zeros(1,2);       % rescaling vector

% Plot parameters
plotop.ax_lim = 3 * [-1 1 -1 1];
plotop.box_sz = 3;
plotop.show_box = false;
plotop.export = false;
plotop.shade_clr = [0.4660    0.6740    0.1880]; % green color

% ODE vector field, operating on rows of a matrix
f_true = @(X) [X(:,1) - X(:,2) - X(:,1) .* X(:,2).^2 + X(:,2).^3, ...
    X(:,2) + X(:,1) - X(:,2) .* X(:,1).^2 - X(:,1).^3];

% Loop over model degrees
for df = df_vals

    % clean up and reset random number generator
    yalmip clear
    rng(0)

    % Learn model
    param.df = df;
    data.train = 'data/data_n100.hdf5';
    data.test  = 'data/data_n100.hdf5';
    [model, lf_initial, lf_final, err] = model_learn(data.train, data.test, param);

    % Evaluate model accuracy
    x = linspace(plotop.ax_lim(1), plotop.ax_lim(2), 501);
    y = linspace(plotop.ax_lim(3), plotop.ax_lim(4), 501);
    [x, y] = meshgrid(x, y);
    XY = [x(:), y(:)];
    F_learn = model_evaluate(XY, model);
    ptwise_2 = vecnorm(f_true(XY) - F_learn, 2, 2);
    ptwise_2 = reshape(ptwise_2, size(x));

    % Make plots
    plot_colorbar(plotop)
    plot_error(x, y, ptwise_2, plotop, param)
    plot_statespace(x, y, lf_final, model, plotop, param, data)

end

% Clean up
rmpath('../utils/')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fig, ax] = make_figure(px_size)
if nargin < 1; px_size = [128 128]; end
fig = figure();
fig.Units = 'pixels';
fig.Position([3 4]) = px_size;
fig.Units = 'centimeters';
fig.PaperUnits = 'centimeters';
fig.PaperSize = fig.Position([3 4]);
fig.PaperPosition = [0 0 fig.Position([3 4])];
ax = gca;
ax.FontSize = 8;
ax.TickLabelInterpreter = 'latex';
ax.Position = [0.200789251342866   0.242341600002993   0.704210748657134   0.682658399997007];
axis square; box on
hold on
xlabel('$x_1$', 'Interpreter','latex', 'FontSize', 8)
ylabel('$x_2$', 'Interpreter','latex', 'FontSize', 8)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_error(x, y, z, opt, param)
[fig, ax] = make_figure();
figure(fig)
pcolor(x, y, z)
shading interp
axis equal
hold on
colormap(viridis)
clim([1e-3, 1e2])
set(gca, 'ColorScale','log')
axis(ax, opt.ax_lim)
ax.XTick = opt.ax_lim(1) : opt.ax_lim(2);
ax.YTick = opt.ax_lim(3) : opt.ax_lim(4);
if opt.export
    print(fig, sprintf('./img/ex01-error-dv%d_df%d.png',param.dv, param.df),'-dpng','-r600')
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_colorbar(opt)
fig = figure();
fig.Units = 'pixels';
fig.Position([3 4]) = [216 128];
colormap(viridis)
ax = gca;
clim([1e-3, 1e2])
set(gca, 'ColorScale','log')
cb = colorbar(ax);
cb.FontSize = 8;
cb.TickLabelInterpreter = 'latex';
cb.Label.Interpreter = 'latex';
cb.Location = 'southoutside';
cb.Ticks = 10.^(-3:2);
cb.TickLabels = {'$10^{-3}$','$10^{-2}$','$10^{-1}$','$10^{0}$','$10^{1}$','$10^{2}$'};
ax.Visible = false;
if opt.export
    print(fig, './img/ex01-colorbar.png','-dpng','-r600')
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_statespace(x, y, lf, model, opt, param, data)
[fig, ax] = make_figure();
% Evaluate lyapunov function on the grid and find sublevel set
XY = [x(:), y(:)];
XY2 = XY * model.rescale.Lambda' + model.rescale.mu;
V = cheb_evaluate(lf.basis, XY2) * lf.coeff(:);
V = reshape(V, size(x));
mask = (V < lf.b / lf.c);
% Plot the absorbing set
RGB = zeros([size(mask), 3]);
RGB(:,:,1) = opt.shade_clr(1);
RGB(:,:,2) = opt.shade_clr(2);
RGB(:,:,3) = opt.shade_clr(3);
h = image(ax, 'XData', opt.ax_lim(1:2), 'YData', opt.ax_lim(3:4), 'CData', RGB);
h.AlphaData = 0.25 * mask;
% streamlines of learned field
Fs  = model_evaluate(XY, model);
F1 = reshape(Fs(:,1), size(x));
F2 = reshape(Fs(:,2), size(x));
h = streamslice(x, y, F1, F2, 1.5);
set(h,'Color',[1 1 1]*0.65)
set(h,'LineWidth',0.25)
% Plot the (noisy) training data
rng(0)
x_train = h5read(data.train, '/x');
if param.noise > 0
    std_noise = sqrt(pi/2) .* param.noise .* x_train;
    x_train = x_train + std_noise .* randn(size(x_train));
end
scatter(x_train(1,:), x_train(2,:), 0.25, 'b', 'o', 'filled');
% plot red box for visual reference
if opt.show_box
    plot(opt.box_sz*[1 -1 -1 1 1], opt.box_sz*[1 1 -1 -1 1], 'r-')
end
% export
if opt.export
    plot_name = './img/ex01-state-space-dv%d_df%d.png';
    print(fig,sprintf(plot_name, param.dv, param.df),'-dpng','-r600')
end
end