function plot_pde_solution(model, param, data, save_plot)

% Plot PDE solution recovered from the ROM and compare to the simulation
% data

% Load the data
xi_train = h5read(data.train, '/x')';
xi_test = h5read(data.test, '/x')';
u_test = h5read(data.test, '/u_test')';
Y_mean = h5read(data.test, '/Y_mean')';
Phi = h5read(data.test, '/Phi')';
t = h5read(data.test, '/t');
x_grid = h5read(data.test, '/x_grid');

% Some sizes
Nx = length(x_grid);
Ntrain = size(xi_train, 1); % 5000
Ntest = Ntrain + 1;
t_test = t(Ntest:end);

% Simulate model
xi0 = xi_test(1,:);
f_learned = @(tt, xi) model_evaluate(xi.', model).';
opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
[~, xi_learned] = ode45(f_learned, t_test, xi0, opts);
Y_learned = Y_mean + Phi * xi_learned.';   % (2*Nx) x Nt
u_learned = Y_learned(1:Nx, :);

% Compute errors and axis limits
u_err = abs(u_test - u_learned);
u_min = min([u_test, u_learned], [], 'all');
u_max = max([u_test, u_learned], [], 'all');

% Make figure
fig = figure('Position',[100 100 1200 240]);
tiledlayout(1,3,'TileSpacing','compact','Padding','compact')
ax_fontsize = 16;
cb_fontsize = 16;

% Exact solution
ax = nexttile;
pcolor(ax, t_test, x_grid, u_test)
shading(ax, 'interp')
axis(ax, 'tight')
box(ax, 'on')
ax.FontSize = ax_fontsize;
ax.LineWidth = 1;
ax.Layer = 'top';
ax.TickLabelInterpreter = 'latex';
xlabel(ax, '$t$', 'Interpreter', 'latex', 'FontSize', ax_fontsize)
ylabel(ax, '$\xi$', 'Interpreter', 'latex', 'FontSize', ax_fontsize)
cb2 = colorbar(ax);
cb2.FontSize = cb_fontsize;
cb2.Label.Interpreter = 'latex';
cb2.TickLabelInterpreter = 'latex';
colormap(ax, 'viridis')
clim(ax, [u_min, u_max])

% Predicted solution
ax = nexttile;
pcolor(ax, t_test, x_grid, u_learned)
shading(ax, 'interp')
axis(ax, 'tight')
box(ax, 'on')
ax.FontSize = ax_fontsize;
ax.LineWidth = 1;
ax.Layer = 'top';
ax.TickLabelInterpreter = 'latex';
xlabel(ax, '$t$', 'Interpreter', 'latex', 'FontSize', ax_fontsize)
ylabel(ax, '$\xi$', 'Interpreter', 'latex', 'FontSize', ax_fontsize)
cb2 = colorbar(ax);
cb2.FontSize = cb_fontsize;
cb2.Label.Interpreter = 'latex';
cb2.TickLabelInterpreter = 'latex';
colormap(ax, 'viridis')
clim(ax, [u_min, u_max])

% Pointwise absolute error
ax = nexttile;
pcolor(ax, t_test, x_grid, u_err)
shading(ax, 'interp')
axis(ax, 'tight')
box(ax, 'on')
ax.FontSize = ax_fontsize;
ax.LineWidth = 1;
ax.Layer = 'top';
ax.TickLabelInterpreter = 'latex';
xlabel(ax, '$t$', 'Interpreter', 'latex', 'FontSize', ax_fontsize)
ylabel(ax, '$\xi$', 'Interpreter', 'latex', 'FontSize', ax_fontsize)
colormap(ax, 'turbo')
cb2 = colorbar(ax);
clim(ax, [0, 1])
cb2.FontSize = cb_fontsize;
cb2.Label.Interpreter = 'latex';
cb2.TickLabelInterpreter = 'latex';

% Save plot?
if save_plot
    print(sprintf('./img/rd-pde-state-d=%d_dv%d_df%d_beta%d.png', model.dim, param.dv, param.df, param.beta),'-dpng','-r600')
    close(fig)
end