function plot_attractors(model, lf_final, param, data, save_plot)

% Plot POD-coefficient attractors for the ROM and the original PDE data

% Load the data
xi_train = h5read(data.train, '/x')';
xi_test = h5read(data.test, '/x')';
t = h5read(data.test, '/t');

% Some sizes
Ntrain = size(xi_train, 1); % 5000
Ntest = Ntrain + 1;
t_test = t(Ntest:end);

% Simulate model
xi0 = xi_test(1,:);
f_learned = @(tt, xi) model_evaluate(xi.', model).';
opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
[~, xi_learned] = ode45(f_learned, t_test, xi0, opts);

% Plot parameters
axlim = [-100 20 0 300];
axPos = [0.33, 0.24, 0.55, 0.65];
name = {'true', 'learned', 'trap'};

% Loop over plots
for i = 1:3

    % Make figure
    fig = figure();
    fig.Units = 'pixels';
    fig.Position([3 4]) = [128 128];
    axis square; box on
    ylabel('$x_2$','Interpreter','latex','FontSize',8)
    ax = gca;
    ax.FontSize = 8;
    ax.TickLabelInterpreter = 'latex';
    axis(axlim)
    ax.Position = axPos;

    % Fix figure on paper
    fig.Units = 'centimeters';
    fig.PaperUnits = 'centimeters';
    fig.PaperSize = fig.Position([3 4]);
    fig.PaperPosition = [0 0 fig.Position([3 4])];

    % Set limits
    if i < 3
        absSetLim = 1 * [-0.4 0.3 -0.02 0.07]; 
    else
        % Plot the absorbing set, too
        patch_clr = [0.4660    0.6740    0.1880];
        XY = lyap_plotProjDeg2(lf_final, model, [1 2], 400);
        p = patch(XY(:,1), XY(:,2), patch_clr);
        p.FaceAlpha = 0.25;
        p.LineStyle = 'none';
        hold on
        xbox = [absSetLim(1) absSetLim(2) absSetLim(2) absSetLim(1) absSetLim(1)];
        ybox = [absSetLim(3) absSetLim(3) absSetLim(4) absSetLim(4) absSetLim(3)];
        plot(xbox, ybox, 'r-')
        absSetLim = 10 * [-0.4 0.4 -0.06 0.06]; % for n = 5
    end

    % Plot attractors
    hold on
    ax.Children = flip(ax.Children);
    axis(absSetLim)
    ax.Position = axPos;
    if i == 1 % PDE projection
        s = scatter(xi_test(:,1),xi_test(:,2), 0.75, 'o', 'filled');
        s.MarkerEdgeAlpha = 0.1;
        s.MarkerFaceAlpha = 0.1;
    else % ROM
        s = scatter(xi_learned(:,1),xi_learned(:,2), 0.75, 'o', 'filled');
        s.MarkerEdgeAlpha = 0.1;
        s.MarkerFaceAlpha = 0.1;
        s.MarkerFaceColor = [0 0.4470 0.7410];
        xlabel('$x_1$','Interpreter','latex','FontSize',8)
    end

    if save_plot
        fname = sprintf('./img/rd-attr-%s-d=%d_dv%d_df%d_beta%d.png', name{i}, model.dim, param.dv, param.df, param.beta);
        print(fname,'-dpng','-r600')
        close(fig)
    end
end
