% Simulate discovered ODEs and plot projections of the exact attractor,
% discovered attractor, and discovered absorbing set on the (x1,x3) plane.
% Note: the simulation and absorbing set computations for the discovered
% dynamics are currently rather inefficient, but they work.
%
% Authors: Albert Alcalde, Giovanni Fantuzzi

% clean up
clear, clc
close all

% Path
addpath("../utils/")
patch_clr = [0.4660    0.6740    0.1880];
axPos = [0.300938010215744, ...
         0.207760418752829, ...
         0.604061989784256, ...
         0.717239581247171];

% Parameters
model_dir = 'models-refined';
data_dir  = 'data-train-refined';
odesolver = @ode89;
opts = odeset('RelTol', 1e-10, 'AbsTol',1e-10);
T = 1;
plot_AS = true;
asPoints = 100;
figsz = 128;
transparency = 0.1;

% Axis limits for figures without absorbing set
ax_lim = [
     0 1 -1 1;      % AtmosphericRegime
    -4 4 -0.2 0.2;  % Bouali 2
    -5 5 -2 1;      % GuckenheimerHolmes
    -100 20 0 300;  % HyperRossler
    -3 3 -3 3;      % SprottA
    -5 5 -5 5;      % Thomas
    -3 3 0 60;      % YuWang
    ];

% Axis limits for figures with absorbing sets
absSetLim = [
    -1 2 -2 2.5;        % AtmosphericRegime
    -20 20 -0.5 0.5;    % Bouali 2
    -25 15 -5 2;        % GuckenheimerHolmes
    -500 400 -600 1000; % HyperRossler
    -15 15 -15 15;      % SprottA
    -10 10 -10 10;      % Thomas
    -5 5 -100 150;      % YuWang
    ];

% System names
sys_names = {
    'AtmosphericRegime'
    'Bouali2'
    'GuckenheimerHolmes'
    'HyperRossler'
    'SprottA'
    'Thomas'
    'YuWang'
    };

% System ODE functions for exact simulation
sys_fun = {
    @(t,x) exact_AtmosphericRegime(x)
    @(t,x) exact_Bouali2(x)
    @(t,x) exact_GuckenheimerHolmes(x)
    @(t,x) exact_HyperRossler(x)
    @(t,x) exact_SprottA(x)
    @(t,x) exact_Thomas(x)
    @(t,x) exact_YuWang(x)
    };

% Loop over available data files
data_files = dir([model_dir, filesep, '*.hdf5']);
nFiles = numel(data_files);
nNames = numel(sys_names);
for name_idx = 1:nNames
    for k = 1:nFiles

        current_name = ['_', sys_names{name_idx}, '.hdf5'];

        if contains(data_files(k).name, current_name, 'IgnoreCase', true)
            % Display header
            current_name = data_files(k).name;
            fprintf('MODEL: %s (%i/%i)\n', current_name, name_idx, nNames)

            % Get model training data for IC
            x_train = h5read([data_dir, filesep, current_name], '/x')';
            x0 = x_train(1,:);

            % Load the model
            model = model_load([model_dir, filesep, current_name]);
            lf = lyap_load([model_dir, filesep, current_name]);
            
            % Simulate the exact model
            fprintf('\t Simulating exact model...\n')
            [texat, xexact] = odesolver(sys_fun{name_idx}, [0, T], x0, opts);

            % Simulate learned model with ode45
            fprintf('\t Simulating discovered model (patience!)...\n')
            make_row = @(x) reshape(x, [], model.dim);
            make_col = @(x) reshape(x, model.dim, []);
            fun = @(t, x) make_col( model_evaluate(make_row(x), model) );
            [tsim, xsim] = odesolver(fun, [0, T], x0, opts);

            % Plot exact attractor
            close all
            fig = figure();
            fig.Units = 'pixels';
            fig.Position([3 4]) = [figsz figsz];
            axis(ax_lim(name_idx,:))
            s = scatter(xexact(:,1),xexact(:,3), 0.25, 'o', 'filled');
            s.MarkerEdgeAlpha = transparency;
            s.MarkerFaceAlpha = transparency;
            axis square; box on
            xlabel('$x_1$','Interpreter','latex','FontSize',8)
            ylabel('$x_3$','Interpreter','latex','FontSize',8)
            ax = gca;
            ax.FontSize = 8;
            ax.TickLabelInterpreter = 'latex';
            axis(ax_lim(name_idx,:))
            ax.Position = axPos;
            fig.Units = 'pixels';
            fig.Position([3 4]) = [figsz figsz];
            fig.Units = 'centimeters';
            fig.PaperUnits = 'centimeters';
            fig.PaperSize = fig.Position([3 4]);
            fig.PaperPosition = [0 0 fig.Position([3 4])];
            drawnow
            fname = sprintf('img/%s_exact',sys_names{name_idx});
            print(fname,'-dpng','-r600')

            % Plot discovered attractor
            close all
            fig = figure();
            fig.Units = 'pixels';
            fig.Position([3 4]) = [figsz figsz];
            axis(ax_lim(name_idx,:))
            s = scatter(xsim(:,1),xsim(:,3), 0.25, 'o', 'filled');
            s.MarkerEdgeAlpha = transparency;
            s.MarkerFaceAlpha = transparency;
            axis square; box on
            xlabel('$x_1$','Interpreter','latex','FontSize',8)
            ylabel('$x_3$','Interpreter','latex','FontSize',8)
            ax = gca;
            ax.FontSize = 8;
            ax.TickLabelInterpreter = 'latex';
            axis(ax_lim(name_idx,:))
            ax.Position = axPos;
            fig.Units = 'pixels';
            fig.Position([3 4]) = [figsz figsz];
            fig.Units = 'centimeters';
            fig.PaperUnits = 'centimeters';
            fig.PaperSize = fig.Position([3 4]);
            fig.PaperPosition = [0 0 fig.Position([3 4])];
            drawnow
            fname = sprintf('img/%s_model',sys_names{name_idx});
            print(fname,'-dpng','-r600')

            % Plot absorbing set from Lyapunov function
            % We evaluate it on the hypercube and scale back
            if plot_AS && model.dim <= 4 && lf.c > 0
                fprintf('\t Computing absorbing set projection (patience!)...\n')
                XY = plot_lyap(current_name, model_dir, asPoints);
                p = patch(XY(:,1), XY(:,2), patch_clr);
                p.FaceAlpha = 0.25;
                p.LineStyle = 'none';
                ax.Children = flip(ax.Children);
                axis(absSetLim(name_idx,:))
                ax.Position = axPos;
                fname = sprintf('img/%s_absSet',sys_names{name_idx});
                drawnow
                print(fname,'-dpng','-r600')
            end

        end

    end
end


% clean up path
rmpath("../utils/")