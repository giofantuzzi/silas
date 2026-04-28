function [model, lf_initial, lf_final, err] = model_learn(file_train, file_test, param)

% Learn a stable model for a one of the dysts dynamical systems
% INPUTS:
% * file_train = name of file with training data (HDF5 with /x and /dt fields)
% * file_test = name of file with test data to evaluate model errors (HDF5 with /x and /f fields)
% * param = structure with various parameters
%
% METHOD:
% 1) We approximate the Lie derivative using EDMD on noisy data
% 2) We find an approximate Lyapunov function and trapping region
% 3) We learn a stable model with the given trapping region

% Clean up
yalmip clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load the training / test data and perform initial checks
% Load the training data from file
dt = h5read(file_train,'/dt')';
x_train = h5read(file_train,'/x')';

% Small enough dynamics?
model_dimension = size(x_train, 2);
if model_dimension > 6
    % too big, give up
    fprintf('Dimension = %i: skipping...\n', model_dimension)
    model = model_initialize(file_train);
    lf_initial = lyap_initialize();
    lf_final = lyap_initialize();
    err = [];
    return
% elseif model_dimension > 4
%     fprintf('Dimension = %i: limit to df <= 7 and dv <= 2\n', model_dimension)
%     param.dv = param.dv( param.dv <= 2 );
end

% Load test data
x_test = h5read(file_test,'/x')';
f_test = h5read(file_test,'/f')';
f_test_norm = vecnorm(f_test, 2, 2);

% Initialize model
[~,model_name,~] = fileparts(file_train);
model = model_initialize(file_train);
lf_initial = lyap_initialize();
lf_final = lyap_initialize();
err.dv = param.dv;
err.df = param.df;
err.max = Inf(numel(err.df), numel(err.dv));
err.avg = Inf(numel(err.df), numel(err.dv));
err.std = Inf(numel(err.df), numel(err.dv));

% Add noise to the data
% Mean zero and standard deviation equal to param.noise * std(data)
if param.noise > 0
    std_noise = sqrt(pi/2) .* param.noise .* x_train;
    x_train = x_train + std_noise .* randn(size(x_train));
end

% Rescale the data to make sure the datapoints are inside the unit box
if isempty(param.Lambda) || isempty(param.mu)
    % use default rescaling
    ub = max(x_train); ub = ub + 0.25*abs(ub);
    lb = min(x_train); lb = lb - 0.25*abs(lb);
    param.Lambda = diag( 2./(ub - lb) );
    param.mu = -0.5*(ub + lb) * param.Lambda';
end
data.x = x_train * param.Lambda' + param.mu;
model.rescale.Lambda = param.Lambda;    % save inside model for evaluation
model.rescale.mu = param.mu;            % save inside model for evaluation

% Approximate the vector field (FD scheme of order 8)
wts = (1/dt)*[1/280, -4/105, 1/5, -4/5, 0, 4/5, -1/5, 4/105, -1/280]';
data.y = zeros(size(data.x,1)-8, size(data.x,2));
for j = 1:numel(wts)
    data.y = data.y + wts(j) * data.x(j:end-9+j,:);
end
data.x = data.x(5:end-4,:);
[data.n, model.dim] = size(data.x);
param.nx = model.dim;
model.num_data = data.n;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extract hyperparameters and loop
% We extract the hyperparameters over which we may want to loop
coerc_vals = param.eps_coerc;
penalties  = param.eps;
numOpt = numel(param.dv) * numel(param.df) * numel(coerc_vals) * size(penalties, 2);
model_iter = 0;
% fprintf("***********************************************************\n")
% fprintf('MODEL: %s\n', model_name)
% fprintf("***********************************************************\n")
fprintf('|    Iter   | df | dv |   Error  | Updated |\n')
fprintf("|-----------|----|----|----------|---------|\n")
progress_str_1 = '| %4i/%4i | %2i | %2i ';
progress_str_2 = '| %6.2e |';
for i_dv = 1:numel(err.dv)

    % Select LF degree
    param.dv = err.dv(i_dv);

    for i_df = 1:numel(err.df)

        % Select model degree
        param.df = err.df(i_df);

        for i_coerc = 1:numel(coerc_vals)

            % Select the coercivity constant
            param.eps_coerc = coerc_vals(i_coerc);

            for i_pen = 1:size(penalties, 2)

                % Select the ROA penalty parameter
                param.eps = penalties(:, i_pen);
        
                % Update iteration count
                model_iter = model_iter + 1;
                fprintf(progress_str_1, model_iter, numOpt, param.df, param.dv);

                % Build initial Lyapunov function
                % Try the constrained problem
                % Revert to penalized problem if this fails
                lf_initial = build_initial_lyap_constrained(data, param);
                if ~lf_initial.feasible || lf_initial.num_prob
                    lf_initial = build_initial_lyap_penalized(data, param);
                end

                % Did we succeed or did we have numerical problems?
                if lf_initial.num_prob
                    % Stop the computation, we did rubbish
                    fprintf('|   NaN    |  FALSE  |\n')
                    break
                end

                % Build the model and re-optimize the Lyapunov function
                new_lf = lf_initial;
                new_model = model;
                for i = 1:param.maxIter
                    % Try to iterate over model and LF
                    candidate_model = build_model_lyap(data, new_lf, param);
                    candidate_lf = build_lyap_model(candidate_model, param);
                    if (candidate_model.num_prob) || (candidate_lf.num_prob)
                        % We encountered numerical problems, stop
                        break
                    else
                        % Update and move to next iteration
                        new_model = candidate_model;
                        new_lf = candidate_lf;
                        new_model.iterations = i;
                    end
                end

                % Did we succeed at all? We check again to catch cases when
                % no iteration runs.
                if isempty(new_model.coeff)
                    % Stop the computation, we did rubbish
                    fprintf('|   NaN    |  FALSE  |\n')
                    break
                end

                % Evaluate the model's accuracy
                % Compare learned and exact vector fields on test data
                f_model = model_evaluate(x_test, new_model);
                err_norm = vecnorm(f_test-f_model, 2, 2);
                new_model.error.raw = err_norm ./ f_test_norm;
                new_model.error.max = max( new_model.error.raw );
                new_model.error.avg = mean( new_model.error.raw );
                new_model.error.std = std( new_model.error.raw );
                fprintf(progress_str_2, new_model.error.avg);

                % Record model error for diagnostics
                update_err = new_model.error.avg < err.avg(i_df, i_dv);
                update_err = update_err && ~(new_model.num_prob);
                update_err = update_err && ~(new_lf.num_prob);
                if update_err
                    err.max(i_df, i_dv) = new_model.error.max;
                    err.avg(i_df, i_dv) = new_model.error.avg;
                    err.std(i_df, i_dv) = new_model.error.std;
                end

                % Do we keep this model?
                update = new_model.error.avg < model.error.avg;
                update = update && ~(new_model.num_prob);
                update = update && ~(new_lf.num_prob);
                if update
                    % Update model, add name and recompute degree
                    model = new_model;
                    model.name = model_name;
                    model.num_data = data.n;
                    lf_final = new_lf;
                    fprintf('  TRUE   |\n')
                else
                    fprintf('  FALSE  |\n')
                end

            end % loop over ROA penalty
        end % loop over coercivity constant
    end % loop over model degrees
end % loop over LF degrees

% We are done
fprintf("|------------------------------------------|\n")


%%% THE END
end