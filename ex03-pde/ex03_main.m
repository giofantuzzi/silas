% Data-driven discovery of a stable model for ROM model of a reaction diffusion PDE 
% Authors: Albert Alcalde, Giovanni Fantuzzi

% Clean up
clear, clc
yalmip clear
close all
addpath('../utils/')

% SILAS parameters
param = make_parameters();
param.dv = 2;               % degree of Lyapunov dictionary *MUST BE EVEN and >=2*
param.df = 5;              % degree of model *>=1* (Better if odd but the code handles all cases)
param.beta = 1e6;

% Other parameters
% use N = 45000 for n_pod_modes = 3,4 and N = 270000 for n_pod_modes = 5
n_pod_modes = 3;
N = 45000; 
save_plots = true;

% Seed random number generator
rng(0)

% Learn model
data.train = ['./data', filesep, sprintf('PODcoeffs_d=%d_N=%4d_train.hdf5', n_pod_modes, N)];
data.test  = ['./data', filesep, sprintf('PODcoeffs_d=%d_N=%4d_test.hdf5', n_pod_modes, N)];
[model, lf_initial, lf_final, err] = model_learn(data.train, data.test, param);
save_name = sprintf('./models/reacDiff_d=%d_dv%d_df%d.hdf5', n_pod_modes, param.dv, param.df);
save_results(save_name, model, lf_initial, lf_final, err, data)

% Plot
plot_pde_solution(model, param, data, save_plots)
plot_attractors(model, lf_final, param, data, save_plots)

% Clean up
rmpath('../utils/')