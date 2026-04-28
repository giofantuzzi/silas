% Learn stable models for systems from the dysts database
% Data is available in the data-* subfolders
%
% Authors: Albert Alcalde and Giovanni Fantuzzi

% clean up
clear, clc
close all

% Set path
addpath('./../utils')

% Data and system parameters
dir_train = './data-train-refined';
dir_test  = './data-test-refined';
dir_save  = './models-refined';
sys_indx = 1:3;

% Model learning parameters
param = make_parameters();
param.df = 2:7;
param.dv = [2 4];
param.alpha = 1e-1;
param.eps(3) = 1e-2;

% Get available data files
data_files = dir([dir_train, filesep, '*.hdf5']);
nFiles = numel(data_files);

% Loop over systems
nSystems = min(numel(sys_indx), nFiles);
for k = 1:nSystems
    current_name = data_files(sys_indx(k)).name;
    save_name = [dir_save, filesep, current_name];
    if 1%~exist(save_name, 'file')
        fprintf("******************************************************************\n")
        fprintf('MODEL: %s (%i/%i)\n', current_name, k, nSystems)
        fprintf("******************************************************************\n")
        data.train = [dir_train, filesep, current_name];
        data.test  = [dir_test , filesep, current_name];
        [model, lf_initial, lf_final, err] = model_learn(data.train, data.test, param);
        if ~isempty(model.basis)
            save_results(save_name, model, lf_initial, lf_final, err, data)
        end
    end
end
 
% Clean path
rmpath('./../utils')