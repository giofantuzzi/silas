function model = model_load(model_name)

% Load model from hdf5 file

% % Initialize
% if nargin < 2; group_name = '/model'; end
% model = struct();
% 
% % Get model info
% info = h5info(model_name, group_name);
% 
% % Load dataset in loop
% nDataset = numel(info.Datasets);
% for i = 1:nDataset
%     data_name = [group_name, '/', info.Datasets(i).Name];
%     data = h5read(model_name,data_name);
%     model = setfield(model, info.Datasets(i).Name, data);
% end
% 
% % Load groups in loop
% nGroups = numel(info.Groups);
% for i = 1:nGroups
%     group_data = model_load(model_name, info.Groups(i).Name);
%     [~, data_name] = fileparts(info.Groups(i).Name);
%     model = setfield(model, data_name, group_data);
% end

% Basic model info
model = model_initialize();
model.name = h5read(model_name, '/model/name');
model.dim = h5read(model_name, '/model/dim');
model.degree = h5read(model_name, '/model/degree');

% Model basis and coefficients
for i = 1:model.dim
    model.basis{i} = h5read(model_name, sprintf('/model/basis/f%i', i));
    model.coeff{i} = h5read(model_name, sprintf('/model/coeff/f%i', i));
end

% Rescaling stuff
model.rescale.Lambda = h5read(model_name, '/model/rescale/Lambda');
model.rescale.mu = h5read(model_name, '/model/rescale/mu');

% Load model error
model.error.max = h5read(model_name, '/model/error/max');
model.error.avg = h5read(model_name, '/model/error/avg');
model.error.std = h5read(model_name, '/model/error/std');
model.error.raw = h5read(model_name, '/model/error/raw');