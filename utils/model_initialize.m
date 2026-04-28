function model = model_initialize(model_name)

% Initialize an empty model as a structure with all of the fields we will
% need throughout the computation.
if nargin > 0
    model.name = model_name;
else
    model.name = [];
end
model.dim   = [];
model.degree = [];
model.basis = [];
model.coeff = [];
model.rescale.Lambda = [];
model.rescale.mu = [];
model.error.max = +Inf;
model.error.avg = +Inf;
model.error.std = +Inf;
model.error.raw = [];
model.num_prob = [];
model.num_data = [];
model.opt = [];
model.iterations = 0;