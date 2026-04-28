function f = model_evaluate(x, model)

% Evaluate a model at the points listed by the rows of x. The output has
% the same size as x.

% Check size
[nPts, nx] = size(x);
assert(nx==model.dim, 'Input x must have %i columns', model.dim)

% Rescale inputs to the unit box
x = x * model.rescale.Lambda' + model.rescale.mu;

% Initialize output
f = zeros(nPts, nx);
for i = 1:nx
    f(:,i) = cheb_evaluate(model.basis{i}, x) * model.coeff{i};
end

% Rescale vector field to original coordinates since we have a model for
% the scaled coordinates only
f = f / model.rescale.Lambda';