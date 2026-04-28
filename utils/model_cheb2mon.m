function mon_model = model_cheb2mon(model)

% Transform a model in scaled chebyshev polynomials to a model in the
% monomial basis. This helps for comparison with ground truths.

% Too many variables?
assert(model.dim<=5, "Conversion too expensive!")

% % Proceed for each coordinate
% for i = 1:model.dim
%     df = max(sum(model.basis{i},2));
%     d = max(df, 10);
%     [X{1:model.dim}] = ndgrid(chebpts(d));
%     X = cellfun(@(X) reshape(X,[],1), X, 'UniformOutput', false);
%     pts_cheb = cell2mat(X);
%     pts_cheb = select_fekete(model.basis{i}, pts_cheb);
%     pts_mons = ( pts_cheb - model.rescale.mu ) / model.rescale.Lambda';
%     nPts = size(pts_cheb,1);
%     fval = cheb_evaluate(model.basis{i}, pts_cheb) * model.coeff{i};
%     Q = zeros(nPts, nPts);
%     for j = 1:nPts
%         Q(:,j) = prod(pts_mons.^model.basis{i}(j,:), 2);
%     end
%     mon_model.coeff{i} = ( ( Q \ fval ) .* model.rescale.a(i) );
%     mon_model.basis{i} = model.basis{i};
% end

% Create union of bases
mon_model.basis = zeros(0,model.dim);
for i = 1:model.dim
    mon_model.basis = union(mon_model.basis, model.basis{i}, 'rows');
end

% Create points for evaluation
df = max(sum(mon_model.basis,2));
d = max(df, 10);
[X{1:model.dim}] = ndgrid(chebpts(d));
X = cellfun(@(X) reshape(X,[],1), X, 'UniformOutput', false);
pts_cheb = cell2mat(X);
pts_cheb = select_fekete(model.basis{i}, pts_cheb);
pts_mons = ( pts_cheb - model.rescale.mu ) / model.rescale.Lambda';
nPts = size(pts_mons,1);

% Monomial evaluations
Q = zeros(nPts, nPts);
for j = 1:nPts
    Q(:,j) = prod(pts_mons.^mon_model.basis(j,:), 2);
end

% Evalue model and convert to monomial basis
fval = model_evaluate(pts_mons, model);
mon_model.coeff = Q \ fval;

end