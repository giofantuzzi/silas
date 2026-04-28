function [basis_d, coeff_d] = cheb_gradient(basis_v, coeff_v)

% Differentiate a chebyshev polynomial with coefficients coeff_v
% corresponding to the multivariate chebyshev basis basis_v

% Get polynomial info
[nTerms, nVars] = size(basis_v);

% Initialize the basis for the derivative
dv = max(sum(basis_v, 2));
basis_d = multi_indices(nVars, dv-1);
coeff_d = zeros(size(basis_d,1), nVars, 'like', coeff_v); % make sure the datatype matches!

for alpha = 1:nTerms
    for i = 1:nVars
        % Get coeffs of derivative of alpha-th basis term wrt x_i
        m = basis_v(alpha, i);
        if rem(m, 2)==0 % even degree
            beta = (1:2:m-1)';
            wts  = 2*ones(size(beta));
        else % odd degree
            beta = (0:2:m-1)';
            wts = 2*ones(size(beta));
            wts(1) = 1;
        end
        nBeta = numel(beta);
        pow = [repmat(basis_v(alpha, 1:i-1), nBeta, 1), ...
            beta, ...
            repmat(basis_v(alpha, i+1:end), nBeta, 1)];

        % Put the in the right place
        [~, idx] = ismember(pow, basis_d, 'rows');
        if ~isempty(idx)
            coeff_d(idx, i) = coeff_d(idx, i) + (coeff_v(alpha) * m) .* wts;
        end
    end
end