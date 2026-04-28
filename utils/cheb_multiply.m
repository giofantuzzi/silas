function gamma = cheb_multiply(alpha, beta)
% Function to multiply two multivariate chebyshev "monomials"

% Check inputs
n = numel(alpha);
m = numel(beta);
assert(m==n, 'Inputs should be vectors of equal size');
% Multiply
if n > 1
    g1 = cheb_multiply(alpha(1:n-1), beta(1:n-1));
    g2 = cheb_multiply(alpha(n), beta(n));
    gamma   = [ g1, repmat(g2(1), size(g1,1), 1); g1, repmat(g2(2), size(g1,1), 1)];
else
    gamma   = [alpha+beta; abs(alpha-beta)];
end
end