function [A, bases] = cheb_mommat(n, d, in_box)

% Create the "Lambda" operator for an n-dimensional hypercube using
% multivariate chebyshev basis for the upper and lower spaces. We use the
% standard description for the hypercube. We write A instead of \Lambda for
% the \Lambda operator because it is the closest standard letter :-)

% box constraints?
if nargin < 3
    in_box = false;
end

% Set up bases
p = multi_indices(n, d);
q = multi_indices(n, 2*d);
L = size(p, 1);
U = size(q, 1);
A = cell(n+1, 1);

% output the bases:
if nargout > 1
    bases.lower = p;
    bases.upper = q;
end

%-------------------------------------------------------------------------%
% The moment matrix
fun = @(i,j) cheb_multiply(p(i,:), p(j,:));
[I, J] = meshgrid(1:L, 1:L);
gamma  = arrayfun(fun, I, J, 'UniformOutput',false);
hash = rand(n,1);
gamma_hash = reshape(vertcat(gamma{:})*hash, 2^n, []);
qhash = q*hash;
iA = []; jA = []; vA = [];
for i = 1:U
    [~, ind] = find(abs(gamma_hash - qhash(i))<=1e-12);
    iA = [iA; ind];
    jA = [jA; repmat(i, numel(ind), 1)];
    vA = [vA; repmat(1/2^n, numel(ind), 1)];
end

% Are we done or do we need the constraints?
if ~in_box
    A = sparse(iA, jA, vA, L^2, U);
    return
end

%-------------------------------------------------------------------------%
% If we get here we need the box constraints
A{1} = sparse(iA, jA, vA, L^2, U);

% Extract the moment matrix of degree d-1
s = find(sum(p,2) <= d-1);
[col,row] = meshgrid(s,s);
s = sub2ind([L, L], row, col);
B = A{1}(s, :);
[iB, jB, vB] = find(B);

% The localizing matrices, one per variable
% In the chebyshev basis, the constraints are  0.5*T0(x_i) - 0.5*T2(x_i)
r = multi_indices(n, 2*d-2);
M = nchoosek(n+d-1, n);
hash = rand(n,1);
qhash = q*hash;
for i = 1:n
    % Component coming from 0.5*T0(x_i)
    A{1+i} = 0.5*B;
    % Components coming from -0.5*T2(x_i)
    % There are two of these, each with 0.5 weight
    % This means we need to multiply B by 0.25
    T = zeros(1,n);
    T(i) = 2;
    alph = r + T;
    beta = abs(r - T);
    [~, col_alph] = ismembertol(alph*hash, qhash, 1e-12);
    [~, col_beta] = ismembertol(beta*hash, qhash, 1e-12);
    A{1+i} = A{1+i} - sparse(iB, col_alph(jB), 0.25*vB, M^2, U);
    A{1+i} = A{1+i} - sparse(iB, col_beta(jB), 0.25*vB, M^2, U);
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
