function lf = build_lyap_model(model, param)

%% Optimize a Lyapunov function giving a trapping region for an ODE model

% We clear yalmip (bad practice but OK for our intended use)
yalmip clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract stuff
dv = param.dv;
nx = param.nx;

% Initialize the lyapunov function object
lf = lyap_initialize();
lf.basis = multi_indices(nx, dv);
nv = size(lf.basis, 1);
lf.coeff = sdpvar(nv, 1);
lf.num_data = NaN;

% Find model degree
df = model.degree;

% extract constant c from model
c = model.opt.c;

% Optimization variables
slack = sdpvar(nv,1); 
sdpvar b b_slack

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% THE OBJECTIVE FUNCTION
% obj = param.eps(3)*|b| + param.eps(4)*sum(slack);
obj = param.eps(3)*b_slack + param.eps(4)*sum(slack);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% THE SIMPLE CONSTRAINTS
cnstr = [slack - lf.coeff >= 0; ...
         slack + lf.coeff >= 0; ...
         param.beta*c - b >= 0; ...
         param.beta*c + b >= 0];
cnstr = [cnstr; b_slack + b >= 0; b_slack - b >= 0];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% THE LYAPUNOV INEQUALITY CONSTRAINT
% Create a multivariate chebyshev basis to represent the polynomial
% b - c*v - alpha * ( f * grad(v) )
dlyap = max([1, dv, df + dv - 1]);        % the correct polynomial degree    
basis_lyap = multi_indices(nx, dlyap);    % the basis
n_extended = size(basis_lyap, 1);         % the size of the basis
LYAP = zeros(n_extended, 1, 'like', b);

% The coefficients and chebyshev basis for f * grad(v)
[basis_gradv, coeff_gradv] = cheb_gradient(lf.basis, lf.coeff);
n_gradv = size(basis_gradv,1);
for k = 1:nx
    fun = @(i,j) cheb_multiply(model.basis{k}(i,:), basis_gradv(j,:));
    nf = size(model.basis{k}, 1);
    [I, J] = meshgrid(1:nf, 1:n_gradv);
    gamma  = arrayfun(fun, I, J, 'UniformOutput',false);
    hash = rand(nx,1);
    gamma_hash = reshape(vertcat(gamma{:})*hash, 2^nx, []);
    qhash = basis_lyap*hash;
    iA = []; jA = []; vA = [];
    for i = 1:n_extended
        [~, ind] = find(abs(gamma_hash - qhash(i))<=1e-12);
        iA = [iA; ind];
        jA = [jA; repmat(i, numel(ind), 1)];
        vA = [vA; repmat(1/2^nx, numel(ind), 1)];
    end
    At = sparse(jA, iA, vA, n_extended, n_gradv*nf);
    B = reshape(coeff_gradv(:,k) * model.coeff{k}', [], 1);
    LYAP = LYAP + param.alpha * (At * B);
end
LYAP = -LYAP; % fix the sign...

% Add the polynomial kappa - v(x) to the Lyapunov expression, expanded
% with respect to the multivariate chebyshev basis basis_lyap
[~, i1]  = ismember(zeros(1,nx), basis_lyap, 'rows');
LYAP(i1) = LYAP(i1) + b;
[~, iv]  = ismember(lf.basis, basis_lyap, 'rows');
LYAP(iv) = LYAP(iv) - c * lf.coeff;

% We have constructed b - c*v - alpha * ( f * grad(v) ), which is to be SOS
% We now construct the SOS polynomial, written in a gram matrix representation
ds = 2* floor(0.5*dlyap);               % admissible degree of the lie derivative
Qs = sdpvar(nchoosek(nx + 0.5*ds, nx));   % PSD matrix for s

% Match p = b - c*v - alpha* ( f * grad(v) ) to the SOS polynomial
% Note that p is expressed in the basis basis_lyap, while the sos
% polynomial is expressed in the basis bases.upper. Both are multivariate
% chebyshev bases and baes.upper \subseteq basis_lyap by construction
% We take this into account below.
[As, bases_s] = cheb_mommat(nx, 0.5*ds);
sos_poly_s = As'*reshape(Qs,[],1);          
[~, isos]  = ismember(bases_s.upper, basis_lyap, 'rows');
LYAP(isos) = LYAP(isos) - sos_poly_s;

% The actual SOS constraint: matching coefficients and PSD Gram matrix in
% the SOS polynomial. This automatically sets up the right cancellation if
% the degree of b - c*v - alpha * ( f * grad(v) ) is odd
cnstr = [cnstr; LYAP == 0; Qs>=0];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% THE COERCIVITY CONSTRAINT
% Finally, the SOS constraint on v - ee*|x|^2 where x is the rescaled
% coordinate. We assume WLOG that v has degree \geq 2 otherwise no hope
% We compute the coefficients of v - ee*|x|^2 w.r.t the multivariate
% chebyshev basis used to construct v. For this, we note that
%   x^2 = 0.5 * T_0(x) + 0.5 * T_2(x)
% where T_i is the i-th Chebyshev polynomial.
coerc = lf.coeff;
for i = 1:nx
    pow = zeros(2,nx);  % T0(xi_i)
    pow(2,i) = 2;       % T2(xi_i)
    [~, idx] = ismember(pow, lf.basis, 'rows');
    coerc(idx) = coerc(idx) - 0.5*param.eps_coerc;
end
[Ac, bases_c] = cheb_mommat(nx, 0.5*dv);
Qc = sdpvar(nchoosek(nx + 0.5*dv, nx));
sos_poly_c = Ac' * reshape(Qc,[],1);          
[~, isos]  = ismember(bases_c.upper, lf.basis, 'rows');
coerc(isos) = coerc(isos) - sos_poly_c;
cnstr = [ cnstr; coerc==0; Qc>=0 ];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% OPTIMIZE
yalmip_opts = sdpsettings();
yalmip_opts.dualize = 1; % cannot dualize rcone constraints...
yalmip_opts.verbose = param.verbose;
yalmip_opts.solver = param.solver;
lf.yalmip = optimize(cnstr, obj, yalmip_opts);
lf.obj = value(obj);
lf.b = value(b);
lf.c = c;
lf.coeff = value(lf.coeff);
lf.feasible = double( lf.yalmip.problem ~= 1 );
lf.num_prob = double( lf.yalmip.problem ~= 0 );
lf.degree = max(sum(lf.basis,2));