function lf = build_initial_lyap_constrained(data, param)

% Find an approximate Lyapunov function
% We use EDMD to approximate the Lie derivative
% NOTE: We assume the data has been rescaled to [-1,1]^n

% We clear yalmip (bad practice but OK for our intended use)
yalmip clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract stuff
dv = param.dv;
df = param.df;
nx = param.nx;

% Initialize the lyapunov function object
lf = lyap_initialize();
lf.num_data = data.n;

% Intialize a bivariate chebyshev basis for the Lyapunov function v(x)
basis_v = multi_indices(nx, dv);
nv = size(basis_v, 1);

% Multivariate chebyshev basis w(x) for generator EDMD
dw = dv - 1 + df;
L = zeros(nv, nchoosek(dw + nx, nx));
A = zeros(data.n, nv);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% EDMD
% Are we performing weighted EDMD? (https://arxiv.org/pdf/2511.17772)
if param.wtEDMD
    wt = zeros(data.n,1);
    wt(1) = 0;
    wt(2:end-1) = 142.250375*exp(-1./((2:data.n-1)/(data.n).*(1 - (2:data.n-1)/(data.n)) ));
    wt(end) = 0;
    wt = sqrt(wt);  % column vector of sqrt weights
end

% Actually run EDMD
% NOTE: since we know we want to learn a polynomial model of degree df, we
% can fit polynomials of the right degree to the Lie derivative.
basis_w = multi_indices(nx, dw);
for i = 1:nv
    % Construct and evaluate the local w basis
    dw_loc = sum(basis_v(i,:)) + df - 1;
    basis_w_loc = multi_indices(nx, dw_loc);
    [~, iw] = ismember(basis_w_loc, basis_w, 'rows');
    B = cheb_evaluate(basis_w_loc, data.x);
    % Evaluate the gradient of the v basis element
    % We use CHEBFUN to keep this simple
    grad_v = zeros(data.n, nx);
    for j = 1:nx
        idx = [1:j-1, j+1:nx];
        P1 = cheb_evaluate(basis_v(i, idx), data.x(:,[1:j-1, j+1:nx]));
        dT = diff( chebpoly(basis_v(i,j)) ); % thanks chebfun
        grad_v(:,j) = dT(data.x(:,j)) .* P1;   % thanks again
    end
    A(:,i) = sum(grad_v .* data.y, 2);

    % Weighted EDMD?
    if param.wtEDMD == true
        B = B .* wt;   % each row k multiplied by sqrt(wt_k)
        A(:,i) = A(:,i) .* wt; % same weights for a
    end

    % Compute
    L(i,iw) = (A(:,i).' * B) / (B.' * B + param.EDMDreg*eye(size(B,2)) );
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% OPTIMIZATION
% Set up the optimization variables
dlyap = max([1, dv, dw]);                 % degree of polynomial for lyapunov condition
ds = 2*floor(0.5*dlyap);                  % degree of SOS polynomial s
Qs = sdpvar(nchoosek(nx + 0.5*ds, nx));   % PSD matrix for s
Qc = sdpvar(nchoosek(nx + 0.5*dv, nx));   % PSD matrix for coercivity
vc = sdpvar(nv, 1);                       % coefficients of v
sdpvar b;                                 % constant in Lyapunov inequality
sdpvar obj;                               % slack variable for objective

% Initial constraints: the Gram matrices are PSD
cnstr = [Qs>=0; Qc>=0];

% Construct the coefficients of the SOS polynomial s
[As, bases] = cheb_mommat(nx, 0.5*ds);
sc = As'*reshape(Qs,[],1);

% Set up the polynomial for the Lyapunov inequality 
% NOTE: basis_lyap will certainly contain basis_v and basis_w
basis_lyap = multi_indices(nx, dlyap);
[~,iz] = ismember(zeros(1,nx), basis_lyap, 'rows');
[~,iv] = ismember(basis_v, basis_lyap, 'rows');
[~,iw] = ismember(basis_w, basis_lyap, 'rows');
[~,is] = ismember(bases.upper, basis_lyap, 'rows');

% Initalize the polynomial as the constant b
LYAP = zeros(size(basis_lyap,1), 1, 'like', b);
LYAP(iz) = b;

% Add coefficients of -vc*v - alpha*vc'L*w
LYAP(iw) = LYAP(iw) - (param.alpha .* L)' * vc;
LYAP(iv) = LYAP(iv) - vc;

% We subtract the SOS polynomial
LYAP(is) = LYAP(is) - sc;

% We set the SOS constraint
cnstr = [cnstr; LYAP==0];

% We create the cost function |b|
% This is handled using slack variables
cnstr = [cnstr; obj + b >=0; obj - b >= 0];

% Penalty on lyapunov function coefficients
% vc_slack = sdpvar(nv, 1);
% obj = obj + param.eps(4)/param(eps(3))*sum(vc_slack);
% cnstr = [cnstr; vc_slack + vc >= 0; vc_slack - vc >= 0];

% Finally, the SOS constraint on v - ee*|x|^2 where x is the rescaled
% coordinate. We assume WLOG that v has degree \geq 2 otherwise no hope
% We compute the coefficients of v - ee*|x|^2 w.r.t the multivariate
% chebyshev basis used to construct v. For this, we note that
%   x^2 = 0.5 * T_0(x) + 0.5 * T_2(x)
% where T_i is the i-th Chebyshev polynomial.
coerc = vc;
for i = 1:nx
    pow = zeros(2,nx);  % T0(xi_i)
    pow(2,i) = 2;       % T2(xi_i)
    [~, idx] = ismember(pow, basis_v, 'rows');
    coerc(idx) = coerc(idx) - 0.5*param.eps_coerc;
end
[Ac, bases] = cheb_mommat(nx, 0.5*dv);
sos_poly_c = Ac'*reshape(Qc,[],1);          
[~, isos]  = ismember(bases.upper, basis_v, 'rows');
coerc(isos) = coerc(isos) - sos_poly_c;
cnstr = [ cnstr; coerc==0];

% Solve (finally)
yalmip_opts = sdpsettings('dualize', 1);
yalmip_opts.verbose = param.verbose;
yalmip_opts.solver = param.solver;
lf.yalmip = optimize(cnstr, obj, yalmip_opts);
lf.obj = value(obj)^2;
lf.b = value(b);
lf.c = 1;
lf.coeff = value(vc);
lf.basis = basis_v;
lf.feasible = double( ~ismember(lf.yalmip.problem, [1, 2]) );
lf.num_prob = double( lf.yalmip.problem ~= 0 );
lf.degree = max(sum(lf.basis,2));
end