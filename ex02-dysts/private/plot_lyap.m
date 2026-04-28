function XY = plot_lyap(model_name, model_dir, nPoints)

% plot lyapunov function in 3D, projected on the (x,z) plane
% We proceed by finding on a ray the largest point for which there exists
% y such that v(x,y,z) < level.

if nargin < 3; nPoints = 20; end

% Load lyapunov function
% model_dir = 'models-v2';
% model_name = '03_AtmosphericRegime.hdf5';
lf = lyap_load([model_dir, filesep, model_name]);
model = model_load([model_dir, filesep, model_name]);

% check that it is a 3D system o a bivariate quartic
nx = size(lf.basis, 2);
isok = nx==3;
isok = isok || (nx==4 && lf.degree<=4);
assert(isok, 'Must have trivariate or 4-variate quartic')

% Level set pf Lyapunov function
level = lf.b / lf.c;

% Ensure the Lyapunov function has a constant term, even if it has zero
% coefficient. We will need this later
lf.basis = [zeros(1,nx); lf.basis];
lf.coeff = [0; lf.coeff];

% Loop over angular coordinates to find boundary of projection of trapping
% region. This will take a while because we are solving a number of
% univariate SOS problems, which are fast but we go through yalmip every
% time
th = linspace(0,2*pi,nPoints)';
r = NaN(size(th));
r0 = 1;
for i = 1:numel(th)-1
    if i > 1; r0 = r(i-1); end
    fun = @(r) compute_vmin(r, th(i), lf) - level;
    r(i) = fzero(fun,r0);
end
r(end) = r(1);

% Get the boundary of the projection
XY = zeros(numel(th),nx);
XY(:,[1 3]) = [r.*cos(th), r.*sin(th)];
XY = ( XY - model.rescale.mu ) / model.rescale.Lambda';
XY = XY(:,[1 3]);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vmin = compute_vmin(r, th, lf)

% Find x, z coordinates
x = r*cos(th);
z = r*sin(th);

% coefficients of polynomial in y for fixed y, z
% This has repetition, we will remove it later
C = lf. coeff .* cheb_evaluate(lf.basis(:,[1 3]), [x,z])';

% minimize this polynomial over y -- we can look for the largest SOS
% lower bound because it is univariate. We can use the chebyshev basis
% code we already have for this.
nx = size(lf.basis, 2);
varID = setdiff(1:nx, [1 3]);
[v_basis, ~, idx] = unique(lf.basis(:,varID), 'rows'); % this is sorted by default
C = accumarray(idx, C, [size(v_basis,1), 1]);     % remove the repetions
dv = max(sum(v_basis, 2));
nx = numel(varID);

if nx==1 && dv == 2
    % we are minimizing a univariate quadratic, which ought to be
    % coercive by construction. We compute its minimum analytically using
    % the fact that
    % p(y) = C(1) + C(2)*x + C(3)*(2x^2 - 1)
    %      = [C(1) - C(3)] + C(2)*x + 2C(3)*x^2
    y = -C(2) / (2 * C(3));
    vmin = C(1) - C(3) + C(2)*y + 2*C(3)*y^2;

else
    % SOS problem (we should handle quadratics separately but OK for now)
    sdpvar t
    d = floor( dv / 2 );
    L = nchoosek(nx + d, d);
    Q = sdpvar(L, L);
    [A, sos_bases] = cheb_mommat(nx, d);
    sos_poly = A' * reshape(Q,[],1);
    sos_coeff = zeros(size(v_basis,1), 1, 'like', t);
    [~, isos]  = ismember(zeros(1,nx), v_basis, 'rows');
    sos_coeff(isos) = t;
    sos_coeff = C - sos_coeff;
    [~, isos]  = ismember(sos_bases.upper, v_basis, 'rows');
    sos_coeff(isos) = sos_coeff(isos) - sos_poly;
    cnstr = [sos_coeff == 0; Q>=0];
    optimize(cnstr, -t, sdpsettings('dualize', 1, 'verbose', 0, 'cachesolvers',1));

    % Get value of v
    vmin = value(t);    
end

end