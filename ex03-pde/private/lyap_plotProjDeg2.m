function XY = lyap_plotProjDeg2(lf_final, model, proj_idx, npts)
% lyap_plotProjDeg2
%
% Plot the projection onto coordinates proj_idx of the quadratic Lyapunov
% sublevel set
%
%   { x : V(x) <= level }
%
% where V is stored in lf_final in Chebyshev basis, with degree 2.
%
% INPUTS
%   lf_final   struct with fields:
%                - degree
%                - basis   (N x nx)
%                - coeff   (N x 1)
%                - b, c    (optional, if level not given)
%   model      struct with field model.rescale containing:
%                - mu
%                - Lambda
%   proj_idx   two coordinates to keep, default [1 2]
%   npts       number of boundary points, default 400
%
% OUTPUTS
%   XY         boundary points in ORIGINAL (unscaled) coordinates
%
% Notes
%   - Assumes lf_final is degree 2 in tensor-product Chebyshev basis:
%         T0(x)=1, T1(x)=x, T2(x)=2x^2-1
%   - Assumes projection is taken after minimizing over the remaining
%     scaled coordinates.
%   - Output XY is mapped back to original coordinates using model.rescale.

    if nargin < 3 || isempty(proj_idx)
        proj_idx = [1 2];
    end
    if nargin < 4 || isempty(npts)
        npts = 400;
    end

    assert(lf_final.degree == 2, 'This function only supports lf_final.degree = 2.')

    basis = lf_final.basis;
    coeff = lf_final.coeff(:);
    level = lf_final.b / lf_final.c;
    [N, nx] = size(basis);

    assert(numel(coeff) == N, 'lf_final.coeff and lf_final.basis size mismatch.')
    assert(numel(proj_idx) == 2, 'proj_idx must have length 2.')
    assert(all(proj_idx >= 1 & proj_idx <= nx), 'proj_idx out of range.')
    assert(numel(unique(proj_idx)) == 2, 'proj_idx entries must be distinct.')
    assert(all(sum(basis,2) <= 2), 'This function only supports total degree <= 2.')

    %--------------------------------------------------------------
    % 1) Extract V(x) = x'Qx + 2 q'x + r from Chebyshev representation
    %--------------------------------------------------------------
    Q = zeros(nx, nx);
    q = zeros(nx, 1);
    r = 0;

    for k = 1:N
        alpha = basis(k,:);
        c = coeff(k);
        deg = sum(alpha);
        idx = find(alpha > 0);

        if deg == 0
            % T0...T0 = 1
            r = r + c;

        elseif deg == 1
            % T1(x_i) = x_i
            i = idx(1);
            q(i) = q(i) + c/2;

        elseif deg == 2
            if numel(idx) == 1 && alpha(idx) == 2
                % T2(x_i) = 2 x_i^2 - 1
                i = idx(1);
                Q(i,i) = Q(i,i) + 2*c;
                r = r - c;

            elseif numel(idx) == 2 && all(alpha(idx) == 1)
                % T1(x_i) T1(x_j) = x_i x_j
                i = idx(1);
                j = idx(2);
                Q(i,j) = Q(i,j) + c/2;
                Q(j,i) = Q(j,i) + c/2;

            else
                error('Unexpected degree-2 term in lf_final.basis.')
            end
        else
            error('Only total degree <= 2 is supported.')
        end
    end

    %--------------------------------------------------------------
    % 2) Project analytically by minimizing over remaining coordinates
    %--------------------------------------------------------------
    elim_idx = setdiff(1:nx, proj_idx, 'stable');

    Qzz = Q(proj_idx, proj_idx);
    qz  = q(proj_idx);

    if isempty(elim_idx)
        Qproj = Qzz;
        qproj = qz;
        rproj = r;
    else
        Qzy = Q(proj_idx, elim_idx);
        Qyy = Q(elim_idx, elim_idx);
        qy  = q(elim_idx);

        % Need bounded minimization in eliminated variables
        [~, pchol] = chol(0.5*(Qyy+Qyy.'), 'lower');
        assert(pchol == 0, 'Eliminated block Qyy is not positive definite.')

        Qproj = Qzz - Qzy * (Qyy \ Qzy.');
        qproj = qz  - Qzy * (Qyy \ qy);
        rproj = r   - qy.' * (Qyy \ qy);
    end

    %--------------------------------------------------------------
    % 3) Parameterize projected ellipse in scaled coordinates
    %--------------------------------------------------------------
    [~, pchol] = chol(Qproj, 'lower');
    assert(pchol == 0, 'Projected quadratic is not positive definite.')

    ctr = -Qproj \ qproj;
    rho = level - rproj + qproj.' * (Qproj \ qproj);

    assert(rho >= -1e-12, 'Projected level set is empty.')

    if rho < 0
        rho = 0;
    end

    [U, D] = eig(Qproj);
    lam = diag(D);
    assert(all(lam > 0), 'Projected quadratic is not positive definite.')

    th = linspace(0, 2*pi, npts).';
    C = [cos(th), sin(th)].';          % 2 x npts
    A = U * diag(sqrt(rho ./ lam));    % 2 x 2
    XY_scaled = (ctr + A*C).';         % npts x 2

    %--------------------------------------------------------------
    % 4) Map projected coordinates back to original variables
    %--------------------------------------------------------------
    if nargin >= 2 && isstruct(model) && isfield(model, 'rescale') ...
            && isfield(model.rescale, 'mu') && isfield(model.rescale, 'Lambda')

        mu = model.rescale.mu(:);

        % Build the needed 2x2 diagonal block
        if isvector(model.rescale.Lambda)
            Lambda_proj = diag(model.rescale.Lambda(proj_idx));
        else
            Lambda_proj = model.rescale.Lambda(proj_idx, proj_idx);
        end

        mu_proj = mu(proj_idx);

        % Same map as before: x_phys = Lambda^{-1} (x_scaled - mu)
        XY = (Lambda_proj \ (XY_scaled.' - mu_proj)).';
    else
        XY = XY_scaled;
    end

end