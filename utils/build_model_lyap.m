function [model] = build_model_lyap(data, lf, param)

%% Build a stable model for a given lyapunov function
% NOTE: Cancellations of potential odd highest-order terms are
% automatically accounted for when we set p(x) = sos(x) and the sos
% polynomial has a lower degree than p(x).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract stuff
dv = lf.degree;
df = param.df;
nx = param.nx;

% Thresholding tolerance for SINDY-like sparsification
cleanTol = param.sindy.tol;

% initialize empty model and set initial params
model = model_initialize();
[model.basis{1:nx,1}] = deal( multi_indices(nx, df) );
model.coeff = cell(nx,1);
model.num_prob = false;
model.num_data = data.n;
model.degree = df;
model.dim = nx;
model.rescale.Lambda = param.Lambda;
model.rescale.mu = param.mu;

% Loop untile done
done = false;
iter = 0;
while ~done && iter < param.sindy.maxIter

    % clear yalmip (bad but OK for our intended use)
    iter = iter + 1;
    yalmip clear
    sdpvar b c b_slack c_slack mse_slack

    % The polynomial model
    % We set up a basis for each coordinate
    fc = cell(nx,1);
    fv = zeros(data.n, nx, 'like', b);
    df = 0;
    for i = 1:nx
        df = max(df, max(sum(model.basis{i},2)));
        nf = size(model.basis{i}, 1);
        fc{i} = sdpvar(nf,1);
        fv(:,i) = cheb_evaluate(model.basis{i}, data.x) * fc{i};
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% CONSTRAINTS ON b AND c
    cnstr = [c >= 0; c <= 1; ...
             param.beta*c + b >= 0; ...
             param.beta*c - b >= 0; ...
             c_slack + (c - 1) >= 0; ...
             c_slack - (c - 1) >= 0; ...
             b_slack + b >= 0; ...
             b_slack - b >= 0];


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% THE OBJECTIVE
    % We handle the quadratic penalties with a rotated second-order cone
    % constraint. Yalmip CANNOT dualize this.
    mse = reshape( fv - data.y, [], 1 );
    cnstr = [cnstr; cone(mse, mse_slack)];        % sqrt of norms in the cost
    obj = mse_slack + param.eps(1)*c_slack + param.eps(3)*b_slack;

    % slacks for l1 penalty on model coefficients
    for i = 1:numel(fc)
        slack{i} = sdpvar(size(fc{i},1),1);
        cnstr = [cnstr; slack{i}>=fc{i}; fc{i}>=-slack{i}];
        obj = obj + param.eps(2) * sum(slack{i}(:));
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% THE LYAPUNOV INEQUALITY CONSTRAINT
    % Create a multivariate chebyshev basis to represent the polynomial
    % b - c*v - alpha * ( f * grad(v) )
    dlyap = max([1, dv, df + dv - 1]);        % the correct polynomial degree
    basis_lyap = multi_indices(nx, dlyap);    % the basis
    n_extended = size(basis_lyap, 1);         % the size of the basis
    LYAP = zeros(n_extended, 1, 'like', fc{1});

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
        B = reshape(coeff_gradv(:,k) * fc{k}', [], 1);
        LYAP = LYAP + param.alpha * (At * B);
    end
    LYAP = - LYAP; % fix the sign (dumb toolbox)

    % Add the polynomial b - c*v(x) to the Lyapunov expression, expanded
    % with respect to the multivariate chebyshev basis basis_lyap
    [~, i1]  = ismember(zeros(1,nx), basis_lyap, 'rows');
    LYAP(i1) = LYAP(i1) + b;
    [~, iv]  = ismember(lf.basis, basis_lyap, 'rows');
    LYAP(iv) = LYAP(iv) - c * lf.coeff;

    % We have constructed b - c*v - alpha * ( f * grad(v) ), which is to be SOS
    % We now construct the SOS polynomial, written in a gram matrix representation
    ds = 2* floor(0.5*dlyap);
    Qs = sdpvar(nchoosek(nx + 0.5*ds, nx));   % PSD matrix for s

    % Match p = b - c*v - alpha * ( f * grad(v) ) to the SOS polynomial
    % Note that p is expressed in the basis basis_lyap, while the sos
    % polynomial is expressed in the basis bases.upper. Both are multivariate
    % chebyshev bases and baes.upper \subseteq basis_lyap by construction
    % We take this into account below.
    [As, bases] = cheb_mommat(nx, 0.5*ds);
    sos_poly = As'*reshape(Qs,[],1);
    [~, isos]  = ismember(bases.upper, basis_lyap, 'rows');
    LYAP(isos) = LYAP(isos) - sos_poly;

    % The actual SOS constraint: matching coefficients and PSD Gram matrix in
    % the SOS polynomial. This automatically sets up the right cancellation if
    % the degree of kappa - c*v - param.alpha * ( f * grad(v) ) is odd
    cnstr = [cnstr; LYAP==0; Qs>=0];


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% OPTIMIZE (finally!)
    yalmip_opts = sdpsettings();
    yalmip_opts.dualize = 1;
    yalmip_opts.verbose = param.verbose;
    yalmip_opts.solver = param.solver;
    model.opt.sol = optimize(cnstr, obj, yalmip_opts);
    model.opt.mse = value(mse).^2; % term in cost function is sqrt(sum(mse))
    model.opt.b = value(b);
    model.opt.c = value(c);
    model.num_prob = model.num_prob || ( model.opt.sol.problem ~= 0 );

    % Check the solution and remove small coefficients
    % If we remove stuff, we are not actually done!
    done = true;
    for i = 1:nx
        model.coeff{i} = value(fc{i});
        tol = cleanTol * max(abs(model.coeff{i}));
        keep = abs(model.coeff{i}) >= tol;
        if ~all(keep)
            done = false;
            model.basis{i} = model.basis{i}(keep,:);
            model.coeff{i} = model.coeff{i}(keep,:);
        end
    end
end

% Update model output: fix degree and numerical problem
model.degree = 0;
for i = 1:param.nx
    di = max(sum(model.basis{i},2));
    model.degree = max( model.degree, di );
end
model.num_prob = double(model.num_prob);