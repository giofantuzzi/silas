function param = make_parameters()

% Make structure with default parameter values
% Parameter names are as in the paper. NOTE: param.eps must be a 4-by-M array
param.dv = 2;               % degree of Lyapunov dictionary *MUST BE EVEN and >=2*
param.df = 3;               % degree of model *>=1*
param.Lambda = [];          % invertible matrix, must be specified by the user or determined from data
param.mu = [];              % vector, must be specified by the user or determined from data
param.eps(1,1) = 1e-2;      % penalty for |c-1| in first model learning problem
param.eps(2,1) = 1e-2;      % penalty for \|F\|_1 in model learning problems
param.eps(3,1) = 1e-2;      % penalty for |b| in model learning problem (refinement)
param.eps(4,1) = 1e-6;      % penalty parameter for |u|_1 in Lyapunov function learning (refinement)
param.alpha = 1e+0;         % parameter in the Lyapunov inequality
param.beta  = 1e+6;         % bound to force b=0 if c=0
param.kappa = 1e-1;         % penalty parameter for SOS mismatch
param.eps_coerc = 1e-2;     % coercivity parameter for Lyapunov function (could be absorbed in Lambda and mu, but better kept separate)
param.maxIter = 1;          % number of model training / trapping region training iterations

% Other parameters
param.wtEDMD = true;        % Use weighted EDMD?
param.EDMDreg = 1e-8;       % Regularization for EDMD
param.sindy.tol = 0;        % tolerance for sindy-like cleaning in model learning
param.sindy.maxIter = 10;   % max number of sindy-like iterations
param.nx = [];              % to be specified by the user / determined from data
param.noise = 0;            % fraction of the data's standard deviation for "measurement noise"
param.verbose = 0;          % run optimization in silent mode or not?
param.solver = [];          % conic solver (use default if empty)
