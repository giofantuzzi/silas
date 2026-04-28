function lf = lyap_initialize()

% Initialize an empty model as a structure with all of the fields we will
% need throughout the computation.
lf.degree = [];
lf.basis = [];
lf.coeff = [];
lf.b = [];          % absorbing level set is v <= b/c
lf.c = [];          % absorbing level set is v <= b/c
lf.num_prob = [];
lf.yalmip = [];
lf.others = [];
lf.num_data = [];