function lf = lyap_load(model_name)

% Load lyapunov function from hdf5 file

% Basic lf info
lf = lyap_initialize();
lf.b = h5read(model_name, '/lf_final/b');
lf.c = h5read(model_name, '/lf_final/c');
lf.degree = h5read(model_name, '/lf_final/degree');

% Model basis and coefficients
lf.basis = h5read(model_name, '/lf_final/basis');
lf.coeff = h5read(model_name, '/lf_final/coeff');


