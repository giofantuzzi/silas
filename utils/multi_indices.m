function A = multi_indices(n, d)

% Create all n-dimensional multi-indices of length <= d.
% EXTERNAL FUNCTIONS:
% * partitions
L = nchoosek(n+d,  d);
A = ones(L,n);
shift = 0;
for t = 0:d
    allDegs = partitions(t, ones(1,n));
    nrDegs = size(allDegs,1);
    A(shift+1:shift+nrDegs, :) = allDegs;
    shift = shift + nrDegs;
end