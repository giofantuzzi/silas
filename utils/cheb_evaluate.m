function B = cheb_evaluate(basis, pts)

% Evaluate multivariate chebyshev specified by basis matrix at the points
% given by the matrix pts. Each column of pts gives a coordinate of the point.
% Each column of B corresponds to evaluation of a chebyshev basis element.

% Initialize
[U,n] = size(pts);
L = size(basis,1);
B = ones(U,L);

% Compute a table of all Chebyshev values that we will need
dmax = max(basis,[],'all');
pts = reshape(pts, [], 1);
T = chebpolyval(fliplr(eye(dmax+1)),pts);
idx = reshape(1:U*n, U, n);

% Now operate 
for col = 1:L
    for j = 1:n
        dj = basis(col,j);
        B(:,col) = B(:,col).*T(idx(:,j),dj+1);
    end
end
end

% Clean
% B(abs(B)<1e-14) = 0;