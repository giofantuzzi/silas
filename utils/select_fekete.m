function [pts, wts] = select_fekete(basis, pts)
% Select points from the matrix pts to have good interpolation for the
% multivariate chebyshev basis specified by the input "basis"

% parameters
[nrPoints,n] = size(pts);
U = size(basis,1);

% Evaluate chebyshev basis at initial points
P = ones(nrPoints,U);
m = ones(U,1);
for i=1:U
    bi = basis(i,:);
    for j = 1:n
        coef = zeros(bi(j)+1,1); 
        coef(1) = 1;
        P(:,i) = P(:,i).*chebpolyval(coef,pts(:,j));
        if bi(j) == 1; m(i) = 0;
        else; m(i) = m(i)*(((-1)^bi(j)+1)/(1-bi(j)^2)); end
    end
end

% extracts the positive entries of w and the corresponding points
wts = P'\ rand(U,1);%  m;
ind = abs(wts)>0;
wts = wts(ind);
pts = pts(ind,:);
return