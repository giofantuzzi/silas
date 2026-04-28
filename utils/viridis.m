
function cmap = viridis(n)

if nargin < 1
    n = 256;
end

% Anchor points
cm_data = [
    68, 1, 84
    71, 44, 122
    59, 81, 139
    44, 113, 142
    33, 144, 141
    39, 173, 129
    92, 200, 99
    170, 220, 50
    253, 231, 37
] / 255;

% Interpolate smoothly
x = linspace(0,1,size(cm_data,1));
xi = linspace(0,1,n);

cmap = interp1(x, cm_data, xi, 'pchip');
end