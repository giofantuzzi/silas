function f = exact_SprottA(x)

% vector field
f = [0; 0; 0];
f(1) = x(2);
f(2) = -x(1) + x(2) * x(3);
f(3) = 1 - x(2)^2;
end