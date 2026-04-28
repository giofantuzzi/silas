function f = exact_Thomas(x)

% Parameters from dtsys
a = 1.85;
b = 10;

% vector field
f = [0; 0; 0];
f(1) = -a * x(1) + b * sin(x(2));
f(2) = -a * x(2) + b * sin(x(3));
f(3) = -a * x(3) + b * sin(x(1));

end