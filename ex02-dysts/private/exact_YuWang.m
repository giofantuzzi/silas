function f = exact_YuWang(x)

% parameters
a = 10;
b = 40;
c = 2;
d = 2.5;

% exact dynamics
f = [0; 0; 0];
f(1) = a * (x(2) - x(1));
f(2) = b * x(1) - c * x(1) * x(3);
f(3) = exp(x(1) * x(2)) - d * x(3);

end