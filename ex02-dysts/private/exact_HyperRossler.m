function fode = exact_HyperRossler(xyzw)

% parameters
a = 0.25;
b = 3.0;
c = 0.5;
d = 0.05;

% dynamics
fode = [0; 0; 0; 0];
fode(1) = -xyzw(2) - xyzw(3);
fode(2) = xyzw(1) + a * xyzw(2) + xyzw(4);
fode(3) = b + xyzw(1) * xyzw(3);
fode(4) = -c * xyzw(3) + d * xyzw(4);