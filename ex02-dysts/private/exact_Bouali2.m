function fode = exact_Bouali2(xyz)

% parameters
a = 3.0;
b = 2.2;
bb = 0;
c = 0;
g = 1.0;
m = -0.0026667;
y0 = 1.0;

% dynamics
fode = [0; 0; 0];
fode(1) = a * y0 * xyz(1) - a * xyz(1) * xyz(2) - b * xyz(3);
fode(2) = -g * xyz(2) + g * xyz(2) * xyz(1)^2;
fode(3) = -1.5 * m * xyz(1) + m * bb * xyz(1) * xyz(3) - c * xyz(3);