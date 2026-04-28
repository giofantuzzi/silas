function fode = exact_GuckenheimerHolmes(xyz)

% parameters
a = 0.4;
b = 20.25;
c = 3;
d = 1.6;
e = 1.7;
f = 0.44;

% exact dynnamics
fode = [0; 0; 0];
fode(1) = a * xyz(1) - b * xyz(2) + c * xyz(3) * xyz(1) + d * xyz(3) * xyz(1)^2 + d * xyz(3) * xyz(2)^2;
fode(2) = a * xyz(2) + b * xyz(1) + c * xyz(3) * xyz(2);
fode(3) = e - xyz(3)^2 - f * xyz(1)^2 - f * xyz(2)^2 - a * xyz(3)^3;