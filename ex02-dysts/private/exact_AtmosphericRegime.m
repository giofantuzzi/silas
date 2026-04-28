function fode = exact_AtmosphericRegime(xyz)

% parameters
alpha = -2.0;
beta = -5.0;
mu1 = 0.05;
mu2 = -0.01;
omega = 3.0;
sigma = 1.1;

% exact dynamics
fode = [0;0;0];
fode(1) = mu1 * xyz(1) + sigma * xyz(1) * xyz(2);
fode(2) = mu2 * xyz(2) + omega * xyz(3) + alpha * xyz(2) * xyz(3) + beta * xyz(3)^2 - sigma * xyz(1)^2;
fode(3) = mu2 * xyz(3) - omega * xyz(2) - alpha * xyz(2)^2 - beta * xyz(2) * xyz(3);