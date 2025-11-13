% Define x-axis
x = linspace(-10, 10, 1000);

% Parameters for Gaussian
sigma = 1;      % width
mu = 0;         % center

% Gaussian curve
G = gaussmf(x, [sigma mu]);

% Plot
figure;
plot(x, G, 'LineWidth', 2);
grid on;

xlabel('x');
ylabel('Amplitude');
title('Gaussian Curve Using gaussmf');
