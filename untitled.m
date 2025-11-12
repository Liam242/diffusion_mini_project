%% Atomistic Diffusion Model - Finite Difference Implementation
% Based on: mini project brief (Berwick) and Eq. 7.40 (constant D)
clear; clc; close all;

%% Common parameters
N = 100;                 % number of grid points
L = 1e-4;                % total depth (cm), ~1 µm
x = linspace(0, L, N);   % depth vector
dx = x(2) - x(1);

D = 3e-14;               % diffusion coefficient (cm^2/s) - choose reasonable value

% Stability parameter for explicit scheme:
% alpha = D * dt / dx^2  must be <= 0.5 for stability
alpha = 0.4;                     % choose stable value
dt = alpha * dx^2 / D;           % time step from chosen alpha

Cs = 2e19;               % constant surface concentration (cm^-3)
nSteps_const = 800;      % iterations for constant source
nSteps_drive = 800;      % iterations for drive-in

%% 1. CONSTANT SOURCE DIFFUSION (Eq. 7.40, constant D)
% Initial condition: "delta-like" at surface, as in the brief
C_const = zeros(1, N);
C_const(1:2) = Cs;

for n = 1:nSteps_const
    C_new = C_const;
    
    % Left boundary: constant source (Dirichlet)
    C_new(1) = Cs;
    
    % Interior points: explicit finite difference for dC/dt = D d2C/dx2
    for i = 2:N-1
        C_new(i) = C_const(i) + alpha * (C_const(i+1) - 2*C_const(i) + C_const(i-1));
    end
    
    % Right boundary: zero-flux (Neumann) -> dC/dx = 0
    C_new(N) = C_new(N-1);
    
    C_const = C_new;
end

t_const = nSteps_const * dt;

% Analytical constant-source (erfc) for comparison
C_erfc = Cs * erfc( x ./ (2*sqrt(D * t_const)) );

% Plot: Constant source diffusion
figure;
plot(x * 1e4, C_const, 'LineWidth', 2); hold on;
plot(x * 1e4, C_erfc, '--', 'LineWidth', 2);
xlabel('Depth (\mum)');
ylabel('Concentration (cm^{-3})');
title('Constant Source Diffusion: Numerical vs Analytical (erfc)');
legend('Numerical (FD)', 'Analytical erfc', 'Location', 'northeast');
grid on;


%% 2. DRIVE-IN DIFFUSION (from limited source)
% Use the final constant-source profile as the initial "dose" for drive-in.
% Now remove the constant source: total dopant should be (approximately) conserved.
C_drive = C_const;      

% Initial total dose Q (cm^-2)
Q_initial = sum(C_drive) * dx;

for n = 1:nSteps_drive
    C_new = C_drive;
    
    % Left boundary: zero-flux (no new dopant in or out at surface)
    % dC/dx = 0  -> C(1) = C(2)
    C_new(1) = C_drive(2);
    
    % Interior points
    for i = 2:N-1
        C_new(i) = C_drive(i) + alpha * (C_drive(i+1) - 2*C_drive(i) + C_drive(i-1));
    end
    
    % Right boundary: zero-flux
    C_new(N) = C_drive(N-1);
    
    C_drive = C_new;
end

t_drive = nSteps_drive * dt;

% Total dose after drive-in (should be ~constant)
Q_final = sum(C_drive) * dx;

% Analytical Gaussian (limited source) using conserved Q
Q = Q_initial;  % use numerical dose from pre-deposition
C_gauss = (Q ./ sqrt(pi * D * t_drive)) .* exp( - x.^2 ./ (4 * D * t_drive) );

% Plot: Drive-in diffusion
figure;
plot(x * 1e4, C_drive, 'LineWidth', 2); hold on;
plot(x * 1e4, C_gauss, '--', 'LineWidth', 2);
xlabel('Depth (\mum)');
ylabel('Concentration (cm^{-3})');
title('Drive-in Diffusion: Numerical vs Analytical (Gaussian)');
legend('Numerical (FD)', 'Analytical Gaussian', 'Location', 'northeast');
grid on;

% Display dose conservation in Command Window
fprintf('Drive-in diffusion dose check:\n');
fprintf('  Q_initial = %.3e cm^-2\n', Q_initial);
fprintf('  Q_final   = %.3e cm^-2\n', Q_final);


%% 3. INSTABILITY DEMONSTRATION (alpha > 0.5)
% Show that the explicit scheme blows up when D*dt/dx^2 > 0.5

alpha_unstable = 0.6;           % > 0.5 -> unstable
dt_unstable = alpha_unstable * dx^2 / D;

C_unstable = zeros(1, N);
C_unstable(1:2) = Cs;

nSteps_unstable = 80;           % fewer steps: it will blow quickly

for n = 1:nSteps_unstable
    C_new = C_unstable;
    
    % Left boundary: constant source
    C_new(1) = Cs;
    
    % Interior points with unstable alpha
    for i = 2:N-1
        C_new(i) = C_unstable(i) + alpha_unstable * (C_unstable(i+1) - 2*C_unstable(i) + C_unstable(i-1));
    end
    
    % Right boundary: zero-flux
    C_new(N) = C_new(N-1);
    
    C_unstable = C_new;
end

figure;
plot(x * 1e4, C_unstable, 'LineWidth', 2);
xlabel('Depth (\mum)');
ylabel('Concentration (cm^{-3})');
title('Unstable Scheme Example: \alpha > 0.5 (Profile Blows Up)');
grid on;

% The profile will show oscillations/unphysical negative/huge values,
% illustrating numerical instability.

