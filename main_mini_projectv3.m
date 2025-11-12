
% Clear everything
clear; clc; close all;
iterations=800;
depth_size=100;
C_vector = zeros(iterations,depth_size);
C_vector(1,1:2) = [2e19 2e19];



%% Constant Source Diffusion
% Simplify using the assumtion 
% DΔt/(Δx^2)=1/2
aCSD=0.5;

% C_vectorCSD=C_vector;
% for i = 2:iterations
%     for j = 2:depth_size %size of array -1
% 
%         if j==depth_size
%             C_vectorCSD(i,j)=aCSD*(C_vectorCSD(i-1,j)+C_vectorCSD(i-1,j-1));
%         else
%             C_vectorCSD(i,j)=aCSD*(C_vectorCSD(i-1,j-1)+C_vectorCSD(i-1,j+1));
%         end
%     end
%     C_vectorCSD(i,1:2) = [2e19 2e19];
% 
% 
% end





% figure;
% plot(C_vectorCSD(end, :));
% hold on;

C_vectorCSD=C_vector;
for i = 2:iterations
    for j = 2:depth_size %size of array -1

        if j==depth_size
            C_vectorCSD(i,j)=aCSD*(C_vectorCSD(i-1,j-1));
        else
            C_vectorCSD(i,j)=aCSD*(C_vectorCSD(i-1,j-1)+C_vectorCSD(i-1,j+1));
        end
    end
    C_vectorCSD(i,1:2) = [2e19 2e19];


end
plot(C_vectorCSD(end, :),'LineWidth', 3);
xlabel('Depth index');
ylabel('Concentration');
title('Constant Source Diffusion Profile');

grid on;

%% erfx Comparison
% --- PHYSICAL PARAMETERS: set to your case ---
dz = 1e-7;         % [m] grid spacing between depth nodes (example: 100 nm)
D  = 1e-14;        % [m^2/s] diffusion coefficient
Cs = 2e19;         % [units of concentration] surface concentration (from BC)

% Time step mapping from your explicit scheme:
% aCSD = D*dt/dz^2  =>  dt = aCSD * dz^2 / D
dt = aCSD * dz^2 / D;

% ---- Your existing solver ----
C_vectorCSD = C_vector;
for i = 2:iterations
    for j = 2:depth_size
        if j == depth_size
            C_vectorCSD(i, j) = aCSD * (C_vectorCSD(i-1, j-1));
        else
            C_vectorCSD(i, j) = aCSD * (C_vectorCSD(i-1, j-1) + C_vectorCSD(i-1, j+1));
        end
    end
    C_vectorCSD(i, 1:2) = [Cs Cs];   % keep BC in sync with Cs
end

% ---- Build x-grid and analytical erfc profile at t = iterations*dt ----
x  = (0:depth_size-1) * dz;          % [m] depth array (x=0 at the surface)
T  = iterations * dt;                 % [s] physical time
Cerfc = Cs * erfc( x ./ (2*sqrt(D*T)) );

% ---- Plot numerical vs analytical ----
figure; hold on; grid on;
plot(x, C_vectorCSD(end, :), 'LineWidth', 3, 'DisplayName', 'Numerical (last iter)');
plot(x, Cerfc,               '--', 'LineWidth', 2, 'DisplayName', 'Analytical: Cs·erfc(x/(2√(DT)))');

xlabel('Depth x [m]');
ylabel('Concentration');
title(sprintf('Constant Source Diffusion: Numerical vs. erfc at t = %.3g s', T));
legend('Location','best');
hold off;

 %% LOg version
% C_vectorCSD = C_vector;
% 
% for i = 2:iterations
%     for j = 2:depth_size
%         if j == depth_size
%             C_vectorCSD(i, j) = aCSD * (C_vectorCSD(i-1, j-1));
%         else
%             C_vectorCSD(i, j) = aCSD * (C_vectorCSD(i-1, j-1) + C_vectorCSD(i-1, j+1));
%         end
%     end
%     % Boundary conditions
%     C_vectorCSD(i, 1:2) = [2e19 2e19];
% end
% 
% % --- Plot in log scale ---
% semilogy(C_vectorCSD(end, :), 'LineWidth', 3);
% xlabel('Depth index');
% ylabel('Concentration (log scale)');
% title('Constant Source Diffusion Profile (Log Output)');
% grid on;


%% Drive-In Diffusion
% Simplify using the assumtion 
% DΔt/(Δx^2)=1/2
aDID=0.5;
iterations=1000;
C_vectorDID=zeros(iterations,depth_size);
C_vectorDID(1, :)=C_vectorCSD(end, :);

for i = 2:iterations

    for j = 1:depth_size %size of array -1





        if j==depth_size
            C_vectorDID(i,j)=aDID*(C_vectorDID(i-1,j)+C_vectorDID(i-1,j-1));
        elseif j==1
            C_vectorDID(i,j)=aDID*(C_vectorDID(i-1,j)+C_vectorDID(i-1,j+1));
        else
            C_vectorDID(i,j)=aDID*(C_vectorDID(i-1,j-1)+C_vectorDID(i-1,j+1));
        end
    end



end
% Plot: Constant source diffusion
figure;
plot(C_vectorDID(end, :),'LineWidth', 3);
xlabel('Depth index');
ylabel('Concentration');
title('Final Diffusion Profile');
grid on;

%% Gaussian Comparison

%% Parameters you must set (consistent units)
dz = 1e-7;       % [m] grid spacing (e.g., 100 nm)
D  = 1e-14;      % [m^2/s] diffusivity for drive-in step

aDID = 0.5;                % you already chose this
iterations = 1000;         % you already chose this
dt = aDID * dz^2 / D;      % scheme-to-physics mapping
T  = iterations * dt;      % physical time of drive-in

%% Drive-In Diffusion (your code as given)
C_vectorDID = zeros(iterations, depth_size);
C_vectorDID(1, :) = C_vectorCSD(end, :);   % initial (from pre-deposition)

for i = 2:iterations
    for j = 1:depth_size
        if j == depth_size
            C_vectorDID(i, j) = aDID * (C_vectorDID(i-1, j) + C_vectorDID(i-1, j-1));
        elseif j == 1
            C_vectorDID(i, j) = aDID * (C_vectorDID(i-1, j) + C_vectorDID(i-1, j+1));
        else
            C_vectorDID(i, j) = aDID * (C_vectorDID(i-1, j-1) + C_vectorDID(i-1, j+1));
        end
    end
end

%% Build x-grid and Gaussian analytic for finite dose
x   = (0:depth_size-1) * dz;                % depth [m], surface at x=0
Q   = trapz(x, C_vectorDID(1, :));          % total dose from initial profile
Cgauss = (Q ./ sqrt(pi * D * T)) .* exp(-(x.^2) ./ (4 * D * T));

%% Plot: numerical vs Gaussian
figure; hold on; grid on;
plot(x, C_vectorDID(end, :), 'LineWidth', 3, 'DisplayName', 'Numerical (drive-in)');
plot(x, Cgauss, '--', 'LineWidth', 2, 'DisplayName', 'Gaussian: Q/√(πDt)·exp(-x^2/4Dt)');
xlabel('Depth x [m]');
ylabel('Concentration');
title(sprintf('Drive-In Diffusion vs Gaussian (t = %.3g s)', T));
legend('Location','best');
hold off;

%% Optional: log-y view (helps compare tails)
% figure; hold on; grid on;
% semilogy(x, C_vectorDID(end, :), 'LineWidth', 3, 'DisplayName', 'Numerical');
% semilogy(x, Cgauss, '--', 'LineWidth', 2, 'DisplayName', 'Gaussian');
% xlabel('Depth x [m]'); ylabel('Concentration (log scale)');
% title(sprintf('Log view (t = %.3g s)', T)); legend('Location','best'); hold off;

%% Quick error metrics
num = C_vectorDID(end, :);
rmse = sqrt(mean((num - Cgauss).^2));
relL1 = trapz(x, abs(num - Cgauss)) / max(Q, eps);
fprintf('RMSE = %.3e,  Relative L1 (vs dose) = %.3e\n', rmse, relL1);
