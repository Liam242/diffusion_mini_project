clear; clc; close all;

% the main equation we should be focusing on is  DΔt/(Δx^2)=1/2



%% Initialization
% Constant Source Diffusion
% We make DΔt/(Δx^2)=1/2
% T=1200K
% R=1.987 (Universal Gas Constant)
%D=5.275027924x10^-19m^2/s
% Δt= 947867.2988 s
% Δx = 1x10^-7 m
D=5.275027924e-19;
deltaT = 9478.672988; % Time step
deltaX = 1e-7; % Spatial step

% D cant' change aslong as temp doesn't change, but what we can 
%change is the time step. And make it more accurate
% Not that deltaT is initially wrong, just that whe want it to be exactly
% 0.5
deltaT= (0.5*(deltaX^2))/(D);

a=(D*deltaT)/(deltaX^2);
% a=0.5;

%% Answering Question 1

iterations=800;
depth_size=100;
C_vector = zeros(iterations,depth_size);
C_vector(1,1:2) = [2e19 2e19];

C_vectorCSD=C_vector;
for i = 2:iterations
    for j = 2:depth_size %size of array -1

        if j==depth_size
            C_vectorCSD(i,j)=a*(C_vectorCSD(i-1,j-1));
        else
            C_vectorCSD(i,j)=a*(C_vectorCSD(i-1,j-1)+C_vectorCSD(i-1,j+1));
        end
    end
    C_vectorCSD(i,1:2)=[2e19 2e19];


end
figure; 
hold on; 
grid on;
x=(0:depth_size-1)*deltaX;
plot(x, C_vectorCSD(end,:), 'LineWidth', 3);
xlabel('Depth (m)');
ylabel('Concentration (cm^-3)');
title('Constant Source Diffusion Profile at 1200 K');
grid off;

fprintf("Constant Source Diffusion Time Taken = %.3e s\n", (iterations * deltaT)/(3600*24));

%% 1.2 Showing that as time increases, so does concentration versus depth



figure;
hold on;
x=(0:depth_size-1)*deltaX;

for k = 1:10
    iterations = k * 100;        % Number of iterations
    C_vector = zeros(iterations, depth_size);
    C_vector(1, 1:2) = [2e19 2e19];

    % Initialize C_vectorCSD
    C_vectorCSD = C_vector;

    % Time iteration
    for i = 2:iterations
        for j = 2:depth_size
            if j == depth_size
                C_vectorCSD(i, j) = a * (C_vectorCSD(i-1, j-1));
            else
                C_vectorCSD(i, j) = a * (C_vectorCSD(i-1, j-1) + C_vectorCSD(i-1, j+1));
            end
        end
        % Apply boundary condition at the start of each iteration
        C_vectorCSD(i, 1:2) = [2e19 2e19];
    end

    % Plot the last iteration for each case
    plot(x,C_vectorCSD(end, :), 'DisplayName', sprintf('%d days', (iterations * deltaT)/(3600*24)));
end

xlabel('Depth Index');
ylabel('Concentration (C)');
title('Concentration Profile for Different Time Amounts');
legend show;
grid on;
hold off;

%% Answering Question 2


iterations=800;
depth_size=100;
C_vector = zeros(iterations,depth_size);
C_vector(1,1:2) = [2e19 2e19];

C_vectorCSD=C_vector;
for i = 2:iterations
    for j = 2:depth_size %size of array -1

        if j==depth_size
            C_vectorCSD(i,j)=a*(C_vectorCSD(i-1,j-1));
        else
            C_vectorCSD(i,j)=a*(C_vectorCSD(i-1,j-1)+C_vectorCSD(i-1,j+1));
        end
    end
    C_vectorCSD(i,1:2)=[2e19 2e19];


end
figure; 
hold on; 
grid on;

x=(0:depth_size-1)*deltaX;
Cerfc=(2e19)*erfc(x ./ (2*sqrt(D*deltaT*iterations)));
plot(x, C_vectorCSD(end, :), 'LineWidth', 3, 'DisplayName', 'Numerical');
plot(x, Cerfc,'--', 'LineWidth', 2, 'DisplayName', 'Analytical');
xlabel('Depth (m)');
ylabel('Concentration (cm^-3)');
title('Constant Source Diffusion Profile at 1200 K');
hold off;
legend('show', 'Location', 'northeast')
fprintf("Constant Source Diffusion Time Taken = %.3e s\n", (iterations * deltaT)/(3600*24));

%% Answering Question 3

figure; 
hold on; 
grid on;

x=(0:depth_size-1)*deltaX;

iterations=800;
depth_size=100;
C_vector = zeros(iterations,depth_size);
C_vector(1,1:2) = [2e19 2e19];

C_vectorCSD=C_vector;
for i = 2:iterations
    for j = 2:depth_size %size of array -1

        if j==depth_size
            C_vectorCSD(i,j)=a*(C_vectorCSD(i-1,j-1));
        else
            C_vectorCSD(i,j)=a*(C_vectorCSD(i-1,j-1)+C_vectorCSD(i-1,j+1));
        end
    end
    C_vectorCSD(i,1:2)=[2e19 2e19];


end
plot(x, C_vectorCSD(end, :), 'LineWidth', 3, 'DisplayName', 'Equation 7.40');


C_vector = zeros(iterations,depth_size);
C_vector(1,1:2) = [2e19 2e19];

C_vectorCSD=C_vector;
for i = 2:iterations
    for j = 2:depth_size 

        if j==depth_size
            C_vectorCSD(i,j)=C_vectorCSD(i,j)+a*(C_vectorCSD(i-1,j-1)-2*(C_vectorCSD(i,j)));
        else
            C_vectorCSD(i,j)=C_vectorCSD(i,j)+a*(C_vectorCSD(i-1,j-1)+C_vectorCSD(i-1,j+1)-2*(C_vectorCSD(i,j)));
        end
    end
    C_vectorCSD(i,1:2)=[2e19 2e19];


end

plot(x, C_vectorCSD(end, :),'--', 'LineWidth', 2, 'DisplayName', 'Equation 7.38');
xlabel('Depth (m)');
ylabel('Concentration (cm^-3)');
title('Constant Source Diffusion Profile at 1200 K');
hold off;
legend('show', 'Location', 'northeast')
fprintf("Constant Source Diffusion Time Taken = %.3e s\n", (iterations * deltaT)/(3600*24));
