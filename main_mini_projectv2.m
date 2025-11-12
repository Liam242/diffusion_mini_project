
% Clear everything
clear; clc; close all;

depth_size=100;




%% Constant Source Diffusion
% Simplify using the assumtion 
% DΔt/(Δx^2)=1/2
aCSD=0.5;


figure;
hold on;

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
                C_vectorCSD(i, j) = aCSD * (C_vectorCSD(i-1, j-1));
            else
                C_vectorCSD(i, j) = aCSD * (C_vectorCSD(i-1, j-1) + C_vectorCSD(i-1, j+1));
            end
        end
        % Apply boundary condition at the start of each iteration
        C_vectorCSD(i, 1:2) = [2e19 2e19];
    end

    % Plot the last iteration for each case
    plot(C_vectorCSD(end, :), 'DisplayName', sprintf('%d iterations', iterations));
end

xlabel('Depth Index');
ylabel('Concentration (C)');
title('Concentration Profile for Different Iteration Counts');
legend show;
grid on;
hold off;

% %% Drive-In Diffusion
% % Simplify using the assumtion 
% % DΔt/(Δx^2)=1/2
% aDID=0.5;
% iterations=1000;
% C_vectorDID=zeros(iterations,depth_size);
% C_vectorDID(1, :)=C_vectorCSD(end, :);
% 
% for i = 2:iterations
%     print("iterations")
%     for j = 1:depth_size %size of array -1
% 
% 
% 
% 
% 
%         if j==depth_size
%             C_vectorDID(i,j)=aDID*(C_vectorDID(i-1,j)+C_vectorDID(i-1,j-1));
%         elseif j==1
%             C_vectorDID(i,j)=aDID*(C_vectorDID(i-1,j)+C_vectorDID(i-1,j+1));
%         else
%             C_vectorDID(i,j)=aDID*(C_vectorDID(i-1,j-1)+C_vectorDID(i-1,j+1));
%         end
%     end
% 
% 
% 
% end
% % Plot: Constant source diffusion
% figure;
% plot(C_vectorDID(end, :));
% xlabel('Depth index');
% ylabel('Concentration');
% title('Final Diffusion Profile');
% grid on;