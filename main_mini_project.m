
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
C_vectorCSD=C_vector;
for i = 2:iterations
    for j = 2:depth_size %size of array -1
       
        if j==depth_size
            C_vectorCSD(i,j)=aCSD*(C_vectorCSD(i-1,j)+C_vectorCSD(i-1,j-1));
        else
            C_vectorCSD(i,j)=aCSD*(C_vectorCSD(i-1,j-1)+C_vectorCSD(i-1,j+1));
        end
    end
    C_vectorCSD(i,1:2) = [2e19 2e19];


end





figure;
plot(C_vectorCSD(end, :));
hold on;

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
plot(C_vectorCSD(end, :), '--');
xlabel('Depth index');
ylabel('Concentration');
title('Constant Source Diffusion Profile: Diffusion into infinite depth vs diffusion into finite depth');
legend('Finite depth', 'Infinite depth', 'Location', 'northeast');
grid on;

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
plot(C_vectorDID(end, :));
xlabel('Depth index');
ylabel('Concentration');
title('Final Diffusion Profile');
grid on;