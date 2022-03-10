%% Part-1 of PBM code: Extraction of beta* at the initial stage
% Created by: Ashok Das (With help of T. De and G. Kaur)
% Date: 09 Aug 2020
clear all;  close all

vol_range = 1:20;

%% Loading collision data for uniform distribution
temp = load('0.txt'); % load data change uniform name
Mat = temp(vol_range,vol_range); %collision matrix
S_Mat = sum(Mat(:));

%% Calculation, Data saving, and Plotting
beta_star = Mat/S_Mat; % Calculation of beta_star

save('beta_star.mat', 'beta_star')  % Saving data

figure
surf(vol_range,vol_range,beta_star);
xlabel('x','fontsize',18); ylabel('y','fontsize',18); zlabel('\beta^*','fontsize',18);

