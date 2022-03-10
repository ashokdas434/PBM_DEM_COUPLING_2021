%% Part-2 of PBM code: PBM code of DEM-PBM coupling part
% Created by: Ashok Das (With help of T. De and G. Kaur)
% Date: 09 Aug 2020
%clear all; close all

function PBM_Main(Loop_counter)
tic
%Loop_counter = 1; % Initialization of loop counter
psi         = 0.001; % Prob. of successful aggregation
vol_range   = 10:10:200; % Volume option of particles in mm^3
t_total     = 10; % Total process time (sec)
PBM_trigger = 1.1; % Trigger criteria of PBM (1.05, 1.1, 1.5, 1.2, etc.)

%% Loading of beta_star
load('beta_star.mat');

load_DEM_input=Loop_counter-1;

%% Loading of total collision freq. and No. of particles from DEM simulations
coll_mat = load(['DEM_collision_mat_loop_',num2str(load_DEM_input),'.txt']); % Loading collision matrix % Change the name  to DEM_coll_mat.txt
tot_freq = 0.5* sum(coll_mat(:));

N_p      = load(['DEM_ntotal_loop_',num2str(load_DEM_input),'.txt']); % Loading tot no of particles % Change the name to DEM_ntotal.txt accordingly.

f_c      = 2*tot_freq/N_p;  % average freq. of collisions per particle

%% Loading process time and PSD from last PBM cycle
if Loop_counter == 1
    t_process = 0; % Initial time
    PSD = zeros(length(vol_range),1); % Initialization of PSD for time 0
    PSD(1) = N_p; % Initially all are monosized particles
else
    load(['PBE_results-',num2str(Loop_counter-1),'.mat']);
    PSD=PSD';
end
    Vol_tot= vol_range*PSD
%% PBE calculations & Saving data
t_span = [0 (t_total -t_process)];  % Temporary time variable
avg_vol = (vol_range*PSD) / sum(PSD);   % avg volume of system initially

options = odeset('Events',@event_function);

[t, N ,t_PBM,PSD] = ode45(@discrete_PBM,t_span,PSD,options,...
                    f_c,beta_star,psi,vol_range,avg_vol,PBM_trigger,t_span(2));  % t_PBM -> trigger time; V -> stopping volume

sprintf('At t = %1.2f seconds the trigger condition is satisfied.', t_PBM)

t_process = t_process + t_PBM; % Updating process time
sim_time=toc
f_call=f_c
N_p = sum(PSD)
Vol_tot= vol_range*PSD'
save(['PBE_results-',num2str(Loop_counter),'.mat'],'t_process', 'PSD','sim_time','f_call') % Saving data

%% PBM to DEM
rad_range = ((3*vol_range/(4*pi)).^(1/3))/1000; % radius range in meter

PBMtoDEM(PSD, length(vol_range), rad_range) % Creation of DEM code

%% set flag to stop loop
if (t_process >= t_total)
	fileID = fopen('text.txt','w');
	fprintf(fileID,'1');
	fclose(fileID);
	type text.txt
else
	fileID = fopen('text.txt','w');
	fprintf(fileID,'0');
	fclose(fileID);
	type text.txt
end

return
