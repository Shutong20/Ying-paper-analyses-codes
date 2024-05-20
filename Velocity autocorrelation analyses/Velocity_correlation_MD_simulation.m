function [VC_norm, VC_sem_norm] = Velocity_correlation_MD_simulation(GEM_traj)

% This function is modified from Velocity_correlation script located at 
% "/Users/tongshu/Documents/Lab project 2020/Mitotisc spindle Gem diffusion program/Tong added/Particle_tracking_vs_Intensity

% This function is to calculate the velocity autocorrelation which can be
% used to distinguish among distinct mechanisms for anomaloous
% sub-diffusion: Cv(tau) = <v(t+tau)*v(t)>; v(t) = (r(t+dt)-r(t))/dt

% Based on the reference: "Vast heterogeneity in cytoplasmic diffusion 
% rates revealed by nanorheology and Doppelgänger simulations"
% & "Analytical Tools To Distinguish the Effects of Localization Error, 
% Confinement, and Medium Elasticity on the Velocity Autocorrelation Function"

Max_traj_length = length(GEM_traj{1});
N_select_traj = length(GEM_traj);
dtau = 1; % This value represents different tau values when calculating velocity autocorrelation, default is 1
VC_time_traj = zeros(N_select_traj, Max_traj_length); % Velocity autocorrelation of each trajectory, autocorrelation array also includes tau=0

index = 0;
for j = 1:N_select_traj % Loop through different trajectories within the file
    index = index+1;
    VC_time_traj_temp = zeros(dtau, floor(Max_traj_length/dtau)-1); % Vc_time_traj_temp is a matrix, each row is the autocorrelation of one velocity array. As an example, row 1 is the velocty autocorrelation array of [x(1+dtau)-x(1),x(1+2*dtau)-x(1+dtau),x(1+3*dtau)-x(1+2*dtau),...]
    traj_x = GEM_traj{j}(:,1);
    traj_y = GEM_traj{j}(:,2);
    traj_z = GEM_traj{j}(:,3);    
    for array_num = 1:dtau % Loop through "dtau" velocity arrays that have intervals of dtau timelag
        traj_temp_x = traj_x(array_num:dtau:end);
        traj_temp_y = traj_y(array_num:dtau:end);
        traj_temp_z = traj_z(array_num:dtau:end);
        vx_traj = diff(traj_temp_x, 1)/dtau; % calculate x direction velocity using subsequent x displacements
        vy_traj = diff(traj_temp_y, 1)/dtau; % calculate y direction velocity using subsequent y displacements
        vz_traj = diff(traj_temp_z, 1)/dtau; % calculate z direction velocity using subsequent y displacements
        if length(vx_traj) > 1 % Only calculate the autocorrelation is there are at least two numbers within the velocity array
            VC_time_traj_temp(array_num, 1:length(vx_traj)) = (autocorr(vx_traj,'NumLags',length(vx_traj)-1) + autocorr(vy_traj,'NumLags',length(vy_traj)-1) + autocorr(vz_traj,'NumLags',length(vz_traj)-1))/3;
%             VC_time_traj_temp(array_num, 1:length(vx_traj)) = ...
%             (autocorr(vx_traj,'NumLags',length(vx_traj)-1) + autocorr(vy_traj,'NumLags',length(vy_traj)-1))/2; % Only calculate velocity autocorrelation on the 2D xy plane
        end
    end
    VC_time_traj_temp(VC_time_traj_temp == 0) = NaN; % Avoid averaging 0 values typically located at the end of the array
    VC_time_traj(index, 1:length(VC_time_traj_temp)) = mean(VC_time_traj_temp,1,'omitnan');
end

% Calculate trajectory-averaged velocity correlation based on individual 
% autocorrelations among all trajectories at different time lag
VC = zeros(1,Max_traj_length);
VC_sem = zeros(1,Max_traj_length);
for i = 1:Max_traj_length
    temp = VC_time_traj(:,i);
    VC(i) = nanmean(temp(temp~=0)); % There exist NAN values in VC_time_traj due to same detection positions over different timepoints
    VC_sem(i) = nanstd(temp(temp~=0))/sqrt(length(temp(temp~=0)));
end
VC_norm = VC/VC(1);
VC_sem_norm = VC_sem/VC(1);

close all
figure
hold on
% errorbar((0:length(VC_norm)-1)*dtau,VC_norm,VC_sem_norm,'o-','Color','[0.93,0.69,0.13]')
errorbar((0:length(VC_norm)-1)*dtau,VC_norm,VC_sem_norm,'o-','Color',c_jet_new(dtau))
xlabel('Time lag / simulation time unit')
ylabel('Normalization velocity autocorrelation')
box on
set(gca,'FontSize',15)

end