close all;clc;clearvars -except Gaitcycle Continuous %removes all variables except for gaitcyle and continuous from the matlab workspace. Useful because loading these variables can take several minutes.
%load the InclineExperiment data from the folder location you specify
% load('Z:\your_file_location_here\InclineExperiment.mat') 
%This file is large and may take several minutes to load. We recommend not
%clearing the variable, and only loading it as many times as necessary.

%% Example: removing strides with outliers from kinematic data, and taking the mean

% Gaitcyle and Continuous store the InclineExperiment data in a MATLAB
% structure array, or 'struct'

% The fields of a struct can be iterated over using a cell array

% A cell array can be defined manually with this notation: 
% sub={'AB01','AB02','AB03','AB04','AB05','AB06','AB07','AB08','AB09','AB10'};
sub={'AB01'};

% The command 'fieldnames' also returns the name of all fields within a
% struct as a cell array

trials = {'s0x8i0','s1i0','s1x2i0'}
leg={'left'};
joint={'foot'};
percent_gait = linspace(0,1,150);

c = {'r','g','b'}

% Initializing a matrix that will contain the mean joint data of all subjects,
%all tasks
% concatenated_data=NaN(150,numel(sub)*numel(trials));

figure(1)
subplot(3,1,1)
hold on
subplot(3,1,2)
hold on
subplot(3,1,3)
hold on

figure(2)
hold on
footAngleZero = 85.46;


figure(3)
subplot(3,1,1)
hold on
subplot(3,1,2)
hold on
subplot(3,1,3)
hold on

figure(4)
hold on


for i = 1:length(sub) %loop through all subjects in 'sub'
    for j = 1:length(trials) %loop through all trials in 'trial'
        % store the kinematic data in a temporary variable
%         sub{i}
        trial = trials{j};
        foot_angle_data = Gaitcycle.(sub{i}).(trials{j}).kinematics.jointangles.(leg{1}).('foot').x;
        % adjust
        foot_angle_data = -(foot_angle_data + footAngleZero);
        
        ankle_angle_data = Gaitcycle.(sub{i}).(trials{j}).kinematics.jointangles.(leg{1}).('ankle').x;
        
        
        time_data = Gaitcycle.(sub{i}).(trials{j}).cycles.(leg{1}).time;

        % delete the strides identified to contain outliers
        foot_angle_data(:,Gaitcycle.(sub{i}).(trials{j}).stepsout.(leg{1})) = [];
        time_data(:,Gaitcycle.(sub{i}).(trials{j}).stepsout.(leg{1})) = [];
        ankle_angle_data(:,Gaitcycle.(sub{i}).(trials{j}).stepsout.(leg{1})) = [];
        
        
        shank_angle_data = foot_angle_data - ankle_angle_data;
        
        switch trial
            
            case 's0x8i0'
                treadmillSpeed = 0.8;
                
            case 's1i0'
                treadmillSpeed = 1;
                
            case 's1x2i0'
                treadmillSpeed = 1.2;
                
        end
        
        
        
        [rows,cols] = size(foot_angle_data);
        
        for k = 1:cols
            
            foot_angle_data_col = foot_angle_data(:,k);
            shank_angle_data_col = shank_angle_data(:,k);
            
            time_data_col = time_data(:,k);
            
            stepDuration = time_data_col(end) - time_data_col(1);
            phaseDot = 1/stepDuration;
            stepLength = treadmillSpeed*stepDuration;
            figure(1)
%             plot_i = plot3(percent_gait,stepLength*ones(size(percent_gait)), temp_data_col, c{j},'LineWidth',2);
            subplot(3,1,1)
            plot_i = plot(percent_gait, foot_angle_data_col, c{j},'LineWidth',2);
            plot_i.Color(4) = 0.2;

            subplot(3,1,2)
            plot_i = plot(stepLength*ones(size(percent_gait)), foot_angle_data_col, c{j},'LineWidth',2);
            plot_i.Color(4) = 0.2;
            
            
            subplot(3,1,3)
            plot_j = plot(phaseDot*ones(size(percent_gait)), foot_angle_data_col, c{j},'LineWidth',2);
            plot_j.Color(4) = 0.2;
            
            figure(2)
            plot_j = plot3(percent_gait,stepLength*ones(size(percent_gait)), foot_angle_data_col, c{j},'LineWidth',2);
            
            
            figure(3)
%             plot_i = plot3(percent_gait,stepLength*ones(size(percent_gait)), temp_data_col, c{j},'LineWidth',2);
            subplot(3,1,1)
            plot_i = plot(percent_gait, shank_angle_data_col, c{j},'LineWidth',2);
            plot_i.Color(4) = 0.2;

            subplot(3,1,2)
            plot_i = plot(stepLength*ones(size(percent_gait)), shank_angle_data_col, c{j},'LineWidth',2);
            plot_i.Color(4) = 0.2;
            
            
            subplot(3,1,3)
            plot_j = plot(phaseDot*ones(size(percent_gait)), shank_angle_data_col, c{j},'LineWidth',2);
            plot_j.Color(4) = 0.2;
            
            figure(4)
            plot_j = plot3(percent_gait,stepLength*ones(size(percent_gait)), shank_angle_data_col, c{j},'LineWidth',2);



            
        end
        
        
        
%         pause
        
        % take the mean of the remaining strides and concatenate them
        
%         nanmean(temp_data,2)
%         concatenated_data(:, numel(trials)*(s-1)+t ) = nanmean(temp_data,2);
    end
end

% The data in Gaitcyle always has 150 points per stride, and they are
% evenly spaced with respect to percent gait.

figure(1)
title([leg{1}, ' ', joint{1}, ' joint angles'])
subplot(3,1,1)

xlabel('percent gait')
% ylabel('Step Length (m)')
ylabel('degrees'); 

subplot(3,1,2)

% ylabel('Step Length (m)')
xlabel('step Length'); 
ylabel('degrees'); 

subplot(3,1,3)

% ylabel('Step Length (m)')
xlabel('phase_dot'); 
ylabel('degrees'); 


figure(2)
title([leg{1}, ' ', joint{1}, ' joint angles'])
xlabel('percent gait')
ylabel('Step Length (1/s)')
zlabel('degrees'); 


figure(3)
title([leg{1}, ' ', 'shank', ' joint angles'])
subplot(3,1,1)

xlabel('percent gait')
% ylabel('Step Length (m)')
ylabel('degrees'); 

subplot(3,1,2)

% ylabel('Step Length (m)')
xlabel('step Length'); 
ylabel('degrees'); 

subplot(3,1,3)

% ylabel('Step Length (m)')
xlabel('phase_dot'); 
ylabel('degrees'); 


figure(4)
title([leg{1}, ' ', 'shank', ' joint angles'])
xlabel('percent gait')
ylabel('Step Length (1/s)')
zlabel('degrees'); 


