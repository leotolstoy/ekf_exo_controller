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
% sub={'AB01','AB05','AB08','AB09'};

sub={'AB01'};

% The command 'fieldnames' also returns the name of all fields within a
% struct as a cell array

% trials = {'s1i10','s1i7x5','s1i5','s1i2x5','s1i0','s1d2x5','s1d5','s1d7x5','s1d10'}
trials = fieldnames(Gaitcycle.AB01);
% The command 'setdiff' is useful for removing cells of a given name
trials = setdiff(trials,'subjectdetails');


leg={'left'};
joint={'foot'};
percent_gait = linspace(0,1,150);

%filter set up

%High Pass filter
s=tf('s');
omega = 0.01 * 2*pi;
zeta= 0.7;
omega_1 = 0.01*2*pi;
zeta_1 = 0.2;
HIGH_PASS_FILTER = (s^2)/(s^2+2*zeta*omega*s+omega^2)*(s^2/(s^2+2*zeta_1*omega_1*s+omega_1^2));


[num,den] = tfdata(HIGH_PASS_FILTER);
[A_HPF,B_HPF,C_HPF,D_HPF] = tf2ss(num{1},den{1});
dT = 1/100;

Ad_HPF = expm(A_HPF*dT);
Bd_HPF = (Ad_HPF - eye(4)) * A_HPF^-1 *B_HPF;
HP_filter_states_x0_heel = [-0.300;-0.300;0;0];


%Low Pass filter
LOW_PASS_FILTER = (1)/(s^2+2*zeta*omega*s+omega^2)*(s^2/(s^2+2*zeta_1*omega_1*s+omega_1^2));
[num,den] = tfdata(LOW_PASS_FILTER)
[A_LPF,B_LPF,C_LPF,D_LPF] = tf2ss(num{1},den{1});

Ad_LPF = expm(A_LPF*dT);
Bd_LPF = (Ad_LPF - eye(4)) * A_LPF^-1 *B_LPF;



LP_filter_states_x0_heel = [-1.50;0.400;0;0];


A_mat_master = [];
b_footAngle_master = [];
b_shankAngle_master = [];

b_heelAccForward_master = [];
b_tibiaAccForward_master = [];


for i = 1:length(sub) %loop through all subjects in 'sub'
    
    switch sub{i}
            
        case 'AB01'
            footAngleZero = 85.46; %subject 1

        case 'AB02'
            footAngleZero = 85.46 + 4.506 + 1.62; %subject 2

        case 'AB03'
            footAngleZero = 85.46 + 5.564; %subject 3
        case 'AB04'
            footAngleZero = 85.46 + 7.681; %subject 4
        case 'AB05'
            footAngleZero = 85.46 + 5.405; %subject 5
        case 'AB06'
            footAngleZero = 85.46 + 4.089; %subject 6
        case 'AB07'
            footAngleZero = 85.46 + 1.523; %subject 7
        case 'AB08'
            footAngleZero = 85.46 + 3.305; %subject 8
        case 'AB09'
            footAngleZero = 85.46 + 4.396; %subject 9
        case 'AB10'
            footAngleZero = 85.46 + 6.555; %subject 10

    end


    
    subjLegLength = Gaitcycle.(sub{i}).subjectdetails{5,2}/1000
%     pause
    for j = 1:length(trials) %loop through all trials in 'trial'
        % store the kinematic data in a temporary variable
%         sub{i}
        trial = trials{j};
        [treadmillSpeed, treadmillIncline] = returnSpeedIncline(trial);
        
        foot_angle_data = Gaitcycle.(sub{i}).(trials{j}).kinematics.jointangles.(leg{1}).('foot').x;
        % adjust
        foot_angle_data = -(foot_angle_data + footAngleZero);
        
        ankle_angle_data = Gaitcycle.(sub{i}).(trials{j}).kinematics.jointangles.(leg{1}).('ankle').x;
        
        % get Heel Marker
        heel_pos_x_data = Gaitcycle.(sub{i}).(trials{j}).kinematics.markers.(leg{1}).('heel').x;
        heel_pos_y_data = Gaitcycle.(sub{i}).(trials{j}).kinematics.markers.(leg{1}).('heel').y/1000;
        heel_pos_z_data = Gaitcycle.(sub{i}).(trials{j}).kinematics.markers.(leg{1}).('heel').z;
        
        
        
        
        %adjust for treadmill frame
        
        if treadmillIncline < 0
            heel_pos_y_data = -heel_pos_y_data;
        end
        
        
        time_data = Gaitcycle.(sub{i}).(trials{j}).cycles.(leg{1}).time;
        

        % delete the strides identified to contain outliers
        foot_angle_data(:,Gaitcycle.(sub{i}).(trials{j}).stepsout.(leg{1})) = [];
        time_data(:,Gaitcycle.(sub{i}).(trials{j}).stepsout.(leg{1})) = [];
        ankle_angle_data(:,Gaitcycle.(sub{i}).(trials{j}).stepsout.(leg{1})) = [];
        heel_pos_x_data(:,Gaitcycle.(sub{i}).(trials{j}).stepsout.(leg{1})) = [];
        heel_pos_y_data(:,Gaitcycle.(sub{i}).(trials{j}).stepsout.(leg{1})) = [];
        heel_pos_z_data(:,Gaitcycle.(sub{i}).(trials{j}).stepsout.(leg{1})) = [];
        
        
    
        
        shank_angle_data = foot_angle_data - ankle_angle_data;
        
        
        incline_scaled = treadmillIncline;

            
        
        [~,cols] = size(foot_angle_data);
        

        for k = 1:cols
            
            foot_angle_data_col = foot_angle_data(:,k);
            shank_angle_data_col = shank_angle_data(:,k);
            heel_pos_y_data_col = heel_pos_y_data(:,k);
            
            %numerically differentiate the heel_pos data
            time_data_col = time_data(:,k);
            heel_vel_y_data_col = diff(heel_pos_y_data_col)./diff(time_data_col);
            
            heel_acc_y_data_col = diff(heel_vel_y_data_col)./diff(time_data_col(2:end));
            
            heel_acc_y_data_col = [heel_acc_y_data_col(1);heel_acc_y_data_col(1);heel_acc_y_data_col];
            
            
            %Highpass filter heel position data
%             [rows,cols1] = size(heel_pos_y_data_col);
            heel_pos_y_data_col_filt = zeros(size(heel_pos_y_data_col));
            
            HPF_states_heel = zeros(4,150);
            
            for ii = 1:150
                if ii == 1
                    HPF_states_heel(:,ii) = (Ad_HPF) * HP_filter_states_x0_heel + Bd_HPF * heel_pos_y_data_col(ii);
                    heel_pos_y_data_col_filt(ii) = heel_pos_y_data_col(ii);
                else
                    HPF_states_heel(:,ii) = (Ad_HPF) * HPF_states_heel(:,ii-1) + Bd_HPF * heel_pos_y_data_col(ii-1);
                    heel_pos_y_data_col_filt(ii) = C_HPF * HPF_states_heel(:,ii-1) + D_HPF * heel_pos_y_data_col(ii-1);

                end

            end
            
            heel_pos_y_data_col_filt(1) = heel_pos_y_data_col_filt(2);
            HP_filter_states_x0_heel = HPF_states_heel(:,end);
            
            %Lowpass filter heel acc data
            heel_acc_y_data_col_filt = zeros(size(heel_acc_y_data_col));
            LPF_states_heel = zeros(4,150);
            
            for ii = 1:150
                if ii == 1
                    LPF_states_heel(:,ii) = (Ad_LPF) * LP_filter_states_x0_heel + Bd_LPF * heel_acc_y_data_col(ii);
                    heel_acc_y_data_col_filt(ii) = heel_acc_y_data_col(ii);
                else
                    LPF_states_heel(:,ii) = (Ad_LPF) * LPF_states_heel(:,ii-1) + Bd_LPF * heel_acc_y_data_col(ii-1);
                    heel_acc_y_data_col_filt(ii) = C_LPF * LPF_states_heel(:,ii-1) + D_LPF * heel_acc_y_data_col(ii-1);

                end

            end
            
            heel_acc_y_data_col_filt(1) = heel_acc_y_data_col_filt(2);
            LP_filter_states_x0_heel = LPF_states_heel(:,end);
            
            
            
            
            
            figure(100)
            subplot(2,1,1)
            plot(time_data_col, heel_pos_y_data_col_filt,'r')
            hold on
            plot(time_data_col, heel_pos_y_data_col,'b')
%             plot(time_data_col, heel_acc_y_data_col,'g')
            plot(time_data_col, heel_acc_y_data_col_filt,'g')
            
            subplot(2,1,2)
            plot(time_data_col, HPF_states_heel,'r')
            hold on
            plot(time_data_col, LPF_states_heel,'g')
            
            pause

            
            
            stepDuration = time_data_col(end) - time_data_col(1);
            phaseDot = 1/stepDuration;
            stepLength = treadmillSpeed*stepDuration;
            
            
            stepLength_scaled = stepLength/subjLegLength;
            phase_col = time_data_col/stepDuration;

            
        end

    end
end

% The data in Gaitcyle always has 150 points per stride, and they are
% evenly spaced with respect to percent gait.



%% Helper function
function [speed, incline] = returnSpeedIncline(trialString)



d_idx = strfind(trialString,'d');
i_idx = strfind(trialString,'i');

isIncline = ~isempty(i_idx) & isempty(d_idx);
isDecline = isempty(i_idx) & ~isempty(d_idx);

mid_idx = [i_idx,d_idx];

speedString = trialString(2:mid_idx-1);

inclineString = trialString(mid_idx+1:end);



switch speedString
            
    case '0x8'
        speed = 0.8;

    case '1'
        speed = 1;

    case '1x2'
        speed = 1.2;

end


switch inclineString
            
    case '10'
        incline = 10;

    case '7x5'
        incline = 7.5;

    case '5'
        incline = 5;
        
    case '2x5'
        incline = 2.5;
        
    case '0'
        incline = 0;
        

end
if isDecline
    incline = incline * -1;
end



end
