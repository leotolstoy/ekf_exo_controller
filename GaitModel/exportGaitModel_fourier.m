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
sub={'AB01','AB02','AB03','AB04','AB05','AB06','AB07','AB08','AB09','AB10'};
% sub={'AB01','AB05','AB08','AB09'};

% sub={'AB01'};

% The command 'fieldnames' also returns the name of all fields within a
% struct as a cell array

% trials = {'s1i10','s1i7x5','s1i5','s1i2x5','s1i0','s1d2x5','s1d5','s1d7x5','s1d10'}
trials = fieldnames(Gaitcycle.AB01);
% The command 'setdiff' is useful for removing cells of a given name
trials = setdiff(trials,'subjectdetails');


leg={'left'};
joint={'foot'};
percent_gait = linspace(0,1,150);

c = {'r','g','b'}
% c1 = {'b','b','b','b','b','r','r','r','r'}
% 
% c2 = [0.5,0.4,0.3,0.2,0.1,0.2,0.3,0.4,0.5]

% Initializing a matrix that will contain the mean joint data of all subjects,
%all tasks
% concatenated_data=NaN(150,numel(sub)*numel(trials));

% figure(1)
% subplot(1,3,1)
% hold on
% subplot(1,3,2)
% hold on
% subplot(1,3,3)
% hold on

figure(2)
subplot(1,3,1)
hold on
subplot(1,3,2)
hold on
subplot(1,3,3)
hold on





% figure(3)
% subplot(1,3,1)
% hold on
% subplot(1,3,2)
% hold on
% subplot(1,3,3)
% hold on

figure(4)
subplot(1,3,1)
hold on
subplot(1,3,2)
hold on
subplot(1,3,3)
hold on



figure(50)
subplot(1,3,1)
hold on
subplot(1,3,2)
hold on
subplot(1,3,3)
hold on

% figure(51)
% subplot(1,3,1)
% hold on
% subplot(1,3,2)
% hold on
% subplot(1,3,3)
% hold on

figure(52)
subplot(1,3,1)
hold on
subplot(1,3,2)
hold on
subplot(1,3,3)
hold on

phaseDelins = [0.1,0.5,0.65,1]

numInclineFuncs = 2;
numStepLengthFuncs = 2;
N_FOURIER = 20
numPhaseFuncs = (length(1:1:N_FOURIER) * 2) + 1;
numFuncs = numInclineFuncs*numStepLengthFuncs*numPhaseFuncs;


%High Pass filter
s=tf('s');
omega = 0.5 * 2*pi;
zeta= 0.9;
% omega_1 = 0.4*2*pi;
% zeta_1 = 0.2;
% HIGH_PASS_FILTER = (s^2)/(s^2+2*zeta*omega*s+omega^2)*(s^2/(s^2+2*zeta_1*omega_1*s+omega_1^2));

HIGH_PASS_FILTER = (s^2)/(s^2+2*zeta*omega*s+omega^2);

[num,den] = tfdata(HIGH_PASS_FILTER);
[A_HPF,B_HPF,C_HPF,D_HPF] = tf2ss(num{1},den{1});
dT = 1/100;

Ad_HPF = expm(A_HPF*dT);
Bd_HPF = (Ad_HPF - eye(2)) * A_HPF^-1 *B_HPF;
HP_filter_states_x0_heel_y = [-0.300;-0.300;];
HP_filter_states_x0_heel_z = [-0.300;-0.300;];
HP_filter_states_x0_tibia = [-0.300;-0.300;];


A_mat_master = [];
b_footAngle_master = [];
b_shankAngle_master = [];

b_heelPosForward_master = [];
b_heelPosUp_master = [];
b_tibiaPosForward_master = [];

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
        heel_pos_x_data = Gaitcycle.(sub{i}).(trials{j}).kinematics.markers.(leg{1}).('heel').x/1000;
        heel_pos_y_data = Gaitcycle.(sub{i}).(trials{j}).kinematics.markers.(leg{1}).('heel').y/1000;
        heel_pos_z_data = Gaitcycle.(sub{i}).(trials{j}).kinematics.markers.(leg{1}).('heel').z/1000;
        
        
        % get Tibia Marker
        tibia_pos_y_data = Gaitcycle.(sub{i}).(trials{j}).kinematics.markers.(leg{1}).('tibia').y/1000;
        
        
        %adjust for treadmill frame
        
        if treadmillIncline < 0
            heel_pos_y_data = -heel_pos_y_data;
            tibia_pos_y_data = -tibia_pos_y_data;
        end
        
        
        time_data = Gaitcycle.(sub{i}).(trials{j}).cycles.(leg{1}).time;
        

        % delete the strides identified to contain outliers
        foot_angle_data(:,Gaitcycle.(sub{i}).(trials{j}).stepsout.(leg{1})) = [];
        time_data(:,Gaitcycle.(sub{i}).(trials{j}).stepsout.(leg{1})) = [];
        ankle_angle_data(:,Gaitcycle.(sub{i}).(trials{j}).stepsout.(leg{1})) = [];
        heel_pos_x_data(:,Gaitcycle.(sub{i}).(trials{j}).stepsout.(leg{1})) = [];
        heel_pos_y_data(:,Gaitcycle.(sub{i}).(trials{j}).stepsout.(leg{1})) = [];
        heel_pos_z_data(:,Gaitcycle.(sub{i}).(trials{j}).stepsout.(leg{1})) = [];
        
        
        
        
        
        
%         tibia_pos_y_data(:,Gaitcycle.(sub{i}).(trials{j}).stepsout.(leg{1})) = [];
        
         %get forward heel marker linear velocity
%         heel_vel_y_data = diff(heel_pos_y_data)./diff(time_data);
%         heel_vel_y_data = [heel_vel_y_data;heel_vel_y_data(end,:)];
%         heel_vel_y_data = lowpass(heel_vel_y_data, 90,100);
%         heel_vel_y_data = highpass(heel_vel_y_data, 0.1,100);
% 
%         %get forward heel marker linear Pos
%         heel_Pos_y_data = diff(heel_vel_y_data)./diff(time_data);
%         heel_Pos_y_data = [heel_Pos_y_data;heel_Pos_y_data(end,:)];
%         heel_Pos_y_data = lowpass(heel_Pos_y_data, 90,100);
%         heel_Pos_y_data = highpass(heel_Pos_y_data, 0.1,100);
        
        
%         %get forward tibia marker linear velocity
%         tibia_vel_y_data = diff(tibia_pos_y_data)./diff(time_data);
%         tibia_vel_y_data = [tibia_vel_y_data;tibia_vel_y_data(end,:)];
%         tibia_vel_y_data = lowpass(tibia_vel_y_data, 90,100);
%         tibia_vel_y_data = highpass(tibia_vel_y_data, 0.1,100);
% 
%         %get forward tibia marker linear Pos
%         tibia_Pos_y_data = diff(tibia_vel_y_data)./diff(time_data);
%         tibia_Pos_y_data = [tibia_Pos_y_data;tibia_Pos_y_data(end,:)];
%         tibia_Pos_y_data = lowpass(tibia_Pos_y_data, 90,100);
%         tibia_Pos_y_data = highpass(tibia_Pos_y_data, 0.1,100);

        shank_angle_data = foot_angle_data - ankle_angle_data;        
        incline_scaled = treadmillIncline;
        
        [~,cols] = size(foot_angle_data);
        
        A_mat = zeros(cols * 150,numFuncs);
        b_footAngle_mat = zeros(cols * 150,1);
        b_shankAngle_mat = zeros(cols * 150,1);
        b_heelPosForward_mat = zeros(cols * 150,1);
        b_heelPosUp_mat = zeros(cols * 150,1);
        b_tibiaPosForward_mat = zeros(cols * 150,1);

        avgHeelPosForward = zeros(150,1);
        avgHeelPosUp = zeros(150,1);
        avgTibiaPosForward = zeros(150,1);
        for k = 1:cols
            
            foot_angle_data_col = foot_angle_data(:,k);
            shank_angle_data_col = shank_angle_data(:,k);
 
            heel_pos_y_data_col = heel_pos_y_data(:,k);
            heel_pos_z_data_col = heel_pos_z_data(:,k);
            tibia_pos_y_data_col = tibia_pos_y_data(:,k);
            
            %Highpass filter heel position data
%             [rows,cols1] = size(heel_pos_y_data_col);
            heel_pos_y_data_col_filt = zeros(size(heel_pos_y_data_col));
            heel_pos_z_data_col_filt = zeros(size(heel_pos_z_data_col));
            tibia_pos_y_data_col_filt = zeros(size(tibia_pos_y_data_col));
            states_heel_y = zeros(2,150);
            states_heel_z = zeros(2,150);
            states_tibia = zeros(2,150);
            for ii = 1:150

                
                if ii == 1
                    
                    states_heel_y(:,ii) = (Ad_HPF) * HP_filter_states_x0_heel_y + Bd_HPF * heel_pos_y_data_col(ii);
                    states_heel_z(:,ii) = (Ad_HPF) * HP_filter_states_x0_heel_z + Bd_HPF * heel_pos_z_data_col(ii);
                    states_tibia(:,ii) = (Ad_HPF) * HP_filter_states_x0_tibia + Bd_HPF * tibia_pos_y_data_col(ii);

                    heel_pos_y_data_col_filt(ii) = heel_pos_y_data_col(ii);
                    heel_pos_z_data_col_filt(ii) = heel_pos_z_data_col(ii);
                    tibia_pos_y_data_col_filt(ii) = tibia_pos_y_data_col(ii);
                else

                    states_heel_y(:,ii) = (Ad_HPF) * states_heel_y(:,ii-1) + Bd_HPF * heel_pos_y_data_col(ii-1);
                    states_heel_z(:,ii) = (Ad_HPF) * states_heel_z(:,ii-1) + Bd_HPF * heel_pos_z_data_col(ii-1);
                    states_tibia(:,ii) = (Ad_HPF) * states_tibia(:,ii-1) + Bd_HPF * tibia_pos_y_data_col(ii-1);
                    
                    heel_pos_y_data_col_filt(ii) = C_HPF * states_heel_y(:,ii-1) + D_HPF * heel_pos_y_data_col(ii-1);
                    heel_pos_z_data_col_filt(ii) = C_HPF * states_heel_z(:,ii-1) + D_HPF * heel_pos_z_data_col(ii-1);
                    tibia_pos_y_data_col_filt(ii) = C_HPF * states_tibia(:,ii-1) + D_HPF * tibia_pos_y_data_col(ii-1);


                end



            end
            heel_pos_y_data_col_filt(1) = heel_pos_y_data_col_filt(2);
            HP_filter_states_x0_heel_y = states_heel_y(:,end);

            heel_pos_z_data_col_filt(1) = heel_pos_z_data_col_filt(2);
            HP_filter_states_x0_heel_z = states_heel_z(:,end);
            
            tibia_pos_y_data_col_filt(1) = tibia_pos_y_data_col_filt(2);
            HP_filter_states_x0_tibia = states_tibia(:,end);
            
            
            
%             tibia_Pos_y_data_col = tibia_pos_y_data(:,k);
            
            time_data_col = time_data(:,k);
            
            figure(100)
        
            plot(time_data_col, heel_pos_y_data_col_filt,'r')
            hold on
            plot(time_data_col, heel_pos_y_data_col,'b')
%             plot(time_data_col, HP_filter_states_x0_heel_y)
            
            pause
            
            heel_pos_y_data_col = heel_pos_y_data_col_filt;
            avgHeelPosForward = avgHeelPosForward + heel_pos_y_data_col;

            heel_pos_z_data_col = heel_pos_z_data_col_filt;
            avgHeelPosUp = avgHeelPosUp + heel_pos_z_data_col;
            
            tibia_pos_y_data_col = tibia_pos_y_data_col_filt;
            avgTibiaPosForward = avgTibiaPosForward + tibia_pos_y_data_col;
            
            
            stepDuration = time_data_col(end) - time_data_col(1);
            phaseDot = 1/stepDuration;
            stepLength = treadmillSpeed*stepDuration;
            
            
            stepLength_scaled = stepLength/subjLegLength;
            phase_col = time_data_col/stepDuration;
            
%             figure(100)
%             
%             
%             if treadmillIncline == 10
%                 cl = 'r';
%             elseif treadmillIncline == 5
%                 cl = 'b';
%             elseif treadmillIncline == 0
%                 cl = 'g';
%             elseif treadmillIncline == -5
%                 cl = 'k';
%                
%             elseif treadmillIncline == -10
%                 cl = 'y';
%                 
%             end
%             subplot(1,3,1)
%             plot(time_data_col, heel_pos_x_data_col)
%             hold on
%             subplot(1,3,2)
%             plot(time_data_col, heel_Pos_y_data_col,cl)
%             hold on
%             plot(time_data_col, heel_vel_y_data_col*1e-1)
%             plot(time_data_col, heel_Pos_y_data_col*1e-2)
%             legend('pos','vel','Pos')
%             subplot(1,3,3)
%             plot(time_data_col, heel_pos_z_data_col)
%             hold on
%             
            
            
            
            
            for phase_idx = 1:length(phase_col)
                mat_idx = phase_idx + (k - 1)*150;
                phase_j = phase_col(phase_idx);

                fourier_coeffs = returnFourierBasis_Eval(phase_j,stepLength_scaled,incline_scaled,N_FOURIER);


                A_mat(mat_idx,:) = fourier_coeffs;
                b_footAngle_mat(mat_idx) = foot_angle_data_col(phase_idx);
                b_shankAngle_mat(mat_idx) = shank_angle_data_col(phase_idx);
                b_heelPosForward_mat(mat_idx) = heel_pos_y_data_col(phase_idx);
                b_heelPosUp_mat(mat_idx) = heel_pos_z_data_col(phase_idx);
                b_tibiaPosForward_mat(mat_idx) = tibia_pos_y_data_col(phase_idx);
             
                
        
            end
            
    
            
            
        end
%         pause
        A_mat_master = [A_mat_master;A_mat];
        b_footAngle_master = [b_footAngle_master;b_footAngle_mat];
        b_shankAngle_master = [b_shankAngle_master;b_shankAngle_mat];
        b_heelPosForward_master = [b_heelPosForward_master;b_heelPosForward_mat];
        b_heelPosUp_master = [b_heelPosUp_master;b_heelPosUp_mat];
        b_tibiaPosForward_master = [b_tibiaPosForward_master;b_tibiaPosForward_mat];
       
        avgStepDuration = mean(time_data(end,:) - time_data(1,:) );
        avgFootAngle = mean(foot_angle_data,2);
        avgShankAngle = mean(shank_angle_data,2);
        
        avgHeelPosForward = avgHeelPosForward /cols;
        avgHeelPosUp = avgHeelPosUp /cols;
        avgTibiaPosForward = avgTibiaPosForward /cols;
%         avgHeelPosForward = mean(heel_pos_y_data,2);
%         avgTibiaPosForward = mean(tibia_Pos_y_data,2);
        
        avgPhase = percent_gait';
        
        
        if strcmp(sub{i}, 'AB01') || true

            figure(2)
            if treadmillIncline == 10
                subplot(1,3,1)
            elseif treadmillIncline == 0
                subplot(1,3,2)
            elseif treadmillIncline == -10
                subplot(1,3,3)
            end
            
            if treadmillSpeed == 1.2
                cl = 'r';
            elseif treadmillSpeed == 1
                cl = 'b';
            elseif treadmillSpeed == 0.8
                cl = 'g';
            end
            
            
            if treadmillIncline == 10 || treadmillIncline == 0 || treadmillIncline == -10
                plot_j = plot3(avgPhase,stepLength_scaled*ones(size(avgPhase)), avgFootAngle, cl,'LineWidth',2);
                plot_j.Color(4) = 1;
            end


            figure(4)
            if treadmillIncline == 10
                subplot(1,3,1)
            elseif treadmillIncline == 0
                subplot(1,3,2)
            elseif treadmillIncline == -10
                subplot(1,3,3)
            end
            if treadmillSpeed == 1.2
                cl = 'r';
            elseif treadmillSpeed == 1
                cl = 'b';
            elseif treadmillSpeed == 0.8
                cl = 'g';
            end
            if treadmillIncline == 10 || treadmillIncline == 0 || treadmillIncline == -10
                plot_j = plot3(avgPhase,stepLength_scaled*ones(size(avgPhase)), avgShankAngle, cl,'LineWidth',2);
                plot_j.Color(4) = 1;
            end
            
            figure(50)
            if treadmillIncline == 10
                subplot(1,3,1)
            elseif treadmillIncline == 0
                subplot(1,3,2)
            elseif treadmillIncline == -10
                subplot(1,3,3)
            end
            if treadmillSpeed == 1.2
                cl = 'r';
            elseif treadmillSpeed == 1
                cl = 'b';
            elseif treadmillSpeed == 0.8
                cl = 'g';
            end
            if treadmillIncline == 10 || treadmillIncline == 0 || treadmillIncline == -10
                plot_j = plot3(avgPhase,stepLength_scaled*ones(size(avgPhase)), avgHeelPosForward, cl,'LineWidth',2);
                plot_j.Color(4) = 1;
            end
            
%             figure(51)
%             if treadmillIncline == 10
%                 subplot(1,3,1)
%             elseif treadmillIncline == 0
%                 subplot(1,3,2)
%             elseif treadmillIncline == -10
%                 subplot(1,3,3)
%             end
%             if treadmillSpeed == 1.2
%                 cl = 'r';
%             elseif treadmillSpeed == 1
%                 cl = 'b';
%             elseif treadmillSpeed == 0.8
%                 cl = 'g';
%             end
%             if treadmillIncline == 10 || treadmillIncline == 0 || treadmillIncline == -10
%                 plot_j = plot3(avgPhase,stepLength_scaled*ones(size(avgPhase)), avgTibiaPosForward, cl,'LineWidth',2);
%                 plot_j.Color(4) = 1;
%             end

            figure(52)
            if treadmillIncline == 10
                subplot(1,3,1)
            elseif treadmillIncline == 0
                subplot(1,3,2)
            elseif treadmillIncline == -10
                subplot(1,3,3)
            end
            if treadmillSpeed == 1.2
                cl = 'r';
            elseif treadmillSpeed == 1
                cl = 'b';
            elseif treadmillSpeed == 0.8
                cl = 'g';
            end
            if treadmillIncline == 10 || treadmillIncline == 0 || treadmillIncline == -10
                plot_j = plot3(avgPhase,stepLength_scaled*ones(size(avgPhase)), avgHeelPosUp, cl,'LineWidth',2);
                plot_j.Color(4) = 1;
            end
            
            
        end

    end
end

% The data in Gaitcyle always has 150 points per stride, and they are
% evenly spaced with respect to percent gait.


figure(2)

subplot(1,3,1)
title('Foot Angle, 10 deg incline')
view(60,60)
subplot(1,3,2)
title('0 deg incline')
view(60,60)
subplot(1,3,3)
title('-10 deg incline')
view(60,60)



figure(4)

subplot(1,3,1)
title('Shank Angle, 10 deg incline')
xlabel('percent gait')
ylabel('Step Length (m)')
zlabel('degrees'); 
view(60,60)
subplot(1,3,2)
title('0 deg incline')
xlabel('percent gait')
ylabel('Step Length (m)')
zlabel('degrees'); 
view(60,60)
subplot(1,3,3)
title('-10 deg incline')
xlabel('percent gait')
ylabel('Step Length (m)')
zlabel('degrees'); 
view(60,60)


figure(50)

subplot(1,3,1)
title('Heel Pos Forward, 10 deg incline')
xlabel('percent gait')
ylabel('Step Length (m)')
zlabel('degrees'); 
view(60,60)
subplot(1,3,2)
title('0 deg incline')
xlabel('percent gait')
ylabel('Step Length (m)')
zlabel('degrees'); 
view(60,60)
subplot(1,3,3)
title('-10 deg incline')
xlabel('percent gait')
ylabel('Step Length (m)')
zlabel('degrees'); 
view(60,60)

% figure(51)
% 
% subplot(1,3,1)
% title('Tibia Pos Forward, 10 deg incline')
% xlabel('percent gait')
% ylabel('Step Length (m)')
% zlabel('degrees'); 
% view(60,60)
% subplot(1,3,2)
% title('0 deg incline')
% xlabel('percent gait')
% ylabel('Step Length (m)')
% zlabel('degrees'); 
% view(60,60)
% subplot(1,3,3)
% title('-10 deg incline')
% xlabel('percent gait')
% ylabel('Step Length (m)')
% zlabel('degrees'); 
% view(60,60)

figure(52)

subplot(1,3,1)
title('Heel Pos Up, 10 deg incline')
xlabel('percent gait')
ylabel('Step Length (m)')
zlabel('degrees'); 
view(60,60)
subplot(1,3,2)
title('0 deg incline')
xlabel('percent gait')
ylabel('Step Length (m)')
zlabel('degrees'); 
view(60,60)
subplot(1,3,3)
title('-10 deg incline')
xlabel('percent gait')
ylabel('Step Length (m)')
zlabel('degrees'); 
view(60,60)


%% Regress

%C0 and %C1 continuity constraints at 0/1

phase = linspace(0,1,100)';
stride_lengths = linspace(0,2,3);
% stride_lengths = [0,1];
% ramps = [0,1];
ramps = linspace(-10,10,3);
A_eq_C0_constraint = [];
A_eq_C1_constraint = [];
for ramp = ramps
    
    A_eq_C0_ramp = [];
    A_eq_C1_ramp = [];
    for j= 1:length(stride_lengths)
        stride_length = stride_lengths(j);
        
        A_eq_fourier1 = [returnFourierBasis_Eval(0,stride_length,ramp, N_FOURIER) - returnFourierBasis_Eval(1,stride_length,ramp, N_FOURIER)];
        A_eq_fourier2 = [returnFourierBasis_DerivEval_dphase(0,stride_length,ramp, N_FOURIER) - returnFourierBasis_DerivEval_dphase(1,stride_length,ramp, N_FOURIER)];


        A_eq_C0_ramp = [A_eq_C0_ramp;A_eq_fourier1];
        A_eq_C1_ramp = [A_eq_C1_ramp;A_eq_fourier2];

    end
    
    A_eq_C0_constraint = [A_eq_C0_constraint;A_eq_C0_ramp];
    A_eq_C1_constraint = [A_eq_C1_constraint;A_eq_C1_ramp];

end   

[rows,~] = size(A_eq_C0_constraint);
b_eq_C0_constraint = zeros(rows,1);

[rows,~] = size(A_eq_C1_constraint);
b_eq_C1_constraint = zeros(rows,1);



% constraint to force constant behavior at sL = 0


phase = linspace(0,1,100)';
lenP = length(phase);


A_eq_fourier_constAtZerosL = [];
% ramps = linspace(-10,10,3);
ramps = [0,10,-10];
for ramp = ramps
    
    A_eq_fourier_ramp = zeros(length(phase), numFuncs);
    
    
    for j= 1:length(phase)
        
        phase_j = phase(j);
        
        fourier_deriv_coeffs = returnFourierBasis_Eval(phase_j,0,ramp, N_FOURIER);


        A_eq_fourier_ramp(j,:) = fourier_deriv_coeffs;

    end
    
    
    A_eq_fourier_constAtZerosL = [A_eq_fourier_constAtZerosL;A_eq_fourier_ramp];

end  

[rows,~] = size(A_eq_fourier_constAtZerosL);
b_eq_fourier_constAtZerosL = zeros(rows,1);

%foot
b_eq_const_foot_at_zero_sL = b_eq_fourier_constAtZerosL;
b_eq_const_foot_at_zero_sL(lenP+1:2*lenP) = 10;
b_eq_const_foot_at_zero_sL(2*lenP+1:3*lenP) = -10;

% additional constraints for foot


A_eq_zero_dsL_at_select_phase = [[returnFourierBasis_DerivEval_dsL(0.2,2,10, N_FOURIER)];...
    [returnFourierBasis_DerivEval_dsL(0.2,2,0, N_FOURIER)];...
    [returnFourierBasis_DerivEval_dsL(0.2,2,-10, N_FOURIER)]];




% force zero dsL at sL = 0
phase = linspace(0,1,100)';

A_eq_fourier_zerodsLatzerosL = [];
ramps = linspace(-10,10,3);
% ramps = [0,10,-10];
for ramp = ramps
    
    A_eq_fourier_shank_ramp = zeros(length(phase), numFuncs);
    for j= 1:length(phase)
        phase_j = phase(j);

        fourier_coeffs = returnFourierBasis_DerivEval_dsL(phase_j,0,ramp, N_FOURIER);

        
        A_eq_fourier_shank_ramp(j,:) = fourier_coeffs;

    end
    
    A_eq_fourier_zerodsLatzerosL = [A_eq_fourier_zerodsLatzerosL;A_eq_fourier_shank_ramp];

end   


A_eq_foot = [A_eq_C0_constraint; A_eq_C1_constraint; A_eq_fourier_constAtZerosL;A_eq_zero_dsL_at_select_phase;];
A_eq_shank = [A_eq_C0_constraint; A_eq_C1_constraint;A_eq_fourier_constAtZerosL;];
A_eq_heel_Pos_forward = [A_eq_C0_constraint; A_eq_C1_constraint;A_eq_fourier_constAtZerosL;];
A_eq_heel_Pos_up = [A_eq_C0_constraint; A_eq_C1_constraint;A_eq_fourier_constAtZerosL;];
A_eq_tibia_Pos_forward = [A_eq_C0_constraint; A_eq_C1_constraint;A_eq_fourier_constAtZerosL;];


% A_eq_foot = [A_eq_C0_constraint; A_eq_C1_constraint;A_eq_zero_dsL_at_select_phase;];
% A_eq_shank = [A_eq_C0_constraint; A_eq_C1_constraint;];
% A_eq_heel_Pos_forward = [A_eq_C0_constraint; A_eq_C1_constraint;];
% A_eq_tibia_Pos_forward = [A_eq_C0_constraint; A_eq_C1_constraint;];



b_eq_zero_dsL_at_select_phase = zeros(3,1);

%shank
b_eq_zero_shank_at_zero_sL = b_eq_fourier_constAtZerosL;

[rows,~] = size(A_eq_fourier_zerodsLatzerosL);
b_eq_zerodsLatzerosL = zeros(rows,1);

%heel Pos Y
b_eq_zero_heel_Pos_forward_at_zero_sL = b_eq_fourier_constAtZerosL;
b_eq_zero_heel_Pos_up_at_zero_sL = b_eq_fourier_constAtZerosL;
%tibia Pos Y
b_eq_zero_tibia_Pos_forward_at_zero_sL = b_eq_fourier_constAtZerosL;


b_eq_foot = [b_eq_C0_constraint; b_eq_C1_constraint; b_eq_const_foot_at_zero_sL;b_eq_zero_dsL_at_select_phase;];
b_eq_shank = [b_eq_C0_constraint; b_eq_C1_constraint;b_eq_zero_shank_at_zero_sL;];
b_eq_heel_Pos_forward = [b_eq_C0_constraint; b_eq_C1_constraint;b_eq_zero_heel_Pos_forward_at_zero_sL;];
b_eq_heel_Pos_up = [b_eq_C0_constraint; b_eq_C1_constraint;b_eq_zero_heel_Pos_up_at_zero_sL;];

b_eq_tibia_Pos_forward = [b_eq_C0_constraint; b_eq_C1_constraint;b_eq_zero_tibia_Pos_forward_at_zero_sL;];

% b_eq_foot = [b_eq_C0_constraint; b_eq_C1_constraint;b_eq_zero_dsL_at_select_phase;];
% b_eq_shank = [b_eq_C0_constraint; b_eq_C1_constraint;];
% b_eq_heel_Pos_forward = [b_eq_C0_constraint; b_eq_C1_constraint;];
% b_eq_tibia_Pos_forward = [b_eq_C0_constraint; b_eq_C1_constraint;];

% scale for regression against step length

best_fit_params_footAngle = lsqlin(A_mat_master,b_footAngle_master,[],[],A_eq_foot, b_eq_foot);
best_fit_params_shankAngle = lsqlin(A_mat_master,b_shankAngle_master,[],[],A_eq_shank, b_eq_shank);
best_fit_params_heelPosForward = lsqlin(A_mat_master,b_heelPosForward_master,[],[],A_eq_heel_Pos_forward, b_eq_heel_Pos_forward);
best_fit_params_heelPosUp = lsqlin(A_mat_master,b_heelPosUp_master,[],[],A_eq_heel_Pos_up, b_eq_heel_Pos_up);

best_fit_params_tibiaPosForward = lsqlin(A_mat_master,b_tibiaPosForward_master,[],[],A_eq_tibia_Pos_forward, b_eq_tibia_Pos_forward);



%% 
stepLength = 0:0.1:2;
phase = linspace(0,1,200)';

[X,Y] = meshgrid(phase,stepLength);

best_fit_footAngle_mat_311 = zeros(size(X));
best_fit_shankAngle_mat_311 = zeros(size(X));
best_fit_heelPosForward_mat_311 = zeros(size(X));
best_fit_heelPosUp_mat_311 = zeros(size(X));
best_fit_tibiaPosForward_mat_311 = zeros(size(X));

best_fit_footAngle_mat_312 = zeros(size(X));
best_fit_shankAngle_mat_312 = zeros(size(X));
best_fit_heelPosForward_mat_312 = zeros(size(X));
best_fit_heelPosUp_mat_312 = zeros(size(X));
best_fit_tibiaPosForward_mat_312 = zeros(size(X));

best_fit_footAngle_mat_313 = zeros(size(X));
best_fit_shankAngle_mat_313 = zeros(size(X));
best_fit_heelPosForward_mat_313 = zeros(size(X));
best_fit_heelPosUp_mat_313 = zeros(size(X));
best_fit_tibiaPosForward_mat_313 = zeros(size(X));

for i = 1:length(phase)
    
    for j = 1:length(stepLength)
%         i
%         j
        phase_i = X(j,i);
        stepLength_j = Y(j,i);
        best_fit_footAngle_mat_311(j,i) = best_fit_params_footAngle' * returnFourierBasis_Eval(phase_i,stepLength_j,10, N_FOURIER)';
        best_fit_shankAngle_mat_311(j,i) = best_fit_params_shankAngle' * returnFourierBasis_Eval(phase_i,stepLength_j,10, N_FOURIER)';
        best_fit_heelPosForward_mat_311(j,i) = best_fit_params_heelPosForward' * returnFourierBasis_Eval(phase_i,stepLength_j,10, N_FOURIER)';
        best_fit_heelPosUp_mat_311(j,i) = best_fit_params_heelPosUp' * returnFourierBasis_Eval(phase_i,stepLength_j,10, N_FOURIER)';        
        best_fit_tibiaPosForward_mat_311(j,i) = best_fit_params_tibiaPosForward' * returnFourierBasis_Eval(phase_i,stepLength_j,10, N_FOURIER)';

        best_fit_footAngle_mat_312(j,i) = best_fit_params_footAngle' * returnFourierBasis_Eval(phase_i,stepLength_j,0, N_FOURIER)';
        best_fit_shankAngle_mat_312(j,i) = best_fit_params_shankAngle' * returnFourierBasis_Eval(phase_i,stepLength_j,0, N_FOURIER)';
        best_fit_heelPosForward_mat_312(j,i) = best_fit_params_heelPosForward' * returnFourierBasis_Eval(phase_i,stepLength_j,0, N_FOURIER)';
        best_fit_heelPosUp_mat_312(j,i) = best_fit_params_heelPosUp' * returnFourierBasis_Eval(phase_i,stepLength_j,0, N_FOURIER)';
        best_fit_tibiaPosForward_mat_312(j,i) = best_fit_params_tibiaPosForward' * returnFourierBasis_Eval(phase_i,stepLength_j,0, N_FOURIER)';

        best_fit_footAngle_mat_313(j,i) = best_fit_params_footAngle' * returnFourierBasis_Eval(phase_i,stepLength_j,-10, N_FOURIER)';
        best_fit_shankAngle_mat_313(j,i) = best_fit_params_shankAngle' * returnFourierBasis_Eval(phase_i,stepLength_j,-10, N_FOURIER)';
        best_fit_heelPosForward_mat_313(j,i) = best_fit_params_heelPosForward' * returnFourierBasis_Eval(phase_i,stepLength_j,-10, N_FOURIER)';
        best_fit_heelPosUp_mat_313(j,i) = best_fit_params_heelPosUp' * returnFourierBasis_Eval(phase_i,stepLength_j,-10, N_FOURIER)';
        best_fit_tibiaPosForward_mat_313(j,i) = best_fit_params_tibiaPosForward' * returnFourierBasis_Eval(phase_i,stepLength_j,-10, N_FOURIER)';


    end
    
    
end

figure(2)
subplot(1,3,1)
surf(phase, stepLength,best_fit_footAngle_mat_311)
subplot(1,3,2)
surf(phase, stepLength,best_fit_footAngle_mat_312)
subplot(1,3,3)
surf(phase, stepLength,best_fit_footAngle_mat_313)

figure(4)
subplot(1,3,1)
surf(phase, stepLength,best_fit_shankAngle_mat_311)

subplot(1,3,2)
surf(phase, stepLength,best_fit_shankAngle_mat_312)
subplot(1,3,3)
surf(phase, stepLength,best_fit_shankAngle_mat_313)

figure(50)
subplot(1,3,1)
surf(phase, stepLength,best_fit_heelPosForward_mat_311)

subplot(1,3,2)
surf(phase, stepLength,best_fit_heelPosForward_mat_312)
subplot(1,3,3)
surf(phase, stepLength,best_fit_heelPosForward_mat_313)

% figure(51)
% subplot(1,3,1)
% surf(phase, stepLength,best_fit_tibiaPosForward_mat_311)
% 
% subplot(1,3,2)
% surf(phase, stepLength,best_fit_tibiaPosForward_mat_312)
% subplot(1,3,3)
% surf(phase, stepLength,best_fit_tibiaPosForward_mat_313)

figure(52)
subplot(1,3,1)
surf(phase, stepLength,best_fit_heelPosUp_mat_311)

subplot(1,3,2)
surf(phase, stepLength,best_fit_heelPosUp_mat_312)
subplot(1,3,3)
surf(phase, stepLength,best_fit_heelPosUp_mat_313)


%% plot against incline

incline = -10:0.1:10;
stepLength_k = 1;


[X,Y] = meshgrid(phase,incline);

best_fit_footAngle_mat_incline = zeros(size(X));
best_fit_shankAngle_mat_incline = zeros(size(X));
best_fit_heelPosForward_mat_incline = zeros(size(X));
best_fit_heelPosUp_mat_incline = zeros(size(X));
best_fit_tibiaPosForward_mat_incline = zeros(size(X));

for i = 1:length(phase)
    
    for j = 1:length(incline)
%         i
%         j
        phase_i = X(j,i);
        incline_j = Y(j,i);
        best_fit_footAngle_mat_incline(j,i) = best_fit_params_footAngle' * returnFourierBasis_Eval(phase_i,stepLength_k,incline_j, N_FOURIER)';
        best_fit_shankAngle_mat_incline(j,i) = best_fit_params_shankAngle' * returnFourierBasis_Eval(phase_i,stepLength_k,incline_j, N_FOURIER)';
        best_fit_heelPosForward_mat_incline(j,i) = best_fit_params_heelPosForward' * returnFourierBasis_Eval(phase_i,stepLength_k,incline_j, N_FOURIER)';
        best_fit_heelPosUp_mat_incline(j,i) = best_fit_params_heelPosUp' * returnFourierBasis_Eval(phase_i,stepLength_k,incline_j, N_FOURIER)';
        best_fit_tibiaPosForward_mat_incline(j,i) = best_fit_params_tibiaPosForward' * returnFourierBasis_Eval(phase_i,stepLength_k,incline_j, N_FOURIER)';


    end
    
    
end



figure(5)
surf(phase, incline,best_fit_footAngle_mat_incline)
title('Foot Angle')

xlabel('percent gait')
ylabel('Incline, scaled (deg)')
zlabel('degrees'); 
view(60,60)

figure(6)
surf(phase, incline,best_fit_shankAngle_mat_incline)
title('Shank Angle')

xlabel('percent gait')
ylabel('Incline, scaled (deg)')
zlabel('degrees'); 
view(60,60)

figure(7)
surf(phase, incline,best_fit_heelPosForward_mat_incline)
title('Heel Pos Forward')

xlabel('percent gait')
ylabel('Incline, scaled (deg)')
zlabel('degrees'); 
view(60,60)

% figure(8)
% surf(phase, incline,best_fit_tibiaPosForward_mat_incline)
% title('Tibia Pos Forward')
% 
% xlabel('percent gait')
% ylabel('Incline, scaled (deg)')
% zlabel('degrees'); 
% view(60,60)

figure(9)
surf(phase, incline,best_fit_heelPosUp_mat_incline)
title('Heel Pos Up')

xlabel('percent gait')
ylabel('Incline, scaled (deg)')
zlabel('degrees'); 
view(60,60)

%% save variables

save('regressionMatrices_gaitModel_fourier_normalizedsL_linearsL','A_mat_master',...
    'b_footAngle_master','b_shankAngle_master','best_fit_params_footAngle',...
    'best_fit_params_shankAngle','best_fit_params_heelPosForward','best_fit_params_tibiaPosForward','best_fit_params_heelPosUp')

M = [best_fit_params_footAngle';...
    best_fit_params_shankAngle';...
    best_fit_params_heelPosForward';...
    best_fit_params_heelPosUp'];
writematrix(M,'gaitModel_fourier_normalizedsL_linearsL.csv')
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
