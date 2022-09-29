close all;clc;clearvars -except Gaitcycle Continuous %removes all variables except for gaitcyle and continuous from the matlab workspace. Useful because loading these variables can take several minutes.
%load the InclineExperiment data from the folder location you specify
% load('Z:\your_file_location_here\InclineExperiment.mat') 
%This file is large and may take several minutes to load. We recommend not
%clearing the variable, and only loading it as many times as necessary.

%% Example: removing strides with outliers from kinematic data, and taking the mean

% load('regressionMatrices_dataport_3P.mat')
load('regressionMatrices_dataport_3P_normalizedsL.mat')


% Gaitcyle and Continuous store the InclineExperiment data in a MATLAB
% structure array, or 'struct'

% The fields of a struct can be iterated over using a cell array

% A cell array can be defined manually with this notation: 
sub={'AB01','AB02','AB03','AB04','AB05','AB06','AB07','AB08','AB09','AB10'};

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

c = {'r','g','b'};
% c1 = {'b','b','b','b','b','r','r','r','r'}
% 
% c2 = [0.5,0.4,0.3,0.2,0.1,0.2,0.3,0.4,0.5]

% Initializing a matrix that will contain the mean joint data of all subjects,
%all tasks
% concatenated_data=NaN(150,numel(sub)*numel(trials));


phaseDelins = [0.1,0.5,0.65,1];

numInclineFuncs = 2;
numStepLengthFuncs = 3;
numPhaseFuncs = 4;
numFuncs = numInclineFuncs*numStepLengthFuncs*numPhaseFuncs;

phase = linspace(0,1,150);

A_covar_mat_master = zeros(150,1);
b_covar_11_master = zeros(150,1);
b_covar_12_master = zeros(150,1);
b_covar_13_master = zeros(150,1);
b_covar_14_master = zeros(150,1);

b_covar_22_master = zeros(150,1);
b_covar_23_master = zeros(150,1);
b_covar_24_master = zeros(150,1);

b_covar_33_master = zeros(150,1);
b_covar_34_master = zeros(150,1);

b_covar_44_master = zeros(150,1);


isPosDef = zeros(150,1);

for phase_idx = 1:150
    
    N = 0;
    X_meas = [];

    
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


       	subjLegLength = Gaitcycle.(sub{i}).subjectdetails{5,2}/1000;

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
            [treadmillSpeed, treadmillIncline] = returnSpeedIncline(trial);
    %         incline_scaled = (treadmillIncline - inclineMin)/(inclineMax - inclineMin);

            incline_scaled = treadmillIncline;


            

            [rows,cols] = size(foot_angle_data);
            
            for k = 1:cols
                N = N + 1;

                foot_angle_data_col = foot_angle_data(:,k);
                shank_angle_data_col = shank_angle_data(:,k);


                
                
                time_data_col = time_data(:,k);
                
                
                % velocities
                foot_vel_data_col = diff(foot_angle_data_col)./diff(time_data_col);
                foot_vel_data_col = [foot_vel_data_col(1);foot_vel_data_col];
                
%                 foot_vel_data_col = lowpass(foot_vel_data_col,25,1/mean(diff(time_data_col)));

                shank_vel_data_col = diff(shank_angle_data_col)./diff(time_data_col);
                shank_vel_data_col = [shank_vel_data_col(1);shank_vel_data_col];

%                 shank_vel_data_col = lowpass(shank_vel_data_col,25,1/mean(diff(time_data_col)));
            
            

                stepDuration = time_data_col(end) - time_data_col(1);
                phaseDot = 1/stepDuration;
                stepLength = treadmillSpeed*stepDuration;

    %             stepLength_scaled = (stepLength - stepMin)/(stepMax - stepMin);

                stepLength_scaled = stepLength/subjLegLength;


                phase_col = time_data_col/stepDuration;


%                 X_meas = [];

%                     mat_idx = phase_idx + (k - 1)*150;
                phase_j = phase_col(phase_idx);
                
                footAngle_measured = foot_angle_data_col(phase_idx);
                shankAngle_measured = shank_angle_data_col(phase_idx);
                
                footAngleVel_measured = foot_vel_data_col(phase_idx);
                shankAngleVel_measured = shank_vel_data_col(phase_idx);

        %         bezier_data = [((1 - phase_i).^3), (3 * ((1 - phase_i).^2 ).*(phase_i) ), ...
        %             (3 * ((1 - phase_i) ).*((phase_i).^2) ), ((phase_i).^3)];

                footAngle_model = returnPiecewiseBezier3D(phase_j,stepLength_scaled,incline_scaled,best_fit_params_footAngle, phaseDelins,numFuncs);
                shankAngle_model = returnPiecewiseBezier3D(phase_j,stepLength_scaled,incline_scaled,best_fit_params_shankAngle, phaseDelins,numFuncs);

                footAngleVel_model = phaseDot * returnPiecewiseBezier3DDeriv_dphase(phase_j,stepLength_scaled,incline_scaled,best_fit_params_footAngle, phaseDelins,numFuncs);
                shankAngleVel_model = phaseDot * returnPiecewiseBezier3DDeriv_dphase(phase_j,stepLength_scaled,incline_scaled,best_fit_params_shankAngle, phaseDelins,numFuncs);


                
                
                footAngle_residual = footAngle_measured - footAngle_model;
                shankAngle_residual = shankAngle_measured - shankAngle_model;
                
                footAngleVel_residual = footAngleVel_measured - footAngleVel_model;
                shankAngleVel_residual = shankAngleVel_measured - shankAngleVel_model;
                
                X_meas = [X_meas;[footAngle_residual, footAngleVel_residual, shankAngle_residual,shankAngleVel_residual]];



            end

        end
    end
    
    N;
    C = ( 1./(N - 1)) * (X_meas' * X_meas);
    
%     C = C + eye(2);
    
    
%     pause
    
    
    b_covar_11_master(phase_idx) = C(1,1);
    b_covar_12_master(phase_idx) = C(1,2);
    b_covar_13_master(phase_idx) = C(1,3);
    b_covar_14_master(phase_idx) = C(1,4);
    
    b_covar_22_master(phase_idx) = C(2,2);
    b_covar_23_master(phase_idx) = C(2,3);
    b_covar_24_master(phase_idx) = C(2,4);
    
    b_covar_33_master(phase_idx) = C(3,3);
    b_covar_34_master(phase_idx) = C(3,4);
    
    b_covar_44_master(phase_idx) = C(4,4);
    
    [~,p] = chol(C);
    
    isPosDef(phase_idx) = p == 0;

end


%%

% 
% b_covar_12_master = lowpass(b_covar_12_master,10,120);
% b_covar_14_master = lowpass(b_covar_14_master,10,120);
% b_covar_22_master = lowpass(b_covar_22_master,10,120);
% b_covar_24_master = lowpass(b_covar_24_master,10,120);
% b_covar_34_master = lowpass(b_covar_34_master,10,120);
% 
% b_covar_44_master = lowpass(b_covar_44_master,10,120);



figure(1)
subplot(4,1,1)
semilogy(phase, b_covar_11_master,'LineWidth',2)
hold on
legend('11')


subplot(4,1,2)
semilogy(phase, b_covar_22_master,'LineWidth',2)
hold on
legend('22')

subplot(4,1,3)
semilogy(phase, b_covar_33_master,'LineWidth',2)
hold on
legend('33')


subplot(4,1,4)
semilogy(phase, b_covar_44_master,'LineWidth',2)
hold on
legend('44')


figure(2)

subplot(4,1,1)
plot(phase, b_covar_11_master,'LineWidth',2)
hold on
legend('11')


subplot(4,1,2)
plot(phase, b_covar_12_master,'LineWidth',2)
hold on
legend('12')

subplot(4,1,3)
plot(phase, b_covar_13_master,'LineWidth',2)
hold on
legend('13')


subplot(4,1,4)
plot(phase, b_covar_14_master,'LineWidth',2)
hold on
legend('14')


figure(3)

subplot(3,1,1)
plot(phase, b_covar_22_master,'LineWidth',2)
hold on
legend('22')


subplot(3,1,2)
plot(phase, b_covar_23_master,'LineWidth',2)
hold on
legend('23')

subplot(3,1,3)
plot(phase, b_covar_24_master,'LineWidth',2)
hold on
legend('24')

figure(4)

subplot(2,1,1)
plot(phase, b_covar_33_master,'LineWidth',2)
hold on
legend('33')


subplot(2,1,2)
plot(phase, b_covar_34_master,'LineWidth',2)
hold on
legend('34')


% subplot(3,1,2)
% plot(phase, b_covar_12_master,'LineWidth',2)
% hold on
% legend('12')


figure(20)
plot(phase, isPosDef,'LineWidth',2)
legend('PD test')



%% regress

% bezier_data_mat_master = [];
% 
% bezier_data_mat = [];
%     
% for j = 1:length(phase)
%     phase_j = phase(j);
% 
% 
%     bezier_data = returnBezierCoeffs(phase_j);
% 
%     if phase_j <= phaseDelins(1)
%         bezier_row = [bezier_data, zeros(1,3*4)];
%     elseif phase_j <= phaseDelins(2)
%         bezier_row = [zeros(1,4),bezier_data, zeros(1,2*4)];
%     elseif phase_j <= phaseDelins(3)
%         bezier_row = [zeros(1,2*4),bezier_data, zeros(1,4)];
%     elseif phase_j <= phaseDelins(4)
%         bezier_row = [zeros(1,3*4),bezier_data];
% 
%     end
% 
% 
%     bezier_data_mat_master = [bezier_data_mat_master; bezier_row];
% 
% 
% end
%     
%     
%     
% 
%     
%     
% A = bezier_data_mat_master;
% 
% % A_eq = [[0,0,0,1,-1,0,0,0,zeros(1,4),zeros(1,4)];...
% %     [0,0,-3,3,3,-3,0,0,zeros(1,4),zeros(1,4)];...
% %     [zeros(1,4),0,0,0,1,-1,0,0,0,zeros(1,4)];...
% %     [zeros(1,4),0,0,-3,3,3,-3,0,0,zeros(1,4)];...
% %     [zeros(1,4),zeros(1,4),0,0,0,1,-1,0,0,0];...
% %     [zeros(1,4),zeros(1,4),0,0,-3,3,3,-3,0,0];...
% %     [-1,0,0,0,zeros(1,4),zeros(1,4),0,0,0,1];...
% %     [3,-3,0,0,zeros(1,4),zeros(1,4),0,0,-3,3]]
% 
% 
% A_eq = [[returnBezierCoeffs(phaseDelins(1)), -returnBezierCoeffs(phaseDelins(1)), zeros(1,4), zeros(1,4)];...
%     [returnBezierDerivCoeffs(phaseDelins(1)), -returnBezierDerivCoeffs(phaseDelins(1)), zeros(1,4), zeros(1,4)];...
%     [zeros(1,4), returnBezierCoeffs(phaseDelins(2)), -returnBezierCoeffs(phaseDelins(2)), zeros(1,4)];...
%     [zeros(1,4), returnBezierDerivCoeffs(phaseDelins(2)), -returnBezierDerivCoeffs(phaseDelins(2)), zeros(1,4)];...
%     [zeros(1,4), zeros(1,4), returnBezierCoeffs(phaseDelins(3)), -returnBezierCoeffs(phaseDelins(3))];...
%     [zeros(1,4), zeros(1,4), returnBezierDerivCoeffs(phaseDelins(3)), -returnBezierDerivCoeffs(phaseDelins(3))];...
%     [-returnBezierCoeffs(0), zeros(1,4), zeros(1,4), returnBezierCoeffs(phaseDelins(4))];...
%     [-returnBezierDerivCoeffs(0), zeros(1,4), zeros(1,4), returnBezierDerivCoeffs(phaseDelins(4))]];
% 
% 
% b_eq = zeros(8,1);
% 
% 
% best_fit_params_11 = lsqlin(A,b_covar_11_master,[],[],A_eq, b_eq);
% best_fit_params_12 = lsqlin(A,b_covar_12_master,[],[],A_eq, b_eq);
% best_fit_params_22 = lsqlin(A,b_covar_22_master,[],[],A_eq, b_eq);

% N = 0:50;
% 
% A_mat = zeros(length(phase),length(N));
% 
% for i = 1:length(N)
%     
%     A_mat(:,i) = cos(2*pi*phase*N(i))';
%     
% end
% 
% best_fit_params_11 = (A_mat' * A_mat)^-1 * A_mat' * b_covar_11_master;
% best_fit_params_12 = (A_mat' * A_mat)^-1 * A_mat' * b_covar_12_master;
% best_fit_params_22 = (A_mat' * A_mat)^-1 * A_mat' * b_covar_22_master;

% best_fit_11 = zeros(size(phase));
% best_fit_12 = zeros(size(phase));
% best_fit_22 = zeros(size(phase));

% for i = 1:length(phase)
%     
%     best_fit_11(i) = returnPiecewiseBezier(phase(i),best_fit_params_11, phaseDelins);
%     best_fit_12(i) = returnPiecewiseBezier(phase(i),best_fit_params_12, phaseDelins);
%     best_fit_22(i) = returnPiecewiseBezier(phase(i),best_fit_params_22, phaseDelins);
%     
% end

% best_fit_11 = A_mat * best_fit_params_11;
% best_fit_12 = A_mat * best_fit_params_12;
% best_fit_22 = A_mat * best_fit_params_22;


% best_fit_params_11 = polyfit(phase', b_covar_11_master,12);
% best_fit_params_12 = polyfit(phase', b_covar_12_master,12);
% best_fit_params_22 = polyfit(phase', b_covar_22_master,12);
% 
% best_fit_11 = polyval(best_fit_params_11,phase);
% best_fit_12 = polyval(best_fit_params_12,phase);
% best_fit_22 = polyval(best_fit_params_22,phase);


best_fit_params_11 = b_covar_11_master;
best_fit_params_12 = b_covar_12_master;
best_fit_params_13 = b_covar_13_master;
best_fit_params_14 = b_covar_14_master;

best_fit_params_22 = b_covar_22_master;
best_fit_params_23 = b_covar_23_master;
best_fit_params_24 = b_covar_24_master;

best_fit_params_33 = b_covar_33_master;
best_fit_params_34 = b_covar_34_master;

best_fit_params_44 = b_covar_44_master;


best_fit_11 = interp1(phase',b_covar_11_master,phase');
best_fit_12 = interp1(phase',b_covar_12_master,phase');
best_fit_22 = interp1(phase',b_covar_22_master,phase');


% figure(1)
% subplot(3,1,1)
% plot(phase, best_fit_11,'LineWidth',2)
% legend('11','fit')
% 
% subplot(3,1,2)
% plot(phase, best_fit_12,'LineWidth',2)
% legend('12','fit')
% 
% subplot(3,1,3)
% plot(phase, best_fit_22,'LineWidth',2)
% legend('22','fit')

%%
save('covar_best_fit','best_fit_params_11','best_fit_params_12','best_fit_params_22')
M = [best_fit_params_11';best_fit_params_12';best_fit_params_13';best_fit_params_14';...
    best_fit_params_22';best_fit_params_23';best_fit_params_24';...
    best_fit_params_33';best_fit_params_34';
    best_fit_params_44';];
writematrix(M,'covar_best_fit_normalizedsL.csv')

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
