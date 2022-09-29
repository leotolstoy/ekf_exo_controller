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

phase_rates = [];
stride_lengths = [];
speeds = [];

figure(1)
hold on

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
        time_data = Gaitcycle.(sub{i}).(trials{j}).cycles.(leg{1}).time;
        
        % delete the strides identified to contain outliers
        time_data(:,Gaitcycle.(sub{i}).(trials{j}).stepsout.(leg{1})) = [];
        incline_scaled = treadmillIncline;
        [~,cols] = size(time_data);


        for k = 1:cols
            
            time_data_col = time_data(:,k);
            stepDuration = time_data_col(end) - time_data_col(1);
            phaseDot = 1/stepDuration;
            stepLength = treadmillSpeed*stepDuration;
            
            
            stepLength_scaled = stepLength/subjLegLength;
            phase_col = time_data_col/stepDuration;

            phase_rates = [phase_rates;phaseDot];
            stride_lengths = [stride_lengths;stepLength_scaled];
            speeds = [speeds,treadmillSpeed];

            
        end

        
            
            
    end

end

A = [ones(length(phase_rates),1),phase_rates];
b = stride_lengths;

phi = (A' * A)^-1 * A' * b
% phi = lsqr(A,b)
k = phi(1);
lambda = phi(2)

xx = linspace(0,max(phase_rates));
% xx = phase_rates;
yy = lambda * xx + k;
% yy = A * phi;


isSlow = speeds == 0.8;
isNorm = speeds == 1.0;
isFast = speeds == 1.2;
figure(1)
plot(phase_rates(isSlow), stride_lengths(isSlow),'go')
hold on
plot(phase_rates(isNorm), stride_lengths(isNorm),'bo')
plot(phase_rates(isFast), stride_lengths(isFast),'ro')
plot(xx,yy,'k--')
xlabel('Phase Rate')
ylabel('Stride Length')
ylim([0,2])
grid minor



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
