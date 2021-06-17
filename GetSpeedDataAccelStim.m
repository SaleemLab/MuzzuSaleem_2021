%% Tomaso Muzzu - UCL - 15/01/2018

% function to prepare data of speed profile of mice during visual stimuli
% presentaiton of grating with changing contrast and possible TF
% perturbation.

function [Array, StimulusParams] = GetSpeedDataAccelStim(ES,varargin)

pnames = {'Norm_mode','log_speed','Smoothing'};

Dflt_Values = {'none', 0, 0.6};

[NormMode, logSpeed, GaussFilter_W] ...
    = internal.stats.parseArgs(pnames,Dflt_Values,varargin{:});

% Default inputs:
%%% NormMode = 'none'; do not normalise the speed. Other options are:
% ------ 'MaxSpeedTrial'= normalise by max speed of trial
% ------ 'MaxSpeedSession'= normalise by max speed of session
% ------ 'z_score'= normalise by z score (X-mean)/std of session
%%% logSpeed = 1; % compute speed in log units. if 0 -> normal units.
%%% Smoothing = 0.2; % width of smoothing window in seconds

% Gaussian filter
BonvisionFR = (ES.Time(end)-ES.Time(1))/length(ES.Time); 
Width = round(GaussFilter_W/BonvisionFR);
Sigma = Width/3; % standard deviation in number of samples (converted from time in seconds)
x_g = linspace(-Width/2, Width/2, Width);
gaussFilter = exp(-x_g.^2/(2*Sigma^2));
gaussFilter_ = gaussFilter / sum (gaussFilter); % normalize
if sum(diff(ES.MousePosition))<0
    SmoothedSpeed = conv(-ES.MouseSpeed, gaussFilter_, 'same');
else 
    SmoothedSpeed = conv(ES.MouseSpeed, gaussFilter_, 'same');
end
SmoothedSpeed(SmoothedSpeed<=0) = 0.01;
% normalisation
switch NormMode 
    case 'none'
        SmoothedSpeed = SmoothedSpeed;
    case 'MaxSpeedSession'
        SmoothedSpeed = SmoothedSpeed/max(SmoothedSpeed);
    case 'z_score'
        SmoothedSpeed = (SmoothedSpeed - mean(SmoothedSpeed))/std(SmoothedSpeed);
    otherwise
        SmoothedSpeed = SmoothedSpeed;
end
% log scale
if logSpeed == 1
    SmoothedSpeed = log(SmoothedSpeed);
end
% save speed, time, TF, and contrast for every trial
trialSide_seconds = 1; % take 60 samples before and after trial
trialSide_samples = trialSide_seconds*round(length(ES.Time)/(max(ES.Time)-min(ES.Time))); % take 60 samples before and after trial
clear SpeedTrial_temp SpeedTrial SpeedTrial_PertStartEnds clear StimulusParams
r = 1;
avg_trial_duration = round(mean(ES.trialStartsEnds(:,2)-ES.trialStartsEnds(:,1)));
for j = 1:size(ES.trialStartsEnds,1)   % scan through the trials
    [v TrialStart_i] = min(abs(ES.Time-ES.trialStartsEnds(j,1)));
    [v TrialEnd_i] = min(abs(ES.Time-ES.trialStartsEnds(j,2)));
    % speed
    clear SpeedTrial_temp
    SpeedTrial_temp = SmoothedSpeed(max(TrialStart_i-trialSide_samples,1):min(length(SmoothedSpeed),TrialEnd_i+trialSide_samples));
    if strcmp(NormMode, 'MaxSpeedTrial')
        SpeedTrial_temp = SpeedTrial_temp/max(SpeedTrial_temp);
    end
    if length(SpeedTrial_temp)<ceil(5/BonvisionFR)
        SpeedTrial_temp = [zeros(ceil(5/BonvisionFR)-length(SpeedTrial_temp),1) ; SpeedTrial_temp];
    end
    % time
    time = linspace(ES.trialStartsEnds(j,1)-trialSide_seconds,ES.trialStartsEnds(j,2)+trialSide_seconds,ceil((avg_trial_duration/BonvisionFR)+(2/BonvisionFR)))'-ES.trialStartsEnds(j,1);
    % temporal frequency
    clear TF
    TF = ES.TF(max(TrialStart_i,1):min(length(ES.TF),TrialEnd_i));
    TF = [zeros(round((ceil(5/BonvisionFR)-length(TF))/2),1); TF];
    TF = [TF; zeros(ceil(5/BonvisionFR)-length(TF),1)];
    
    % contrast
    clear contrast
    contrast = ES.Contrast(max(TrialStart_i,1):min(length(ES.Contrast),TrialEnd_i));
    contrast = [zeros(round((ceil(5/BonvisionFR)-length(TF))/2),1); contrast];
    contrast = [contrast; zeros(ceil(5/BonvisionFR)-length(contrast),1)];
    
    StimulusParams{j,:}= [SpeedTrial_temp(1:ceil(5/BonvisionFR)) time TF contrast];
end
   

% save array with info necessary to plot everything
clear Array
cc = 1;
Array = table(1,1,1,1,[1 1],1,{1},{1},{1},{1},...
        'VariableNames',{'Trial_ID';'Direction';'Starting_TF';'Ending_TF';'Accel_start_end';'dTF';'Speed_in_trial';'time_in_trial'; 'TF_in_trial'; 'contrast_in_trial'});
for i = 1:size(StimulusParams,1)
    Array{i,1} = i; % order of display
    Array{i,2} = 0; % grating direction
    Array{i,3} = mean(StimulusParams{i,:}(trialSide_samples+10:trialSide_samples+20,3));% starting TF
    Array{i,4} = mean(StimulusParams{i,:}(end-trialSide_samples-20:end-trialSide_samples-10,3));% starting TF% ending TF
    if Array{i,3}~=Array{i,4}
        Array{i,5} = [ min(find((diff(StimulusParams{i,:}(trialSide_samples+10:end-trialSide_samples-10,3)))))+trialSide_samples+10 , ...
                       max(find((diff(StimulusParams{i,:}(trialSide_samples+10:end-trialSide_samples-10,3)))))+trialSide_samples+10 ]; % indexes at which TF starts and stops changing
    else
        Array{i,5} = Array{i-1,5};
    end
    Array{i,6} = round((Array{i,4}-Array{i,3})/(diff(Array{i,5})*BonvisionFR)*abs(Array{i,3})/Array{i,3},0); % acceleration
    Array{i,7} = {StimulusParams{i,:}(:,1)}; % speed of every trial
    Array{i,8} = {StimulusParams{i,:}(:,2)}; % time of every trial
    Array{i,9} = {StimulusParams{i,:}(:,3)}; % TF of every trial
    Array{i,10} = {StimulusParams{i,:}(:,4)}; % contrast of every trial
end



%EOF
end
