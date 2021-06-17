%% Tomaso Muzzu - UCL - 16/09/2019

% function to select running trials

function [AM_SpeedResponse AM_SpeedResponseControl] = ComputeSpeedResponse(AM_Param, AM_Speed,SelectedCells,Trials_PertON,Trials_PertOFF,ProjectData)

trialSide_seconds = 1; % take 60 samples before and after trial
trialSide_samples = 60; % take 60 samples before and after trial
BonVision_SR = 60;

% generate 3D matrix for including time responses
t1 = double(repmat(Trials_PertON, 1, 1, size(AM_Speed,3))); % repeating the selection across time
t2 =  double(repmat(Trials_PertOFF, 1, 1, size(AM_Speed,3))); % repeating the selection across time
% nan all the zero values
t1(t1(:)==0) = nan; t2(t2(:)==0) = nan;

% select only the speed traces during selected trials (pert ON and OFF separately)
Trials_Speed_OI = t1.*AM_Speed;
TrialsControl_Speed_OI = t2.*AM_Speed;

% flatten the speed traces onto a 2D matrix while keeping info of the
% recording number and of the unit (not necessary though)
Trials_Speed_OI_2D = reshape(Trials_Speed_OI, size(Trials_Speed_OI,1)*size(Trials_Speed_OI,2),size(Trials_Speed_OI,3));
TrialsControl_Speed_OI_2D = reshape(TrialsControl_Speed_OI, size(TrialsControl_Speed_OI,1)*size(TrialsControl_Speed_OI,2),size(TrialsControl_Speed_OI,3));
RecordingID = reshape(AM_Param(:,:,1),size(AM_Param,1)*size(AM_Param,2),1);

%% compute the mean speed at different point in times for these trials in
% perturbation trials
TimeLine = linspace(min(ProjectData.Trial_data{1,1}.time_in_trial{1}), ...
                    max(ProjectData.Trial_data{1,1}.time_in_trial{1}),...
                    size(AM_Speed,3))-trialSide_seconds;
% -1s to 2s of visual stim. onset
AM_SpeedResponse_t(:,1) = mean(Trials_Speed_OI_2D(:,trialSide_samples:trialSide_samples+2*BonVision_SR),2);
% -0.5s to end of perturbation period
PertOnsets_2D = AM_Param(:,:,4);
PertOnsets_1D = PertOnsets_2D(:);
[v ind(1)] = min(abs(TimeLine-min(PertOnsets_1D)));
[v ind(2)] = min(abs(TimeLine-max(PertOnsets_1D)));
% -0.5s to +1s at perturbation offset time
PertOffsets_2D = AM_Param(:,:,5);
PertOffsets_1D = PertOffsets_2D(:);
[v ind2(1)] = min(abs(TimeLine-min(PertOffsets_1D)));
[v ind2(2)] = min(abs(TimeLine-max(PertOffsets_1D)));
AM_SpeedResponse_t(:,2) = mean(Trials_Speed_OI_2D(:,ind(1)-trialSide_samples:ind2(2)),2);
AM_SpeedResponse_t(:,3) = mean(Trials_Speed_OI_2D(:,ind2(1)-trialSide_samples:ind2(2)+trialSide_samples),2);
% avg speed for the entire trial
AM_SpeedResponse_t(:,4) = mean(Trials_Speed_OI_2D(:,trialSide_samples:end-trialSide_samples),2);

%% compute the mean speed at different point in times for these trials in
% perturbation trials
AM_SpeedResponseControl_t(:,1) = mean(TrialsControl_Speed_OI_2D(:,trialSide_samples:trialSide_samples+2*BonVision_SR),2);
AM_SpeedResponseControl_t(:,2) = mean(TrialsControl_Speed_OI_2D(:,ind(1)-trialSide_samples:ind2(2)),2);
AM_SpeedResponseControl_t(:,3) = mean(TrialsControl_Speed_OI_2D(:,ind2(1)-trialSide_samples:ind2(2)+trialSide_samples),2);
AM_SpeedResponseControl_t(:,4) = mean(TrialsControl_Speed_OI_2D(:,trialSide_samples:end-trialSide_samples),2);

%% select only recordings where there are at least 10 trials for each condition (running and stationary)
for i = 1:size(AM_SpeedResponse_t,2)
    AM_SpeedResponse(:,:,i) = reshape(AM_SpeedResponse_t(:,i),size(Trials_Speed_OI,1),size(Trials_Speed_OI,2),1);
    AM_SpeedResponseControl(:,:,i) = reshape(AM_SpeedResponseControl_t(:,i),size(Trials_Speed_OI,1),size(Trials_Speed_OI,2),1);
end


end





