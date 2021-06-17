%% Tomaso Muzzu - UCL - 17/09/2019

% function to select only trials with pre-specified grating directions


function [AM_AngleResponse AM_AngleResponseControl] = ComputeAngleResponse(AM_Param, AM_Speed,SelectedCells,Trials_PertON,Trials_PertOFF,ProjectData, Angles_OI);


trialSide_seconds = 1; % take 60 samples before and after trial
trialSide_samples = 60; % take 60 samples before and after trial
BonVision_SR = 60;

for i = 1 : size(Angles_OI,2)
    % find all PertON trials with the grating direction of interest
    AngleFilter(:,:,i) =  AM_Param(:,:,2)==Angles_OI(i);
    AM_AngleResponse(:,:,i) = AngleFilter(:,:,i) & Trials_PertON;
    % find all PertOFF trials with the grating direction of interest
    AM_AngleResponseControl(:,:,i) = AngleFilter(:,:,i) & Trials_PertOFF;
end



end