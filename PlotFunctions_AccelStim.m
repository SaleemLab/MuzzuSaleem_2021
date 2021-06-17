%% Tomaso Muzzu - UCL - 24/09/2019

%% Load raw and formatted data
if ~exist('ProjectData','var')
    [ProjectData AM_UnitResponses AM_Param AM_Speed AM_UOI] = LoadData;
end

%% Initial selection of best units
% to be run once in theory, repeated in case params need changing
FR_thres = 0.1;
AM_UOI = UnitsSelection(ProjectData,AM_UnitResponses,FR_thres);
% save([ DataFolder filesep 'UnitsResponse.mat'], 'AM_UnitResponses', 'AM_Speed', 'AM_Param', 'AM_UOI','-v7.3')
% ATM, units with mean FR>0.1Hz and comparable FR in first and last third
% of recording are kept

%% smooth data with Gaussian filter
% 3rd argument is width of Gaussian filter
% 4th OPTIONAL argument is width of std - sigma=1/3 width by defeault;
Gauss_width = 0.3; % in seconds
AM_UnitResponses_smooth = Smooth_AM_FR(ProjectData,AM_UnitResponses,Gauss_width);

%% apply conditions to select responses
% conditions are set to define subset of units (AM_UOI), trials, periods within
% trials based on info store in matrix AM_Param(possible trials, units,
% conditions) where last dimension is 
% conditions = 1 --> nr of recording
% conditions = 2 --> TF at stimulus onset
% conditions = 3 --> TF at simutlus offset
% conditions = 4 --> dTF shown
% conditions = 5 --> index of when Tf starts changing
% conditions = 6 --> index of when Tf finishes changing



%Cond_Sel = {} ; % {optional('run'), optional(grating dir=[0:45:315])}
Units_Sel = AM_UOI(:,1) ; % this could include further selection criteria steming from other analysis
% Units_Sel = logical(ones(size(AM_UOI(:,1),1),1));
Run_TH = 3; % running speed threshold is 1 cm/s to be used like this: ..., 'run', Run_TH)
% Cond = 0; % %% all units together, distinguishing for accel. values
Cond = 1; % %% all units together, distringuishing for accel. values and starting TF
%Cond = 2; % %% plot single units, distringuishing for accel. values
%Cond = 3; % %% plot single units, distringuishing for accel. values and starting TF
SelectedResponses = AccelResponsesSelection(ProjectData,AM_UnitResponses_smooth,AM_Param,AM_Speed(:,:,1:end-1),Units_Sel,Cond)%,'run',Run_TH);
Units_Sel= true(449,1);
Units_Sel(UnitPerRec(7)+1:UnitPerRec(8)) = 0;
Units_Sel(UnitPerRec(9)+1:UnitPerRec(10)) = 0;

%% plot the responses
switch Cond 
    case 0
        % Is there a noticeable population response to acceleration?
        PlotResponses_AccelStim_0(ProjectData,SelectedResponses,AM_Param,AM_Speed,Units_Sel,0);
    case 1
        % Does the initial TF influence the dTF respon ses of the population?
        PlotResponses_AccelStim_1(ProjectData,SelectedResponses,AM_Param,AM_Speed,Units_Sel,0);
    case 2
        % What do the single unit responses look like?
        PlotResponses_AccelStim_2(ProjectData,SelectedResponses,AM_Param,AM_Speed,Units_Sel,0);
    case 3
        % Are the single unit reponses influenced by the initial TF?
end






    