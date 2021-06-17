%% Tomaso Muzzu - UCL - 26/03/2020

% Acceleration stimulus analysis

%% load data
if ~exist('ProjectData','var')
    [ProjectData AM_UnitResponses AM_Param AM_Speed AM_UOI SelectedResponses AM_UnitResponses_smooth] = LoadDataALL;
end

%% options for control experiments in naive mice
% rec_1day = [1 4 7];
% rec_2day = [2 5 8];
% rec_3day = [3 6 9];
% for i = 1:length(rec_1day)
%     units_1day(:,i) = (AM_Param(1,:,1) == rec_1day(i));
%     units_2day(:,i) = (AM_Param(1,:,1) == rec_2day(i));
%     units_3day(:,i) = (AM_Param(1,:,1) == rec_3day(i));
% end
% AM_UOI(:,2) = sum(units_1day,2);
% AM_UOI(:,3) = sum(units_2day,2);
% AM_UOI(:,4) = sum(units_3day,2);

%% first 7 animals
if  size(ProjectData,1)>10
    CTRL_exp = 0;
    Animal_1st_idx = [1 5 7 12 15 24 31];   
    if ~exist('PertResp_units','var')
        % select only perturbation responsive units
        load('AUC_shuffled.mat')
        Sh_responses = AUC_shuffled(:,2:end);
        p_pert_th = prctile(Sh_responses(:),99);
        PertResp_units = (AUC_shuffled(:,1)>p_pert_th);
        % select only pos. modulated perturbation responsive units
        load('DM_pert_shuffled.mat')
        DM = DM_sh(:,1);
        DM_sign_i(:,1) = DM>0;
        DM_sign_i(:,2) = DM<=0;
        % select only pos. modulated perturbation responsive units
        PertRespUnits_pos = PertResp_units & DM_sign_i(:,1);
        PertRespUnits_neg = PertResp_units & DM_sign_i(:,2);
    end
    
else
     CTRL_exp = 1;
    % naive animals
    Animal_1st_idx = [1 4 7];
    % select only perturbation responsive units
    load('AUC_shuffled_CTRL_1.mat')
    Sh_responses = AUC_shuffled(:,2:end);
    p_pert_th = prctile(Sh_responses(:),99);
    PertResp_units = (AUC_shuffled(:,1)>p_pert_th);
    % select only pos. modulated perturbation responsive units
    load('DM_pert_shuffled_CTRL.mat')
    DM = DM_sh(:,1);
    DM_sign_i(:,1) = DM>0;
    DM_sign_i(:,2) = DM<=0;
    % select only pos. modulated perturbation responsive units
    PertRespUnits_pos = PertResp_units & DM_sign_i(:,1);
    PertRespUnits_neg = PertResp_units & DM_sign_i(:,2);
end


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
Cond = 1; % %% all units together, distinguishing for accel. values and starting TF
%Cond = 2; % %% plot single units, distinguishing for accel. values
%Cond = 3; % %% plot single units, distinguishing for accel. values and starting TF
SelectedResponses = AccelResponsesSelection(ProjectData,AM_UnitResponses_smooth,AM_Param,AM_Speed(:,:,1:end-1),Units_Sel,Cond)%,'run',Run_TH);
Recs = unique(AM_Param(1,:,1)); clear UnitPerRec
for k = 1:length(Recs)
    UnitPerRec(k) = max(find(AM_Param(1,:,1)==Recs(k)));
end
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

%% plot only the perturbation responsive units

% have info about the direction and orientation selectivity

% plot the 12 units identified somewhere else

% 









