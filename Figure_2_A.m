%% Tomaso Muzzu - UCL - 06/03/2020

%% Functions to plot Figure 2A - direction tuning
if ~exist('ProjectData','var')
    [ProjectData AM_UnitResponses AM_Param AM_Speed AM_UOI SelectedResponses AM_UnitResponses_smooth] = LoadDataALL;
end

%%
p_value = 95;

%% first 7 animals
thres = 95;
if  size(ProjectData,1)==37
    CTRL_exp = 0;
    Animal_1st_idx = [1 5 7 12 15 24 31];   
    if ~exist('PertResp_units','var')
        % select only perturbation responsive units
        load('AUC_shuffled.mat')
        Sh_responses = AUC_shuffled(:,2:end);
        p_pert_th = prctile(Sh_responses(:),thres);
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
elseif size(ProjectData,1)==10
    CTRL_exp = 1;
    % naive animals
    Animal_1st_idx = [1 4 7];
    % select only perturbation responsive units
    load('AUC_shuffled_CTRL_1.mat')
    Sh_responses = AUC_shuffled(:,2:end);
    p_pert_th = prctile(Sh_responses(:),thres);
    PertResp_units = (AUC_shuffled(:,1)>p_pert_th);
    % select only pos. modulated perturbation responsive units
    load('DM_pert_shuffled_CTRL.mat')
    DM = DM_sh(:,1);
    DM_sign_i(:,1) = DM>0;
    DM_sign_i(:,2) = DM<=0;
    % select only pos. modulated perturbation responsive units
    PertRespUnits_pos = PertResp_units & DM_sign_i(:,1);
    PertRespUnits_neg = PertResp_units & DM_sign_i(:,2);
else
    % select only perturbation responsive units
    load('AUC_shuffled_CTRL_1.mat')
    AUC_shuffled_CTRL = AUC_shuffled;
    AUC_shuffled_pp_CTRL = AUC_shuffled_pp;
    load('AUC_shuffled.mat')
    AUC_shuffled = cat(1,AUC_shuffled,AUC_shuffled_CTRL);
    AUC_shuffled_pp = cat(1,AUC_shuffled_pp,AUC_shuffled_pp_CTRL);
    Sh_responses = AUC_shuffled(:,2:end);
    p_pert_th = prctile(Sh_responses(:),thres);
    PertResp_units = (AUC_shuffled(:,1)>p_pert_th);
    % select only pos. modulated perturbation responsive units
    load('DM_ALL.mat')
    DM_sign_i(:,1) = DM>0;
    DM_sign_i(:,2) = DM<=0;
    % select only pos. modulated perturbation responsive units
    PertRespUnits_pos = PertResp_units & DM_sign_i(:,1);
    PertRespUnits_neg = PertResp_units & DM_sign_i(:,2);
end



%% direction
if size(AM_UOI,2)==1
    SelectedCells = AM_UOI;
else
    SelectedCells = AM_UOI(:,1) & AM_UOI(:,2);
end
% select trials of interest and control trials as well
Trials_PertON = AM_Param(:,:,3)==1; % find indexes where perturbation is on
Trials_PertOFF= AM_Param(:,:,3)==0; % find indexes where perturbation is off

BonvisionFR = 60; %Hz
trialSide_samples = 60;

%% FIGURE 2A 
%%%% with components congruent with running direction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AngleCombo = [0 ; 180]; 
AngleCombo = [45 ; 225];
AngleCombo = [90 ; 270];
AngleCombo = [135 ; 315];

AngleCombo = [0 180; 90 270];
AngleCombo = [0 45 315; 135 180 225];

clear SelectedResponses
for a = 1:size(AngleCombo,1)
    Angles_OI = AngleCombo(a,:);
    clear AngleFilter
    for i = 1:length(Angles_OI)
        % find all PertON trials with the grating direction of interest
        AngleFilter_t(:,:,i) =  AM_Param(:,:,2)==Angles_OI(i);
        if i == length(Angles_OI)
            AngleFilter = sum(AngleFilter_t,3);
        end
    end
    clear AM_AngleResponse AM_AngleResponseControl
    AM_AngleResponse = AngleFilter & Trials_PertON;
    % find all PertOFF trials with the grating direction of interest
    AM_AngleResponseControl = AngleFilter & Trials_PertOFF;  
    
    
    % generate 3D matrix for including time responses
    t1 = double(repmat(AM_AngleResponse, 1, 1, size(AM_UnitResponses_smooth,3))); % repeating the selection across time
    t2 =  double(repmat(AM_AngleResponseControl, 1, 1, size(AM_UnitResponses_smooth,3))); % repeating the selection across time
    % nan all the zero values
    t1(t1(:)==0) = nan; t2(t2(:)==0) = nan;
    % select only the responses during selected trials (pert ON and OFF separately)
    TrialsResp_OI = t1.*AM_UnitResponses_smooth;
    TrialsRespControl_OI = t2.*AM_UnitResponses_smooth;
    % compute the mean responses for each unit during trials with pert ON
    % and OFF
    TrialsResp_OI_2D = squeeze(nanmean(TrialsResp_OI,1));
    TrialsRespControl_OI_2D = squeeze(nanmean(TrialsRespControl_OI,1));
    TrialsResp_max = max(nanmax(TrialsResp_OI_2D,[],2),nanmax(TrialsRespControl_OI_2D,[],2));
    Response = TrialsResp_OI_2D./TrialsResp_max;
    ResponseControl = TrialsRespControl_OI_2D./TrialsResp_max;
    Response_OI = Response(SelectedCells,:);
    ResponseControl_OI = ResponseControl(SelectedCells,:);
    
    SelectedResponses{a,1} = Response_OI;
    SelectedResponses{a,2} = ResponseControl_OI;
    
end

%PlotResponses(ProjectData,{SelectedResponses{1,1}(UnitsPertPos,:) SelectedResponses{2,1}(UnitsPertPos,:)},AM_Param,AM_Speed,UnitsPertPos,0); 
%PlotResponses(ProjectData,{SelectedResponses{1,1}(UnitsPertNeg,:) SelectedResponses{2,1}(UnitsPertNeg,:)},AM_Param,AM_Speed,UnitsPertNeg,0); 

PlotMeanCurveAngle(ProjectData,{SelectedResponses{1,1}(PertRespUnits_pos,:) SelectedResponses{2,1}(PertRespUnits_pos,:)},AM_Param,AM_Speed,PertRespUnits_pos,AngleCombo)
PlotMeanCurveAngle(ProjectData,{SelectedResponses{1,1}(PertRespUnits_neg,:) SelectedResponses{2,1}(PertRespUnits_neg,:)},AM_Param,AM_Speed,PertRespUnits_neg,AngleCombo)


%% evaluate the statistically significant difference of the responses for the various angles as a function of time
AllAnglesCombo = nchoosek(0:45:315,2);

AngleCombo = [0 180; 90 270];
AngleCombo = [0 45 315; 135 180 225];

clear SelectedResponses
for Combo = 1:length(AllAnglesCombo)
    AngleCombo = AllAnglesCombo(Combo,:);
    for a = 1:size(AngleCombo,1)
        Angles_OI = AngleCombo(a,:);
        clear AngleFilter
        for i = 1:length(Angles_OI)
            % find all PertON trials with the grating direction of interest
            AngleFilter_t(:,:,i) =  AM_Param(:,:,2)==Angles_OI(i);
            if i == length(Angles_OI)
                AngleFilter = sum(AngleFilter_t,3);
            end
        end
        clear AM_AngleResponse AM_AngleResponseControl
        AM_AngleResponse = AngleFilter & Trials_PertON;
        % find all PertOFF trials with the grating direction of interest
        AM_AngleResponseControl = AngleFilter & Trials_PertOFF;
        
        
        % generate 3D matrix for including time responses
        t1 = double(repmat(AM_AngleResponse, 1, 1, size(AM_UnitResponses_smooth,3))); % repeating the selection across time
        t2 =  double(repmat(AM_AngleResponseControl, 1, 1, size(AM_UnitResponses_smooth,3))); % repeating the selection across time
        % nan all the zero values
        t1(t1(:)==0) = nan; t2(t2(:)==0) = nan;
        % select only the responses during selected trials (pert ON and OFF separately)
        TrialsResp_OI = t1.*AM_UnitResponses_smooth;
        TrialsRespControl_OI = t2.*AM_UnitResponses_smooth;
        % compute the mean responses for each unit during trials with pert ON
        % and OFF
        TrialsResp_OI_2D = squeeze(nanmean(TrialsResp_OI,1));
        TrialsRespControl_OI_2D = squeeze(nanmean(TrialsRespControl_OI,1));
        TrialsResp_max = max(nanmax(TrialsResp_OI_2D,[],2),nanmax(TrialsRespControl_OI_2D,[],2));
        Response = TrialsResp_OI_2D./TrialsResp_max;
        ResponseControl = TrialsRespControl_OI_2D./TrialsResp_max;
        Response_OI = Response(SelectedCells,:);
        ResponseControl_OI = ResponseControl(SelectedCells,:);
        
        SelectedResponses{Combo,1} = Response_OI;
        SelectedResponses{Combo,2} = ResponseControl_OI;
        
    end
    
end

% select only positive units
for i = 1:size(SelectedResponses,1)
    Responses_pos{i,1} = SelectedResponses{i,1}(PertRespUnits_pos,:);
    Responses_pos{i,2} = SelectedResponses{i,2}(PertRespUnits_pos,:);
    Responses_neg{i,1} = SelectedResponses{i,1}(PertRespUnits_neg,:);
    Responses_neg{i,2} = SelectedResponses{i,2}(PertRespUnits_neg,:);
end
% compute the significance tests for all combinations every 3 samples
% (50ms) during the perturbation period
PertTrial_Ex(1) = 1;
PertTrial_Ex(2) = min(find(ProjectData.Trial_data{1,1}.PerturbationON==1))
TF = ProjectData.Trial_data{PertTrial_Ex(1),1}.TF_in_trial{PertTrial_Ex(2),1}(1:end-1);
PertPeriod = [min(find(TF==0)) max(find(TF==0))];
for i = 1:size(Responses_pos,1)
    r = 1;
    for dt = PertPeriod(1):3:PertPeriod(2)
        [h(i,r) p(i,r)] = ttest2(reshape(Responses_pos{i,1}(:,dt:dt+2),[1,prod(size(Responses_pos{i,1}(:,dt:dt+2)))]), ...
                        reshape(Responses_pos{i,2}(:,dt:dt+2),[1,prod(size(Responses_pos{i,1}(:,dt:dt+2)))]) );
        r = r+1;
    end
end
% reorder the array of the p-values and plot them
AllAnglesCombo
CorrMatrix = nan(length(0:45:315),prod(size(p)));
Directions = 0:45:315;
CorrMatrix_h = nan(8,prod(size(h(AllAnglesCombo(:,1) == Directions(1),:))));
CorrMatrix_p = nan(8,prod(size(p(AllAnglesCombo(:,1) == Directions(1),:))));
for i = 1:length(Directions)
    CorrMatrix_h(i,end-prod(size(h(AllAnglesCombo(:,1) == Directions(i),:)))+1:end) = ...
        reshape(h(AllAnglesCombo(:,1) == Directions(i),:),[1,prod(size(h(AllAnglesCombo(:,1) == Directions(i),:)))]);
    CorrMatrix_p(i,end-prod(size(h(AllAnglesCombo(:,1) == Directions(i),:)))+1:end) = ...
        reshape(p(AllAnglesCombo(:,1) == Directions(i),:),[1,prod(size(h(AllAnglesCombo(:,1) == Directions(i),:)))]);
end

colormap(gray)
figure
imagesc(CorrMatrix_p)
colorbar
colormap(flipud(gray))
set(gca,'box','off','YDir','normal','XDir','normal','TickDir','out','XAxisLocation','bottom')
set(gca,'yTickLabel',0:45:315,'xTick',(0:size(h,2):size(CorrMatrix_h,2)-size(h,2))+size(h,2)/2,'XTickLabel',45:45:315)
hold on
g_x=[0:size(h,2):size(CorrMatrix_h,2)]; % user defined grid X [start:spaces:end]
g_y=[0.5:1:8.5]; % user defined grid Y [start:spaces:end]
for i=1:length(g_x)
   plot([g_x(i) g_x(i)],[g_y(1) g_y(end)],'Color',[0.5 0.5 0.5]) %y grid lines
   hold on    
end
for i=1:length(g_y)
   plot([g_x(1) g_x(end)],[g_y(i) g_y(i)],'Color',[0.5 0.5 0.5]) %x grid lines
   hold on    
end
