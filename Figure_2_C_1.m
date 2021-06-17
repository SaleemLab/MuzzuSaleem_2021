%% Tomaso Muzzu - 09/03/2020 - UCL - Figure 2D


%% Functions to plot Figure 2 - direction tuning
if ~exist('ProjectData','var')
    [ProjectData AM_UnitResponses AM_Param AM_Speed AM_UOI SelectedResponses AM_UnitResponses_smooth] = LoadDataALL;
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
trialSide_seconds = 1;

%%
p_value = 95;
% first 7 animals
if  size(ProjectData,1)>10
    CTRL_exp = 0;
    Animal_1st_idx = [1 5 7 12 15 24 31];
    
    if ~exist('PertResp_units','var')
        % select only perturbation responsive units
        load('AUC_shuffled.mat')
        Sh_responses = AUC_shuffled(:,2:end);
        p_pert_th = prctile(Sh_responses(:),p_value);
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
    p_pert_th = prctile(Sh_responses(:),p_value);
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

%%
% AM_Param :
% conditions = 1 --> nr of recording
% conditions = 2 --> grating direction [0:45:315]
% conditions = 3 --> 1 pert ON, 0 pert OFF
% conditions = 4 --> pert onset
% conditions = 5 --> pert offset

% find angle direction of each trial
Param = AM_Param(:,AM_UOI,2);
UnitResponses = AM_UnitResponses_smooth(:,AM_UOI,:)*60;
% for t = 1:size(UnitResponses,2)
%     OneUnitResponse = squeeze(UnitResponses(:,t,:));
%     UnitResponses_n(:,t,:) = OneUnitResponse./nanmax(OneUnitResponse,[],2);
% end
%UnitResponses = UnitResponses_n;
GratingDirections = unique(Param(~isnan(Param(:,1)),1));
clear Gr_Dir_2D t1 TrialsResp_OI
for i = 1:length(GratingDirections)
   Gr_Dir_2D(:,:,i) = Param==GratingDirections(i); 
   t1 = double(repmat(Gr_Dir_2D(:,:,i), 1, 1, size(UnitResponses,3))); % repeating the selection across time
   t1(t1(:)==0) = nan; 
   TrialsResp_OI{i} = t1.*UnitResponses;
end

%% compute significance of visual responses
% 1) compute MI for visual responses for different angles
ResponseTime = [-1 +4];
BonvisionFR = 60; %Hz
trialSide_samples = 60;
trialSide_seconds = 1; clear FR_vis MI_vis
for i = 1:length(TrialsResp_OI)
    for j = 1:size(TrialsResp_OI{i},2)
        SingleAngleMeanResponse = nanmean(squeeze(TrialsResp_OI{i}(~isnan(TrialsResp_OI{i}(:,j,1)),j,:)));
        FR_vis(j,i,1) = mean(SingleAngleMeanResponse(trialSide_samples:trialSide_samples+ResponseTime(2)*BonvisionFR)) - mean(SingleAngleMeanResponse(1:trialSide_samples));
    end
end
Param = AM_Param(:,AM_UOI,2);
RecordingsNR = unique(AM_Param(1,AM_UOI,1));
for rNR = 1:length(RecordingsNR)
    rNR_i(rNR) = min(find(AM_Param(1,AM_UOI,1)==RecordingsNR(rNR)))
end

% go through responses of single untis
ResponseTime = [-1 +4];
BonvisionFR = 60; %Hz
trialSide_samples = 60;
trialSide_seconds = 1;
clear GratingDirRespV
for k = 1:size(UnitResponses,2)
    clear UnitResp SingleUnitResponse_oneDir
    for i = 1:length(GratingDirections)
        SingleUnitResponse_oneDir = squeeze(TrialsResp_OI{i}(~isnan(TrialsResp_OI{i}(:,k,1)),k,1:trialSide_samples+ResponseTime(2)*BonvisionFR));
        %SingleUnitResponse_oneDir = SingleUnitResponse_oneDir./max(SingleUnitResponse_oneDir);
        if size(SingleUnitResponse_oneDir,2)==1       
            BL = nanmean(SingleUnitResponse_oneDir);
            UnitResp(i) = nanmean(SingleUnitResponse_oneDir(trialSide_samples+1:end))-BL;
            UnitRespSEM(i) = nanstd((SingleUnitResponse_oneDir(trialSide_samples+1:end))-BL)/sqrt(length((SingleUnitResponse_oneDir(trialSide_samples+1:end))-BL));
        else
            BL = nanmean(SingleUnitResponse_oneDir(:,1:trialSide_samples),2);
            UnitResp(i) = nanmean(nanmean(SingleUnitResponse_oneDir(:,trialSide_samples+1:end)-BL));
            UnitRespSEM(i) = nanstd(nanmean(SingleUnitResponse_oneDir(:,trialSide_samples+1:end)-BL))/sqrt(length((SingleUnitResponse_oneDir(:,trialSide_samples+1:end)-BL)));
        end
    end
    UnitResp(isnan(UnitResp)) = 0;
    %UnitResp = UnitResp+abs(min(UnitResp));
    GratingDirRespV(k,:) = UnitResp;
    GratingDirRespV_sem(k,:) = UnitRespSEM;
end

%% check whether the unit is negatively modulated by perturbation
% load info regarding modulation index of the visual stimulus
% load('DM_visresp_shuffled.mat') 
% FR_UOI_VIS(1,p_i,:) % mean FR of first 4 seconds of stim
% FR_UOI_VIS(2,p_i,:) % mean FR of 1s before stim onset
% FR_UOI_VIS(3,p_i,:) % DM
% FR_UOI_VIS(4,p_i,:) % MI
% FR_UOI_VIS(5,p_i,:) % avg baseline FR
% select best units
% SelectedCells = AM_UOI;
% FR_UOI_VIS = FR_UOI_VIS(:,:,SelectedCells);
% DM_vis = squeeze(FR_UOI_VIS(3,1,:));

% % flip the tuning curve if modulation index DM is negative
% for d = 1:size(GratingDirRespV,1)
%     if DM_vis(d)<=0
%        GratingDirRespV(d,:) = -(GratingDirRespV(d,:)-max(GratingDirRespV(d,:)));
%     end
% end

% Compute OSI and DSI for each unit 
clear TuningPropsV Pref_Dir OSI DSI L_ori L_dir
for k = 1:size(UnitResponses,2)
    [R_pref ind] = max(GratingDirRespV(k,:));
    R_null = GratingDirRespV(k,max(mod(ind+4,length(GratingDirections)),1));
    R_ortho1 = GratingDirRespV(k,max(mod(ind+2,length(GratingDirections)),1));
    R_ortho2 = GratingDirRespV(k,max(mod(ind+6,length(GratingDirections)),1));
    Pref_Dir(k,1) = GratingDirections(ind);% perferred direction; 
    % ORI = (R_pref+R_null)-(R_ortho1+R_ortho2)/(R_pref+R_null)+(R_ortho1+R_ortho2)
    OSI(k,1) = ((R_pref+R_null)-(R_ortho1+R_ortho2))/((R_pref+R_null)+(R_ortho1+R_ortho2)); 
    % DSI = (R_pref-R_null)/(R_pref+R_null)
    DSI(k,1) = (R_pref-R_null)/(R_pref+R_null);
    % circular variance or magnitude of the vector L = 1-CirVar, ORI space
    L_ori(k,1) = norm( GratingDirRespV(k,:)*exp(2i*degtorad(GratingDirections))/sum(GratingDirRespV(k,:)) ) ;
    P_ori(k,1) = rad2deg(angle(GratingDirRespV(k,:)*exp(2i*degtorad(GratingDirections))/sum(GratingDirRespV(k,:))))/2;
    if sign(P_ori(k,1))==-1
        P_ori(k,1) = P_ori(k,1)+180;
    end
    % circular variance or magnitude of the vector L = 1-CirVar, DIR space
    L_dir(k,1) = norm( GratingDirRespV(k,:)*exp(1i*degtorad(GratingDirections))/sum(GratingDirRespV(k,:))  );
    P_dir(k,1) = rad2deg(angle(GratingDirRespV(k,:)*exp(1i*degtorad(GratingDirections))/sum(GratingDirRespV(k,:))));
    if sign(P_dir(k,1))==-1
        P_dir(k,1) = 360+P_dir(k,1);
    end 
end
TuningPropsV = table(Pref_Dir,OSI,DSI,L_ori,P_ori,L_dir,P_dir);

%% compute significance of response for Orientation tuning : Hotelling's T-squared test
clear GratingDirResp_trial
for k = 1:size(UnitResponses,2)
    clear UnitResp_all SingleUnitResponse_oneDir
    for i = 1:length(GratingDirections)
        SingleUnitResponse_oneDir = squeeze(TrialsResp_OI{i}(~isnan(TrialsResp_OI{i}(:,k,1)),k,1:trialSide_samples+ResponseTime(2)*BonvisionFR));
        %SingleUnitResponse_oneDir = SingleUnitResponse_oneDir./max(SingleUnitResponse_oneDir);
        if size(SingleUnitResponse_oneDir,2)==1       
            BL = nanmean(SingleUnitResponse_oneDir);
            UnitResp_all{i} = nanmean(SingleUnitResponse_oneDir(trialSide_samples+1:end))-BL;
        else
            BL = nanmean(SingleUnitResponse_oneDir(:,1:trialSide_samples),2);
            UnitResp_all{i} = nanmean(SingleUnitResponse_oneDir(:,trialSide_samples+1:end)-BL,2);
        end
        trials_nr(i) = size(UnitResp_all{i},1);
    end
    clear TrialResp
    for ll = 1:length(GratingDirections)
        TrialResp(ll,:) = UnitResp_all{ll}(1:min(trials_nr));
        TrialResp(ll,isnan(TrialResp(ll,:))) = 0;
        TrialResp(ll,:) = TrialResp(ll,:)+abs(min(TrialResp(ll,:)));
    end
    GratingDirResp_trial{k} = TrialResp;
end
p_thres = 0.01;
for i = 1:length(GratingDirResp_trial)
    clear Ori_vec Ori_vec1
    for j = 1:size(GratingDirResp_trial{i},2)
        Ori_vec(j,:) =  [real(GratingDirResp_trial{i}(:,j)'*exp(2i*degtorad(GratingDirections))) ...
            imag(GratingDirResp_trial{i}(:,j)'*exp(2i*degtorad(GratingDirections)))]/2;
    end
    try
        p_HT2(i) = HotellingT2Test(Ori_vec,p_thres);
    catch
        p_HT2(i) = 1;
    end
    TuningPropsV.p_HT2(i) = p_HT2(i);
end

%% compute significance of response for Direction tuning: direction dot product test
clear p_tt
for k = 1:size(GratingDirResp_trial,2)
    clear OriVec OriAxis Dir_vec ProjMagn
    % step 1 : calculate the average orientation vector onto the orientation
    % axis
    
    OriVec = (GratingDirRespV(k,:)*exp(2i*degtorad(GratingDirections)))/sum(GratingDirRespV(k,:));
    OriAxis = rad2deg(angle(GratingDirRespV(k,:)*exp(2i*degtorad(GratingDirections))))/2;
    if sign(OriAxis)==-1
        OriAxis = OriAxis+180;
    end
    OriAxis = deg2rad(OriAxis);
    % OriAxis = acos(real(OriVec)/imag(OriVec));
    % step 2 : calculate the magnitude of the projection of each direction
    % vector onto the orientation axis
    clear Dir_vec ProjMagn
    for j = 1:size(GratingDirResp_trial{k},2)
        Dir_vec(j,:) =  [real(GratingDirResp_trial{1,k}(:,j)'*exp(1i*degtorad(GratingDirections))) ...
            imag(GratingDirResp_trial{1,k}(:,j)'*exp(1i*degtorad(GratingDirections)))];
        ProjMagn(j) = dot(Dir_vec(j,:),[real(OriVec) imag(OriVec)]);
    end
    % step 3 : compute Student's T-test on the distribution of direction dot
    % products
    [h p_tt(k)] = ttest(ProjMagn); 
    TuningPropsV.p_tt(k) = p_tt(k);
end
p_thres = 0.01;
sum(p_HT2<p_thres)
sum(p_tt<p_thres)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at tuning curves to vis. stim. onset
%% consider only units that perturbation responsive 
Resp_PERT = AUC_shuffled(:,1)>p_pert_th;

%% sort the cells by their response to perturbation (reference for heat map)
Units_Sel = AM_UOI;
[B,I,PertTrial_Ex] = PlotSortingOrder(ProjectData,SelectedResponses,AM_Param,AM_Speed,Units_Sel,0);
% I = descending order taking into account:[pert response; 2s vis stim respo; post-pert reponse]
I = I(~isnan(B(:,1)));
B = reshape(B(~isnan(B(:))),length(B(~isnan(B(:))))/3,3);

%% CIRCULAR PLOTS
%Unit_Inds = find(Combo);
Unit_Inds = find(AM_UOI);
%Unit_Inds = PertRespUnits_pos; % for naive animals (ExperienceEffect.m lines 41-70)
ValidUnits = AM_UOI(SelectedCells);
SP_coords = [95 96 109 110;...
             52 53 66 67; ...
             35 36 49 50; ...
             46 47 60 61; ...
             87 88 101 102; ...
             130 131 144 145; ...
             147 148 161 162; ...
             136 137 150 151];
SP_Polar = [76:79 90:93 104:107 118:121];         

%% tuning of perturbation responses for different directions
% select trials of interest and control trials as well
Trials_PertON = AM_Param(:,:,3)==1; % find indexes where perturbation is on
Trials_PertOFF= AM_Param(:,:,3)==0; % find indexes where perturbation is off
% OSI, DSI, etc. for perturbation responses
% find angle direction of each trial
Param = AM_Param(:,AM_UOI,2);
PertON = Trials_PertON(:,AM_UOI);
UnitResponses = AM_UnitResponses_smooth(:,AM_UOI,:)*60;
% for t = 1:size(UnitResponses,2)
%     OneUnitResponse = squeeze(UnitResponses(:,t,:));
%     UnitResponses_n(:,t,:) = OneUnitResponse./nanmax(OneUnitResponse,[],2);
% end
% UnitResponses = UnitResponses_n;
GratingDirections = unique(Param(~isnan(Param(:,1)),1));
clear Gr_Dir_2D t1 TrialsResp_OI
for i = 1:length(GratingDirections)
   Gr_Dir_2D(:,:,i) = Param==GratingDirections(i); 
   Gr_Dir_2D(:,:,i) = Gr_Dir_2D(:,:,i).*PertON;
   t1 = double(repmat(Gr_Dir_2D(:,:,i), 1, 1, size(UnitResponses,3))); % repeating the selection across time
   t1(t1(:)==0) = nan; 
   TrialsResp_OI{i} = t1.*UnitResponses;
end

% go through responses of single untis
ResponseTimes = AM_Param(:,AM_UOI,4:5);
BonvisionFR = 60; %Hz
trialSide_samples = 60;
trialSide_seconds = 1;
clear GratingDirResp GratingDirResp_BL GratingDirResp_BL1
for k = 1:size(UnitResponses,2)
    clear UnitResp SingleUnitResponse_oneDir
    for i = 1:length(GratingDirections)
        ResponseTime = round(([nanmin(ResponseTimes(:,k,1)) nanmax(ResponseTimes(:,k,2))])*BonvisionFR)+trialSide_samples;
        SingleUnitResponse_oneDir = squeeze(TrialsResp_OI{i}(~isnan(TrialsResp_OI{i}(:,k,1)),k,ResponseTime(1):ResponseTime(2)));
        Baseline = squeeze(TrialsResp_OI{i}(~isnan(TrialsResp_OI{i}(:,k,1)),k,1:trialSide_samples));
        %SingleUnitResponse_oneDir = SingleUnitResponse_oneDir./max(SingleUnitResponse_oneDir);
        if size(SingleUnitResponse_oneDir,2)==1
            BL = nanmean(Baseline); 
            BL1 = nanmean(squeeze(TrialsResp_OI{i}(~isnan(TrialsResp_OI{i}(:,k,1)),k,ResponseTime(1)-trialSide_samples:ResponseTime(1))));
        else
            BL = nanmean(Baseline,2);
            BL1 = nanmean(squeeze(TrialsResp_OI{i}(~isnan(TrialsResp_OI{i}(:,k,1)),k,ResponseTime(1)-trialSide_samples:ResponseTime(1))),2);
        end
        UnitResp(i) = nanmean(nanmean(SingleUnitResponse_oneDir,2))-nanmean(BL);
        UnitResp_sem(i) = (nanstd(nanmean(SingleUnitResponse_oneDir,2))-nanmean(BL))/sqrt(length(nanmean(SingleUnitResponse_oneDir,2))-nanmean(BL));
        UnitResp_BL(i) = nanmean(BL);
        UnitResp_BL1(i) = nanmean(BL1);
    end
    UnitResp(isnan(UnitResp)) = 0;
    %UnitResp = UnitResp+abs(min(UnitResp));
    GratingDirResp(k,:) = UnitResp;
    GratingDirResp_sem(k,:) = UnitResp_sem;
    GratingDirResp_BL(k,:) = UnitResp_BL;
    GratingDirResp_BL1(k,:) = UnitResp_BL1;
end

%% Compute OSI and DSI for each unit            
clear Pref_Dir OSI DSI L_ori P_ori L_dir P_dir MaxResp
for k = 1:size(UnitResponses,2)
    [R_pref ind] = max(GratingDirResp(k,:));
    R_null = GratingDirResp(k,max(mod(ind+4,length(GratingDirections)),1));
    R_ortho1 = GratingDirResp(k,max(mod(ind+2,length(GratingDirections)),1));
    R_ortho2 = GratingDirResp(k,max(mod(ind+6,length(GratingDirections)),1));
    Pref_Dir(k,1) = GratingDirections(ind);% perferred direction; 
    % ORI = (R_pref+R_null)-(R_ortho1+R_ortho2)/(R_pref+R_null)+(R_ortho1+R_ortho2)
    OSI(k,1) = ((R_pref+R_null)-(R_ortho1+R_ortho2))/((R_pref+R_null)+(R_ortho1+R_ortho2)); 
    % DSI = (R_pref-R_null)/(R_pref+R_null)
    DSI(k,1) = (R_pref-R_null)/(R_pref+R_null);
    % circular variance or magnitude of the vector L = 1-CirVar, ORI space
    L_ori(k,1) = norm( GratingDirResp(k,:)*exp(2i*degtorad(GratingDirections))/sum(GratingDirResp(k,:)) ) ;
    P_ori(k,1) = rad2deg(angle(GratingDirResp(k,:)*exp(2i*degtorad(GratingDirections))/sum(GratingDirResp(k,:))))/2;
    if sign(P_ori(k,1))==-1
        P_ori(k,1) = P_ori(k,1)+180;
    end
    % circular variance or magnitude of the vector L = 1-CirVar, DIR space
    L_dir(k,1) = norm( GratingDirResp(k,:)*exp(1i*degtorad(GratingDirections))/sum(GratingDirResp(k,:)) ) ;
    P_dir(k,1) = rad2deg(angle(GratingDirResp(k,:)*exp(1i*degtorad(GratingDirections))/sum(GratingDirResp(k,:))));
    if sign(P_dir(k,1))==-1
        P_dir(k,1) = 360+P_dir(k,1);
    end
end
TuningProps = table(Pref_Dir,OSI,DSI,L_ori,P_ori,L_dir,P_dir);

%% plot units

plotexamples = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Units to plot
% 90 111 240 224
Unit_2_plot = 111;
UI_idx = 1;
while UI_idx<=size(AUC_shuffled,1)
    Target(UI_idx) = find(I==UI_idx) == Unit_2_plot;
    UI_idx = UI_idx + 1;
end
UI_idx = find(Target);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% plot them all
for  UI_idx = 1:size(AUC_shuffled,1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    UI = I(UI_idx) 
%%%%%%%%%%%%%%%%%%%
if DM(UI_idx)
    param_sel = AM_Param(:,SelectedCells,:);
    Pert_lims =[param_sel(min(find(param_sel(:,UI_idx,3)==1)),UI_idx,4) param_sel(min(find(param_sel(:,UI_idx,3)==1)),UI_idx,5)];
    TimeLine = linspace(-1,8.33,size(AM_UnitResponses_smooth,3));
    UnitResponse = squeeze(UnitResponses(:,UI_idx,:));
    Trial_Angle = squeeze(param_sel(:,UI_idx,2));
    Trial_Pert = squeeze(param_sel(:,UI_idx,3));
    figure
    set(gcf,'Position',[100 50 900 675])
    %GratingDirectionsOrdered = [135 90 45 180 0 225 270 315];
    GratingDirectionsOrdered = GratingDirections;
    for i = 1 : length(GratingDirectionsOrdered)
        
        subplot(14,14,SP_coords(i,:))
        % visual
        % plot all trials
        Trial2plot = find(Trial_Angle==GratingDirectionsOrdered(i));
        % plot mean
        [val FirstPart] = min(abs(TimeLine-Pert_lims(1)));
        shadedErrorBar(TimeLine(1:FirstPart-90),...
            mean(UnitResponse(Trial2plot,1:FirstPart-90)),...
            std(UnitResponse(Trial2plot,1:FirstPart-90))/sqrt(size(UnitResponse(Trial2plot,1:FirstPart-90),1)),...
            'lineprop',{'k-','markerfacecolor',[0.5 0.5 0.5]});
        hold on
        xlim([-0.8 max(TimeLine)-0.5]);
        set(gca,'box','off','TickDir','out')
        % perturbation
        Trial2plot_p = find(Trial_Angle==GratingDirectionsOrdered(i) & Trial_Pert==1);
        Trial2plot_np= find(Trial_Angle==GratingDirectionsOrdered(i) & Trial_Pert==0);
        % plot mean
        shadedErrorBar(TimeLine(FirstPart-60:end),...
            mean(UnitResponse(Trial2plot_p,FirstPart-60:end)),...
            std(UnitResponse(Trial2plot_p,FirstPart-60:end))/sqrt(size(UnitResponse(Trial2plot_p,FirstPart-60:end),1)),...
            'lineprop',{'r-','markerfacecolor','r'});
        hold on
        shadedErrorBar(TimeLine(FirstPart-60:end),...
            mean(UnitResponse(Trial2plot_np,FirstPart-60:end)),...
            std(UnitResponse(Trial2plot_np,FirstPart-60:end))/sqrt(size(UnitResponse(Trial2plot_np,FirstPart-60:end),1)),...
            'lineprop',{'k-','markerfacecolor','k'});
        
        Limit_Y_axis(i) = max([ max(mean(UnitResponse(Trial2plot,1:FirstPart-90)))*1.1, ...
            max(mean(UnitResponse(Trial2plot_p,FirstPart-60:end)))*1.1, ...
            max(mean(UnitResponse(Trial2plot_np,FirstPart-60:end)))*1.1]);
        Trial_nr_p(i) = length(Trial2plot_p);
        Trial_nr_np(i) = length(Trial2plot_np);
    end
    for i = 1 : length(GratingDirectionsOrdered)
        %             if i >= 5
        %                 subplot(3,3,i+1)
        %             else
        %                 subplot(3,3,i)
        %             end
        subplot(14,14,SP_coords(i,:))
        ylim([0 max(Limit_Y_axis)]);
        plot([0 0], [0 max(Limit_Y_axis)],'k:'); plot([max(TimeLine)-1 max(TimeLine)-1],[0 max(Limit_Y_axis)],'k:'); % trial start and stop
        plot([Pert_lims(1) Pert_lims(1)], [0,max(Limit_Y_axis)],'r:'); plot([Pert_lims(2) Pert_lims(2)], [0,max(Limit_Y_axis)],'r:'); % pert start and stop
        set(gca,'XTickLabel',[]);
        if i ~=5
            set(gca,'YTickLabel',[]);
            set(gca,'xcolor','none','ycolor','none')
        end
        %legend({['pert. ON n=' num2str((Trial_nr_p(i)))],...
        %   ['pert. OFF n=' num2str((Trial_nr_np(i)))]});
    end
    %subplot(3,3,5)
    %subplot(14,14,SP_Polar);
    hax = axes('Position', [.35, .35, .33, .33]);
    %         bar(hax,y,'EdgeColor','none')
    %         set(hax,'XTick',[])
    polarplot([deg2rad(GratingDirections); 0],max([GratingDirRespV(UI_idx,:) GratingDirRespV(UI_idx,1)],0),'k')
    % hold on
    % polarplot([deg2rad(GratingDirections); 0],[GratingDirResp(UI_idx,:) GratingDirResp(UI_idx,1)],'r')
    hold on
    GR2Plot = GratingDirResp(UI_idx,:);
    GR2Plot(GR2Plot(:)<0) = 0;
    polarplot([deg2rad(GratingDirections); 0],max([GratingDirResp(UI_idx,:) GratingDirResp(UI_idx,1)],0),'r')
    %polarplot([deg2rad(GratingDirections); 0],[GratingDirResp(UI_idx,:)+GratingDirResp_BL(UI_idx,:) GratingDirResp(UI_idx,1)+GratingDirResp_BL(UI_idx,1)],'r')
    hold on
    % plot baseline
    polarplot(linspace(0,2*pi,100), ones(100,1)*mean(GratingDirResp_BL(UI_idx,:)),'k:')
    polarplot(deg2rad(TuningProps.P_ori(UI_idx)),max(GR2Plot),'or')
    polarplot(deg2rad(TuningPropsV.P_ori(UI_idx)),max(GR2Plot),'ok')
    polarplot(deg2rad(TuningProps.P_dir(UI_idx)),max(GR2Plot),'^r')
    polarplot(deg2rad(TuningPropsV.P_dir(UI_idx)),max(GR2Plot),'^k')
    %         if TuningPropsV.p_HT2(UI_idx)<0.05
    %             polarplot(deg2rad(TuningPropsV.P_ori(UI_idx)),max(GR2Plot),'*k')
    %         end
    %         if TuningPropsV.p_tt(UI_idx)<0.05
    %             polarplot(deg2rad(TuningPropsV.P_dir(UI_idx)),max(GR2Plot),'*k')
    %         end
    thetaticks(GratingDirections)
    thetaticklabels({''})
    RadiusTicks = get(gca,'RTick');
    clear RadiusTickLabels;
    for str = 1:length(RadiusTicks)-1
        RadiusTickLabels{str} = '';
    end
    RadiusTickLabels{end+1} = num2str(RadiusTicks(end));
    set(gca,'RTickLabel',RadiusTickLabels)
    suptitle(['Unit ' num2str(UI_idx,3) ', ' ... % num2str(find(I==UI_idx),3)
        'AUC=' num2str(AUC_shuffled(UI_idx,1),2) ', ' ...
        'MI=' num2str(DM(UI_idx),2) ',' ...  %num2str(DM_sh(UI_idx,1),2)
        'Pert.Resp. = ' num2str(AM_UOI(UI_idx))  ', ' ... num2str(Resp_PERT(UI_idx))
        'Ori_a_v=' num2str(TuningPropsV.P_ori(UI_idx),3) ', ' ...
        'Dir_a_v='  num2str(TuningPropsV.P_dir(UI_idx),3) ', ' ...
        'Ori_a_p=' num2str(TuningProps.P_ori(UI_idx),3) ', '...
        'Dir_v_p=' num2str(TuningProps.P_dir(UI_idx),3) ]);
    % % % % % % %% % % % % % %% % % % % % %% % % % % % %
    subplot(14,14, [155 156 157 158 169 170 171 172 183 184 185 186])
    errorbar(GratingDirections, GratingDirRespV(UI_idx,:), GratingDirRespV_sem(UI_idx,:), 'k');
    hold on
    errorbar(GratingDirections, GratingDirResp(UI_idx,:), GratingDirResp_sem(UI_idx,:), 'r');
    ylabel('spikes/s'); 
    ymax = nanmax(max(GratingDirRespV(UI_idx,:))+max(GratingDirRespV_sem(UI_idx,:)), max(GratingDirResp(UI_idx,:))+max(GratingDirResp_sem(UI_idx,:)));
    ymin = min(0,min(min(GratingDirRespV(UI_idx,:))+min(GratingDirRespV_sem(UI_idx,:)), min(GratingDirResp(UI_idx,:))+min(GratingDirResp_sem(UI_idx,:))));
    try
        ylim([ymin ymax]);
    end
    xlim([-45 360]);
    title(['#' num2str(UI_idx) ',P' num2str(PertRespUnits_pos(UI_idx)) ',p_H_T = ' num2str(TuningPropsV.p_HT2(UI_idx)<0.01,3) ',p_t_t = ' num2str(TuningPropsV.p_tt(UI_idx)<0.01,3)])
    set(gca,'XTick',0:45:315,'box','off','TickDir','out')
    
    % % % % % % %% % % % % % %% % % % % % %% % % % % % %
    [v ind_m] = max(GratingDirRespV(UI_idx,:));
    for i = 1 : length(GratingDirectionsOrdered)
        % plot all trials
        Trial2plot = find(Trial_Angle==GratingDirectionsOrdered(i));
        % plot mean
        [val FirstPart] = min(abs(TimeLine-Pert_lims(1)));
        Part_VResponses(i,1) = max(mean(UnitResponse(Trial2plot,trialSide_samples:trialSide_samples+120))); % save max during first 2 seconds
        Part_VResponses(i,2) = mean(mean(UnitResponse(Trial2plot,trialSide_samples+120:trialSide_samples+240))); % save mean during last 2 seconds
    end
    subplot(14, 14, [165 166 167 179 180 181 193 194 195])
    plot(Part_VResponses(:,2),Part_VResponses(:,1),'.k','MarkerSize',20)
    hold on
    plot(Part_VResponses(ind_m,2),Part_VResponses(ind_m,1),'+k','MarkerSize',4)
    hold on
    plot([0 max(Part_VResponses(:))], [0 max(Part_VResponses(:))],':k','LineWidth',2)
    xlim([min(Part_VResponses(:)) max(Part_VResponses(:))]); ylim([min(Part_VResponses(:)) max(Part_VResponses(:))]);
    set(gca,'box','off','TickDir','out')
    xlabel('Mean 2s-4s V.resp.'); ylabel('max 0-2s V.resp.')
    
             saveas(gcf,['E:\Dropbox\UCL\Tomaso-Aman\VisPerturbation\Fig_2\panels\AllUnits' filesep 'Unit_' num2str(UI) '.pdf']);
    %   saveas(gcf,['E:\Dropbox\UCL\Tomaso-Aman\VisPerturbation\Fig_3\panels\Naive\gratingDirExamples' filesep 'Unit_' num2str(UI_idx) '.pdf']);
    %
    %         close all
    %     saveas(gcf,['X:\DATA\PROJECTS\VisPerturbation\Figures\Fig_2\panels\ExampleUnits' filesep FileSuffix '_' num2str(UI) '.fig']);
    %
    %input('');
    close all
end


end

















%% relevant units
Nr = [44 51 52 61 62 63 67 71 76 87 92 129 136 158 160 161 163 185

115


%% plot responses to pert and grating together
% Units to plot
% 90 111 240 224
Unit_2_plot = 111;
UI_idx = 1;
while UI_idx<=830
    Target(UI_idx) = find(I==UI_idx) == Unit_2_plot;
    UI_idx = UI_idx + 1;
end
UI_idx = find(Target);

figure
polarplot([deg2rad(GratingDirections); 0],max([GratingDirResp(UI_idx,:) GratingDirResp(UI_idx,1)],0),'r')
hold on
polarplot([deg2rad(GratingDirections); 0],max([GratingDirRespV(UI_idx,:) GratingDirRespV(UI_idx,1)],0),'k')
figure
plot(GratingDirections,GratingDirResp(UI_idx,:),'r')
hold on
plot(GratingDirections,GratingDirRespV(UI_idx,:),'k')
figure
plot(GratingDirections, GratingDirResp(UI_idx,:)-GratingDirRespV(UI_idx,:),'k')
hold on
plot(GratingDirections, GratingDirResp(UI_idx,:)-GratingDirResp_BL1(UI_idx,:),'r')

[R P] = corrcoef(GratingDirResp(UI_idx,:)-mean(GratingDirResp(UI_idx,:)), ...
GratingDirRespV(UI_idx,:)-mean(GratingDirRespV(UI_idx,:)))

[h , p ] = chi2gof(GratingDirResp(UI_idx,:))

%% Pearson Correlation coefficient for all perturbation resp. units
PertRespUnits_i = find(PertRespUnits_pos);
for i = 1:length(PertRespUnits_i)
    
    UI_idx = PertRespUnits_i(i);
    
    [r p] = corrcoef(GratingDirResp(UI_idx,:)-mean(GratingDirResp(UI_idx,:)), ...
                    GratingDirRespV(UI_idx,:)-mean(GratingDirRespV(UI_idx,:)));
    R(i,1:2) = [r(1,2) p(1,2)]; 
    
end

figure
for i = 1:length(R)
    hold on
    if R(i,2)<=0.05
        plot(R(i,2),R(i,1),'r.')
    else
        plot(R(i,2),R(i,1),'k.')         
    end
end
min(abs(R(find(abs(R(:,2))<=0.05),1)))
% cumulative distribution of R
figure
plot(0:0.05:1,[0 cumsum(histcounts(abs(R(:,2)),0:0.05:1))]./sum(PertResp_units)*100 );
hold on
Min_th = min(abs(R(find(abs(R(:,2))<=0.05),1)));
plot([Min_th Min_th],[0 100], ':b')
xlabel('R'); ylabel('Pert. resp. units')
box off
set(gca,'TickDir','out','XTick',0:0.2:1,'YTick',0:20:100)

%% look at responses in angles where there is no visual response
clear X Y
figure
for i = 1:length(PertRespUnits_i)
    
    UI_idx = PertRespUnits_i(i);
    
    plot(GratingDirRespV(UI_idx,:), GratingDirResp(UI_idx,:), 'r.');
    hold on
    if i == 1
        X = GratingDirRespV(UI_idx,:);
        Y = GratingDirResp(UI_idx,:);
    else
        X = [X GratingDirRespV(UI_idx,:)];
        Y = [Y GratingDirResp(UI_idx,:)];
    end
    
end
plot([-100 100],[-100 100],':')
xlim([-10 40])
ylim([-10 40])
set(gca, 'box', 'off', 'TickDir','out', 'XTick',[-10:10:40], 'YTick',[-10:10:40])
xlabel('Drifting Grating (spikes/s)'); ylabel('Perturbation (spikes/s)')

dlm = fitlm(X,Y)
m_coeff = dlm.Coefficients.Estimate;
x_values = -10:1:40;
plot(x_values,m_coeff(2)+m_coeff(1)*x_values)


