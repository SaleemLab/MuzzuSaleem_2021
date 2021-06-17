%% Tomaso Muzzu - 09/03/2020 - UCL - Figure 2D


%% Functions to plot Figure 2 - direction tuning
if ~exist('ProjectData','var')
    [ProjectData AM_UnitResponses AM_Param AM_Speed AM_UOI SelectedResponses AM_UnitResponses_smooth AM_EyeTracking] = LoadDataALL;
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
%%
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
        else
            BL = nanmean(SingleUnitResponse_oneDir(:,1:trialSide_samples),2);
            UnitResp(i) = nanmean(nanmean(SingleUnitResponse_oneDir(:,trialSide_samples+1:end)-BL));
        end
    end
    UnitResp(isnan(UnitResp)) = 0;
    %UnitResp = UnitResp+abs(min(UnitResp));
    GratingDirRespV(k,:) = UnitResp;
end

%% check whether the unit is negatively modulated by perturbation
% load info regarding modulation index of the visual stimulus
% load('DM_visresp_shuffled.mat') 
% % FR_UOI_VIS(1,p_i,:) % mean FR of first 4 seconds of stim
% % FR_UOI_VIS(2,p_i,:) % mean FR of 1s before stim onset
% % FR_UOI_VIS(3,p_i,:) % DM
% % FR_UOI_VIS(4,p_i,:) % MI
% % FR_UOI_VIS(5,p_i,:) % avg baseline FR
% % select best units
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
    clear Ori_vec 
    for j = 1:size(GratingDirResp_trial{i},2)
        Ori_vec(j,:) =  [real(GratingDirResp_trial{i}(:,j)'*exp(2i*degtorad(GratingDirections))) ...
            imag(GratingDirResp_trial{i}(:,j)'*exp(2i*degtorad(GratingDirections)))];
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
clear GratingDirResp GratingDirResp_BL
for k = 1:size(UnitResponses,2)
    clear UnitResp SingleUnitResponse_oneDir
    for i = 1:length(GratingDirections)
        ResponseTime = round(([nanmin(ResponseTimes(:,k,1)) nanmax(ResponseTimes(:,k,2))])*BonvisionFR)+trialSide_samples;
        SingleUnitResponse_oneDir = squeeze(TrialsResp_OI{i}(~isnan(TrialsResp_OI{i}(:,k,1)),k,ResponseTime(1):ResponseTime(2)));
        Baseline = squeeze(TrialsResp_OI{i}(~isnan(TrialsResp_OI{i}(:,k,1)),k,1:trialSide_samples));
        %SingleUnitResponse_oneDir = SingleUnitResponse_oneDir./max(SingleUnitResponse_oneDir);
        if size(SingleUnitResponse_oneDir,2)==1
            BL = nanmean(Baseline); 
            %BL = nanmean(squeeze(TrialsResp_OI{i}(~isnan(TrialsResp_OI{i}(:,k,1)),k,ResponseTime(1)-trialSide_samples:ResponseTime(1))));
        else
            BL = nanmean(Baseline,2);
            %BL = nanmean(squeeze(TrialsResp_OI{i}(~isnan(TrialsResp_OI{i}(:,k,1)),k,ResponseTime(1)-trialSide_samples:ResponseTime(1))),2);
        end
        UnitResp(i) = nanmean(nanmean(SingleUnitResponse_oneDir,2))-nanmean(BL);
        UnitResp_BL(i) = nanmean(BL);
    end
    UnitResp(isnan(UnitResp)) = 0;
    %UnitResp = UnitResp+abs(min(UnitResp));
    GratingDirResp(k,:) = UnitResp;
    GratingDirResp_BL(k,:) = UnitResp_BL;
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

%% ORIENTATION
figure
BinningDir = -45:15:180;
p_thres = 0.01; 
%set(gcf,'Position',[100 100 650 600])
subplot(4,4,[9 10 13 14])
title('Ori. tuned units (vis.)')
%plot(TuningPropsV.P_ori,TuningProps.P_ori,'.k')
%load('AUC_shuffled.mat')
%Sh_responses = AUC_shuffled(:,2:end);
p_pert_th = prctile(Sh_responses(:),99);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PertResp_units = (AUC_shuffled(:,1)>p_pert_th);
%PertResp_units = UnitsPertPos;
PertResp_units = PertRespUnits_pos;
I_dx_pos = find(PertRespUnits_pos); 
PertResp_units(I_dx_pos(12))=0; PertResp_units(I_dx_pos(19))=0;
PertResp_units_V = PertResp_units & TuningPropsV.p_HT2<p_thres;
%PertResp_units = PertResp_units_V;
%PertResp_units = 1:189; PertResp_units = logical(PertResp_units');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on
% plot(TuningProps.P_ori(PertResp_units),TuningPropsV.P_ori(PertResp_units),'or','MarkerSize',5)
plot(TuningProps.P_ori(PertResp_units_V),TuningPropsV.P_ori(PertResp_units_V),'or','MarkerSize',5,'markerfacecolor','r')
xlim([-45 180])
ylim([-45 180])
% plot higher angles for vis. stim. responses
%plot(TuningPropsV.P_ori(TuningPropsV.P_ori>135)-180,TuningProps.P_ori(TuningPropsV.P_ori>135),'.k')
% plot(TuningProps.P_ori(PertResp_units & (TuningPropsV.P_ori>135)),TuningPropsV.P_ori(PertResp_units & (TuningPropsV.P_ori>135))-180,'or','MarkerSize',5)
plot(TuningProps.P_ori(PertResp_units_V & (TuningPropsV.P_ori>135)),TuningPropsV.P_ori(PertResp_units_V & (TuningPropsV.P_ori>135))-180,'or','MarkerSize',5,'markerfacecolor','r')
% plot higher angles for pert. responses
%plot(TuningPropsV.P_ori(TuningProps.P_ori>135),TuningProps.P_ori(TuningProps.P_ori>135)-180,'.k')
% plot(TuningProps.P_ori(PertResp_units & (TuningProps.P_ori>135))-180,TuningPropsV.P_ori(PertResp_units & (TuningProps.P_ori>135)),'or','MarkerSize',5)
plot(TuningProps.P_ori(PertResp_units_V & (TuningProps.P_ori>135))-180,TuningPropsV.P_ori(PertResp_units_V & (TuningProps.P_ori>135)),'or','MarkerSize',5,'markerfacecolor','r')
% bottom left corner
%plot(TuningPropsV.P_ori(TuningProps.P_ori>135 & TuningPropsV.P_ori>135)-180,TuningProps.P_ori(TuningProps.P_ori>135 & TuningPropsV.P_ori>135)-180,'.k')
% plot(TuningProps.P_ori(PertResp_units & (TuningProps.P_ori>135) & TuningPropsV.P_ori>135)-180,TuningPropsV.P_ori(PertResp_units & (TuningProps.P_ori>135)& TuningPropsV.P_ori>135)-180,'or','MarkerSize',5)
plot(TuningProps.P_ori(PertResp_units_V & (TuningProps.P_ori>135) & TuningPropsV.P_ori>135)-180,TuningPropsV.P_ori(PertResp_units_V & (TuningProps.P_ori>135)& TuningPropsV.P_ori>135)-180,'or','MarkerSize',5,'markerfacecolor','r')
% scatterDiagHist(TuningPropsV.P_ori(PertResp_units & TuningPropsV.p_HT2<=p_thres),...
%      TuningProps.P_ori(PertResp_units & TuningPropsV.p_HT2<=p_thres)) 
ylabel('Visual preferred orientation'); xlabel('Perturbation max resp. orientation'); 
% legend({['Ori. tuned vis. units n=' num2str(length(TuningPropsV.P_ori(TuningPropsV.p_HT2<=p_thres)))], ...
%         ['Pert. responsive units n=' num2str(length(TuningPropsV.P_ori(PertResp_units & TuningPropsV.p_HT2<=p_thres)))]})
set(gca,'YAxisLocation', 'left','XTick',-45:45:180,'XTickLabel',[135 0:45:180],'YTick',-45:45:180,'YTickLabel',[135 0:45:180],'TickDir','out'); box off;

subplot(4,4,[1 2 5 6])
%histogram([TuningProps.P_ori(TuningProps.P_ori>135)-180; TuningProps.P_ori],-45:11.25:180,'EdgeColor','none','FaceColor','k')
hold on
% histogram([TuningProps.P_ori(PertResp_units & TuningProps.P_ori>135)-180; TuningProps.P_ori(PertResp_units)], BinningDir,'EdgeColor','none','FaceColor',[229 44 37]/255)
histogram([TuningProps.P_ori(PertResp_units_V & TuningProps.P_ori>135)-180; TuningProps.P_ori(PertResp_units_V)], BinningDir,'EdgeColor','k','FaceColor',[229 44 37]/255)
set(gca,'YAxisLocation','left','XTick',0:45:180,'XTickLabel',[],'TickDir','out'); box off;
%set(gca, 'XDir','reverse')
ylabel('Units')
xlim([-45 180])

subplot(4,4,[11 12 15 16])
%histogram([TuningPropsV.P_ori(TuningPropsV.P_ori>135)-180; TuningPropsV.P_ori],-45:11.25:180,'EdgeColor','none','FaceColor','k')
hold on
% histogram([TuningPropsV.P_ori(PertResp_units & TuningPropsV.P_ori>135)-180; TuningPropsV.P_ori(PertResp_units)], BinningDir,'EdgeColor','none','FaceColor',[229 44 37]/255)
histogram([TuningPropsV.P_ori(PertResp_units_V & TuningPropsV.P_ori>135)-180; TuningPropsV.P_ori(PertResp_units_V)], BinningDir,'EdgeColor','k','FaceColor',[229 44 37]/255)
xlim([-45 180])
set(gca,'XTick',-45:45:180,'XTickLabel',[],'TickDir','out','XAxisLocation', 'bottom'); box off;
ylabel('Units');
set(gca, 'XDir','reverse')
view(90,90)
hold off

% chance level
Distr1 = TuningPropsV.P_ori(PertResp_units_V); % preferred direction for grating
Distr2 = TuningProps.P_ori(PertResp_units_V);  % preferred direction for perturbation
clear Distr_uni N_sh
for sh = 1:1000
    Distr1_sh = Distr1(randperm(length(Distr1)));
    Distr2_sh = Distr1(randperm(length(Distr2)));
    Distr_uni(sh,:) = rad2deg(angdiff(deg2rad(Distr1_sh*2),deg2rad(Distr2_sh*2)))/2;
    %Distr_uni_(sh,:) = histcounts(Distr1_sh-Distr2_sh,-90:18:90);
    %[N_sh(:,:,sh), Xedges, Yedges] = histcounts2(Distr1_sh,Distr2,0:10:180,0:10:180);
end
subplot(4,4,[3 4 7 8])
hold off
Ori_diff = rad2deg(angdiff(deg2rad(TuningPropsV.P_ori(PertResp_units_V)*2),deg2rad(TuningProps.P_ori(PertResp_units_V)*2)))/2;
histogram(Distr_uni,-90:18:90,'EdgeColor','k','FaceColor',[0.5 0.5 0.5],'Normalization','probability')
hold on
histogram(Ori_diff,-90:18:90,'EdgeColor','r','FaceColor',[229 44 37]/255,'Normalization','probability')
% pd = fitdist(TuningPropsV.P_ori(PertResp_units)-TuningProps.P_ori(PertResp_units),'Normal')
% hold on
% x_values = -180:1.8:180;
% y_values = (1/(pd.sigma*sqrt(2*pi))).*exp((-(x_values-pd.mean).^2)./(2*pd.sigma^2));
% plot(x_values,y_values,'LineWidth',2)
xlim([-90 90])
set(gca,'xTick',-90:45:90,'TickDir','out','View',[45,90])
box off;
xlabel('Orientation'); ylabel('Units')
set(gcf,'Position',[10 10 800 700])

suptitle(['n=' num2str(length(TuningProps.P_ori(PertResp_units_V)))])
% [p h] = signrank(Ori_diff)

sum(abs(Ori_diff)<22.5)


% figure
% histogram(Distr_uni,-90:18:90,'EdgeColor','k','FaceColor',[0.5 0.5 0.5],'Normalization','probability')
% hold on
% histogram(Ori_diff,-90:18:90,'EdgeColor','k','FaceColor',[229 44 37]/255,'Normalization','probability')
% xlim([-180 180])
% set(gca,'xTick',-90:45:90,'TickDir','out','View',[45,90])
% box off;
% xlabel('Orientation'); ylabel('Units')
% set(gcf,'Position',[10 10 800 700])
% title(['Delta ori. angle, n=' num2str(sum(PertResp_units_V))])


[h p] = kstest2(histcounts(Distr_uni,-90:5:90)/length(Distr_uni(:)) , histcounts(Ori_diff,-90:5:90)/length(Distr1-Distr2))

p_1 = cumsum(histcounts(Distr_uni,-180:36:180))/(length(Distr_uni(:)));
p_2 = cumsum(histcounts(Distr1-Distr2,0:36:180))/length(Distr1-Distr2);
figure
plot(p_1,'k')
hold on
plot(p_2,'r')

% make 2D histogram out of the scatter plot
[N, Xedges, Yedges] = histcounts2(Distr1,Distr2,0:10:180,0:10:180);
scal_f = 10; sigma = 3;
I = N/sum(N(:));
clear J
J = imresize(I,scal_f);
J = imgaussfilt(J,sigma);
%J = (J/max(J(:)));
I_sh = mean(N_sh,3)/sum(sum(mean(N_sh,3)));
J_sh = imresize(I_sh,scal_f);
J_sh = imgaussfilt(J_sh,sigma);
%J_sh = (J_sh/max(J_sh(:)));
figure
imagesc(J_sh)
set(gca, 'YDir','normal','TickDir','out','box','off')
ylabel('Grating Pref Ori.');
xlabel('Perturbation Pref Ori.');
title('Shuffled')
title('Observed')
title('Observed-shuffled')
% figure
% histfit(TuningPropsV.P_dir(PertResp_units)-TuningProps.P_dir(PertResp_units),12) %-360:45:360)

% figure
% subplot(1,3,1)
% histogram(TuningPropsV.P_ori(PertResp_units_V)-TuningProps.P_ori(PertResp_units_V),-180:36:180,'EdgeColor','k','FaceColor',[229 44 37]/255,'Normalization','probability')
% xlim([-180 180])
% set(gca,'xTick',-180:45:180,'TickDir','out')
% box off;
% xlabel('Orientation'); ylabel('Units')
% set(gcf,'Position',[10 10 800 700])
% title(['Observed, n=' num2str(sum(PertResp_units_V))])
% ylim([0 0.3])
% subplot(1,3,2)
% histogram(Distr_uni(:),-180:36:180,'EdgeColor','k','FaceColor','k','Normalization','probability')
% xlim([-180 180])
% set(gca,'xTick',-180:45:180,'TickDir','out')
% box off;
% xlabel('Orientation'); ylabel('Units')
% set(gcf,'Position',[10 10 800 700])
% title(['Shuffled, n=' num2str(length(Distr_uni(:)))])
% ylim([0 0.3])
% subplot(1,3,3)
% Bars_p(1,:) = histcounts(TuningPropsV.P_ori(PertResp_units_V)-TuningProps.P_ori(PertResp_units_V),-180:36:180)/length(TuningPropsV.P_ori(PertResp_units_V)-TuningProps.P_ori(PertResp_units_V));
% Bars_p(2,:) = histcounts(Distr_uni(:),-180:36:180)/length(Distr_uni(:));
% bar(-180+18:36:180-18,Bars_p(1,:)-Bars_p(2,:),'EdgeColor','k','FaceColor',[229 44 37]/255,'BarWidth', 1)
% xlim([-180 180])
% set(gca,'xTick',-180:45:180,'TickDir','out')
% box off;
% xlabel('Orientation'); ylabel('Units')
% set(gcf,'Position',[10 10 800 700])
% title(['Observed-Shuffled'])


%% DIRECTION
figure
Binning = -90:30:360;
p_thres = 0.01; PertResp_units = (AUC_shuffled(:,1)>0);
%set(gcf,'Position',[100 100 650 600])
subplot(4,4,[9 10 13 14])
title('Dir. tuned units (vis.)')
%plot(TuningPropsV.P_dir,TuningProps.P_dir,'.k')
%load('AUC_shuffled.mat')
%Sh_responses = AUC_shuffled(:,2:end);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PertResp_units = (AUC_shuffled(:,1)>p_pert_th);
%PertResp_units = UnitsPertPos;
PertResp_units = PertRespUnits_pos;
I_dx_pos = find(PertRespUnits_pos); 
PertResp_units(I_dx_pos(12))=0; PertResp_units(I_dx_pos(19))=0;
PertResp_units_V = PertResp_units & TuningPropsV.p_tt<p_thres;
% PertResp_units = PertResp_units_V; 
% PertResp_units = 1:189; PertResp_units = logical(PertResp_units');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on
% plot(TuningProps.P_dir(PertResp_units),TuningPropsV.P_dir(PertResp_units),'^r','MarkerSize',5)
plot(TuningProps.P_dir(PertResp_units_V),TuningPropsV.P_dir(PertResp_units_V),'^r','MarkerSize',5,'markerfacecolor','r')
xlim([-90 360])
ylim([-90 360])
% plot higher angles for vis. stim. responses
%plot(TuningPropsV.P_dir(TuningPropsV.P_dir>270)-360,TuningProps.P_dir(TuningPropsV.P_dir>270),'.k')
% plot(TuningProps.P_dir(PertResp_units & (TuningPropsV.P_dir>270)),TuningPropsV.P_dir(PertResp_units & (TuningPropsV.P_dir>270))-360,'^r','MarkerSize',5)
plot(TuningProps.P_dir(PertResp_units_V & (TuningPropsV.P_dir>270)),TuningPropsV.P_dir(PertResp_units_V & (TuningPropsV.P_dir>270))-360,'^r','MarkerSize',5,'markerfacecolor','r')
% plot higher angles for pert. responses
%plot(TuningPropsV.P_dir(TuningProps.P_dir>270),TuningProps.P_dir(TuningProps.P_dir>270)-360,'.k')
% plot(TuningProps.P_dir(PertResp_units & (TuningProps.P_dir>270))-360,TuningPropsV.P_dir(PertResp_units & (TuningProps.P_dir>270)),'^r','MarkerSize',5)
plot(TuningProps.P_dir(PertResp_units_V & (TuningProps.P_dir>270))-360,TuningPropsV.P_dir(PertResp_units_V & (TuningProps.P_dir>270)),'^r','MarkerSize',5,'markerfacecolor','r')
% bottom left cornere
%plot(TuningPropsV.P_dir(TuningProps.P_dir>270 & TuningPropsV.P_dir>270)-360,TuningProps.P_dir(TuningProps.P_dir>270 & TuningPropsV.P_dir>270)-360,'.k')
% plot(TuningProps.P_dir(PertResp_units & (TuningProps.P_dir>270) & TuningPropsV.P_dir>270)-360,TuningPropsV.P_dir(PertResp_units & (TuningProps.P_dir>270)& TuningPropsV.P_dir>270)-360,'^r','MarkerSize',5)
plot(TuningProps.P_dir(PertResp_units_V & (TuningProps.P_dir>270) & TuningPropsV.P_dir>270)-360,TuningPropsV.P_dir(PertResp_units_V & (TuningProps.P_dir>270)& TuningPropsV.P_dir>270)-360,'^r','MarkerSize',5,'markerfacecolor','r')
% scatterDiagHist(TuningPropsV.P_dir(PertResp_units & TuningPropsV.p_HT2<=p_thres),...
%      TuningProps.P_dir(PertResp_units & TuningPropsV.p_HT2<=p_thres)) 
ylabel('Visual preferred direction'); xlabel('Perturbation max resp. direction'); 
% legend({['Ori. tuned vis. units n=' num2str(length(TuningPropsV.P_dir(TuningPropsV.p_HT2<=p_thres)))], ...
%         ['Pert. responsive units n=' num2str(length(TuningPropsV.P_dir(PertResp_units & TuningPropsV.p_HT2<=p_thres)))]})
set(gca,'YAxisLocation', 'left','XTick',-90:90:360,'XTickLabel',[270 0:90:360],'YTick',-90:90:360,'YTickLabel',[270 0:90:360],'TickDir','out'); box off;

subplot(4,4,[1 2 5 6])
%histogram([TuningProps.P_dir(TuningProps.P_dir>270)-360; TuningProps.P_dir],-90:22.5:360,'EdgeColor','none','FaceColor','k')
hold on
% histogram([TuningProps.P_dir(PertResp_units & TuningProps.P_dir>270)-360; TuningProps.P_dir(PertResp_units)],Binning,'EdgeColor','none','FaceColor',[229 44 37]/255)
histogram([TuningProps.P_dir(PertResp_units_V & TuningProps.P_dir>270)-360; TuningProps.P_dir(PertResp_units_V)],Binning,'EdgeColor','k','FaceColor',[229 44 37]/255)
set(gca,'YAxisLocation', 'left','XTick',0:90:360,'XTickLabel',[],'TickDir','out'); box off;
%set(gca, 'XDir','reverse')
ylabel('Units')
xlim([-90 360])

subplot(4,4,[11 12 15 16])
%histogram([TuningPropsV.P_dir(TuningPropsV.P_dir>270)-360; TuningPropsV.P_dir],-90:22.5:360,'EdgeColor','none','FaceColor','k')
hold on
% histogram([TuningPropsV.P_dir(PertResp_units & TuningPropsV.P_dir>270)-360; TuningPropsV.P_dir(PertResp_units)],Binning,'EdgeColor','none','FaceColor',[229 44 37]/255)
histogram([TuningPropsV.P_dir(PertResp_units_V & TuningPropsV.P_dir>270)-360; TuningPropsV.P_dir(PertResp_units_V)],Binning,'EdgeColor','k','FaceColor',[229 44 37]/255)
xlim([-90 360])
set(gca,'XTick',-90:90:360,'XTickLabel',[],'TickDir','out','XAxisLocation', 'bottom'); box off;
ylabel('Units');
set(gca, 'XDir','reverse')
view(90,90)

% chance level
Distr1 = TuningPropsV.P_dir(PertResp_units_V);
Distr2 = TuningProps.P_dir(PertResp_units_V);
clear Distr_uni
for sh = 1:1000
    Distr1_sh = Distr1(randperm(length(Distr1)));
    Distr2_sh = Distr2(randperm(length(Distr2)));
    Distr_uni(sh,:) = Distr1_sh-Distr2_sh;
end
subplot(4,4,[3 4 7 8])
hold off
Dir_diff = rad2deg(angdiff(deg2rad(TuningPropsV.P_dir(PertResp_units_V)),deg2rad(TuningProps.P_dir(PertResp_units_V))));
histogram(Dir_diff,-180:22.5:180,'EdgeColor','r','FaceColor',[229 44 37]/255,'Normalization','probability')
hold on
histogram(Distr_uni(:),-180:22.5:180,'EdgeColor','k','FaceColor',[0.5 0.5 0.5],'Normalization','probability')
% pd = fitdist(TuningPropsV.P_dir(PertResp_units)-TuningProps.P_dir(PertResp_units),'Normal')
% hold on
% x_values = -360:1.8:360;
% y_values = (1/(pd.sigma*sqrt(2*pi))).*exp((-(x_values-pd.mean).^2)./(2*pd.sigma^2));
% plot(x_values,y_values,'LineWidth',2)
xlim([-180 180])
set(gca,'xTick',-180:45:180,'TickDir','out','View',[45,90])
box off;
xlabel('Direction'); ylabel('Units')
set(gcf,'Position',[10 10 800 700])

suptitle(['n=' num2str(length(TuningProps.P_dir(PertResp_units_V)))])

% chance level
% figure
% histogram(Distr_uni(:),-360:45:360,'EdgeColor','k','FaceColor',[0.5 0.5 0.5],'Normalization','probability')
% hold on
% histogram(TuningPropsV.P_dir(PertResp_units_V)-TuningProps.P_dir(PertResp_units_V),-360:45:360,'EdgeColor','k','FaceColor',[229 44 37]/255,'Normalization','probability')
% xlim([-360 360])
% set(gca,'xTick',-360:90:360,'TickDir','out','View',[45,90])
% box off;
% xlabel('Direction'); ylabel('Units')
% set(gcf,'Position',[10 10 800 700])
% title(['Delta dir. angle, n=' num2str(sum(PertResp_units_V))])

[p h] = kstest2(histcounts(Distr_uni,-360:45:360)/length(Distr_uni(:)) , histcounts(Distr1-Distr2,-360:45:360)/length(Distr1-Distr2))
p_1 = cumsum(histcounts(Distr_uni,-360:45:360))/(length(Distr_uni(:)));
p_2 = cumsum(histcounts(Distr1-Distr2,-360:45:360))/length(Distr1-Distr2);
figure
plot(p_1,'k')
hold on
plot(p_2,'r')


% (Z, smth_win, M, bins1, bins2, Fcircular1, Fcircular2,FnormXonly)
% PertResp_units = 1:189;
% Z(:,1) = histcounts(TuningProps.P_dir(PertResp_units),0:90:360);
% Z(:,2) = histcounts(TuningPropsV.P_dir(PertResp_units),0:90:360);
% smoothhist2D_corrected2_MM(Z,[2 2],[length(0:90:360) length(0:90:360)],0:90:360,0:90:360,1,1,1)
% 
% figure
% histfit(TuningPropsV.P_dir(PertResp_units)-TuningProps.P_dir(PertResp_units),12) %-360:45:360)
% 
% figure
% histfit(TuningPropsV.P_dir(PertResp_units))


%% Max responses
figure
p_thres = 0.01; PertResp_units = (AUC_shuffled(:,1)>0);
%set(gcf,'Position',[100 100 650 600])
subplot(4,4,[9 10 13 14])
title('Dir. tuned units (vis.)')
%plot(TuningPropsV.P_dir,TuningProps.P_dir,'.k')
%load('AUC_shuffled.mat')
%Sh_responses = AUC_shuffled(:,2:end);
p_pert_th = prctile(Sh_responses(:),99);
%PertResp_units = (AUC_shuffled(:,1)>p_pert_th);
PertResp_units = PertRespUnits_pos;
hold on
plot(TuningProps.Pref_Dir(PertResp_units),TuningPropsV.Pref_Dir(PertResp_units),'^r','MarkerSize',5)
xlim([-90 360])
ylim([-90 360])
% plot higher angles for vis. stim. responses
%plot(TuningPropsV.Pref_Dir(TuningPropsV.Pref_Dir>270)-360,TuningProps.Pref_Dir(TuningPropsV.Pref_Dir>270),'.k')
plot(TuningProps.Pref_Dir(PertResp_units & (TuningPropsV.Pref_Dir>270)),TuningPropsV.Pref_Dir(PertResp_units & (TuningPropsV.Pref_Dir>270))-360,'^r','MarkerSize',5)
% plot higher angles for pert. responses
%plot(TuningPropsV.Pref_Dir(TuningProps.Pref_Dir>270),TuningProps.Pref_Dir(TuningProps.Pref_Dir>270)-360,'.k')
plot(TuningProps.Pref_Dir(PertResp_units & (TuningProps.Pref_Dir>270))-360,TuningPropsV.Pref_Dir(PertResp_units & (TuningProps.Pref_Dir>270)),'^r','MarkerSize',5)
% bottom left cornere
%plot(TuningPropsV.Pref_Dir(TuningProps.Pref_Dir>270 & TuningPropsV.Pref_Dir>270)-360,TuningProps.Pref_Dir(TuningProps.Pref_Dir>270 & TuningPropsV.Pref_Dir>270)-360,'.k')
plot(TuningProps.Pref_Dir(PertResp_units & (TuningProps.Pref_Dir>270) & TuningPropsV.Pref_Dir>270)-360,TuningPropsV.Pref_Dir(PertResp_units & (TuningProps.Pref_Dir>270)& TuningPropsV.Pref_Dir>270)-360,'^r','MarkerSize',5)
% scatterDiagHist(TuningPropsV.Pref_Dir(PertResp_units & TuningPropsV.p_HT2<=p_thres),...
%      TuningProps.Pref_Dir(PertResp_units & TuningPropsV.p_HT2<=p_thres)) 
ylabel('Visual preferred orientation'); xlabel('Perturbation max resp. orientation'); 
% legend({['Ori. tuned vis. units n=' num2str(length(TuningPropsV.Pref_Dir(TuningPropsV.p_HT2<=p_thres)))], ...
%         ['Pert. responsive units n=' num2str(length(TuningPropsV.Pref_Dir(PertResp_units & TuningPropsV.p_HT2<=p_thres)))]})
set(gca,'YAxisLocation', 'left','XTick',-90:90:360,'XTickLabel',[270 0:90:360],'YTick',-90:90:360,'YTickLabel',[270 0:90:360],'TickDir','out'); box off;

subplot(4,4,[1 2 5 6])
%histogram([TuningProps.Pref_Dir(TuningProps.Pref_Dir>270)-360; TuningProps.Pref_Dir],-90:22.5:360,'EdgeColor','none','FaceColor','k')
hold on
histogram([TuningProps.Pref_Dir(PertResp_units & TuningProps.Pref_Dir>270)-360; TuningProps.Pref_Dir(PertResp_units)],-90:22.5:360,'EdgeColor','none','FaceColor',[229 44 37]/255)
set(gca,'YAxisLocation', 'left','XTick',0:90:360,'XTickLabel',[],'TickDir','out'); box off;
%set(gca, 'XDir','reverse')
ylabel('Units')
xlim([-90 360])

subplot(4,4,[11 12 15 16])
%histogram([TuningPropsV.Pref_Dir(TuningPropsV.Pref_Dir>270)-360; TuningPropsV.Pref_Dir],-90:22.5:360,'EdgeColor','none','FaceColor','k')
hold on
histogram([TuningPropsV.Pref_Dir(PertResp_units & TuningPropsV.Pref_Dir>270)-360; TuningPropsV.Pref_Dir(PertResp_units)],-90:22.5:360,'EdgeColor','none','FaceColor',[229 44 37]/255)
xlim([-90 360])
set(gca,'XTick',-90:90:360,'XTickLabel',[],'TickDir','out','XAxisLocation', 'bottom'); box off;
ylabel('Units');
set(gca, 'XDir','reverse')
view(90,90)

subplot(4,4,[3 4 7 8])
scatterDiagHist(TuningProps.Pref_Dir(PertResp_units),TuningPropsV.Pref_Dir(PertResp_units))  


%% plot residuals of preferred angles
clear Z
Z(:,1) = TuningProps.P_dir(PertResp_units); % pref direction to perturbation
Z(:,2) = TuningPropsV.P_dir(PertResp_units); % pref direction to grating stim.
smth_win = 3;
M = [100 100];
bins = -15:30:375;
[nF, F, c1, c2, con, H] = smoothhist2D_corrected2_MM(Z, smth_win, M, bins)
[nF, F, c1, c2, con, H] = smoothhist2D_corrected2_MM(Z, smth_win, M, bins1, bins2, Fcircular1, Fcircular2, FnormXonly)

TuningProps.P_ori(PertResp_units); % pref orientation to perturbation
TuningProps.P_ori(PertResp_units); % pref orientation to grating stim.


%% Orientation color map

X_values = [TuningProps.P_ori(PertResp_units_V) ; ...
            TuningProps.P_ori(PertResp_units_V & (TuningPropsV.P_ori>135)) ; ...
            TuningProps.P_ori(PertResp_units_V & (TuningProps.P_ori>135))-180 ; ...
            TuningProps.P_ori(PertResp_units_V & (TuningProps.P_ori>135) & TuningPropsV.P_ori>135)-180];
Y_values = [TuningPropsV.P_ori(PertResp_units_V) ; ...
            TuningPropsV.P_ori(PertResp_units_V & (TuningPropsV.P_ori>135))-180 ; ...
            TuningPropsV.P_ori(PertResp_units_V & (TuningProps.P_ori>135)) ; ...
            TuningPropsV.P_ori(PertResp_units_V & (TuningProps.P_ori>135)& TuningPropsV.P_ori>135)-180];

edges = {-45:1:180 -45:1:180};
[N, c] = hist3( ([X_values, Y_values]) ,'CdataMode','auto','Edges',edges)

scal_f = 1; sigma = 9;
I = N;
J = imresize(I,scal_f);
J = imgaussfilt(J,sigma);
J = (J/max(J(:)));
figure
surf(J,'FaceColor','interp','EdgeColor','none')
view(2)
ylabel('Visual preferred orientation'); xlabel('Perturbation max resp. orientation'); 
xlim([0 length(J)-1])
ylim([0 length(J)-1])
set(gca,'YAxisLocation','left','XTick',0:(length(J)-1)/5:length(J)-1,'XTickLabel',[135 0:45:180],...
        'YTick',0:(length(J)-1)/5:length(J)-1,'YTickLabel',[135 0:45:180],'TickDir','out'); box off;
colorbar

%% Direction color map

X_values = [TuningProps.P_dir(PertResp_units_V) ; ...
            TuningProps.P_dir(PertResp_units_V & (TuningPropsV.P_dir>270)) ; ...
            TuningProps.P_dir(PertResp_units_V & (TuningProps.P_dir>270))-360 ; ...
            TuningProps.P_dir(PertResp_units_V & (TuningProps.P_dir>270) & TuningPropsV.P_dir>270)-360];
Y_values = [TuningPropsV.P_dir(PertResp_units_V) ; ...
            TuningPropsV.P_dir(PertResp_units_V & (TuningPropsV.P_dir>270))-360 ; ...
            TuningPropsV.P_dir(PertResp_units_V & (TuningProps.P_dir>270)) ; ...
            TuningPropsV.P_dir(PertResp_units_V & (TuningProps.P_dir>270)& TuningPropsV.P_dir>270)-360];

edges = {-90:2:360 ; -90:2:360}; 
[N, c] = hist3( ([X_values, Y_values]) ,'CdataMode','auto','Edges',edges)


scal_f = 1; sigma = 18;
I = N;
J = imresize(I,scal_f);
J = imgaussfilt(J,sigma);
J = (J/max(J(:)));
figure
surf(J,'FaceColor','interp','EdgeColor','none')
view(2)
ylabel('Visual preferred direction'); xlabel('Perturbation max resp. direction'); 
xlim([0 length(J)-1])
ylim([0 length(J)-1])
set(gca,'YAxisLocation','left','XTick',0:(length(J)-1)/5:length(J)-1,'XTickLabel',[-270 0:90:360],...
        'YTick',0:(length(J)-1)/5:length(J)-1,'YTickLabel',[-270 0:90:360],'TickDir','out'); box off;
colorbar




ylabel('Visual preferred direction'); xlabel('Perturbation max resp. direction'); 
% legend({['Ori. tuned vis. units n=' num2str(length(TuningPropsV.P_dir(TuningPropsV.p_HT2<=p_thres)))], ...
%         ['Pert. responsive units n=' num2str(length(TuningPropsV.P_dir(PertResp_units & TuningPropsV.p_HT2<=p_thres)))]})
set(gca,'YAxisLocation', 'left','XTick',-90:90:360,'XTickLabel',[270 0:90:360],'YTick',-90:90:360,'YTickLabel',[270 0:90:360],'TickDir','out'); box off;

