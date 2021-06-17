%% Functions to plot Figure 2 - acceleration / speed tuning


%% load main stimulus daata
if ~exist('ProjectData','var')
    [ProjectData AM_UnitResponses AM_Param AM_Speed AM_UOI SelectedResponses AM_UnitResponses_smooth_acc] = LoadDataALL;
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
%%
CTRL_exp = 1;
% naive animals
Animal_1st_idx = [1 4 7];
% select only perturbation responsive units
load('AUC_shuffled_CTRL_1.mat')
Sh_responses = AUC_shuffled(:,2:end);
p_pert_th = prctile(Sh_responses(:),p_value);
PertResp_units = (AUC_shuffled(:,1)>p_pert_th);
% select only pos. modulated perturbation responsive units
load('DM_CTRL.mat')
%DM = DM_sh(:,1);
DM_sign_i(:,1) = DM>0;
DM_sign_i(:,2) = DM<=0;
% select only pos. modulated perturbation responsive units
PertRespUnits_pos = PertResp_units & DM_sign_i(:,1);
PertRespUnits_neg = PertResp_units & DM_sign_i(:,2);

%% Load raw and formatted data of acceleration stimulus
[ProjectData_acc AM_UnitResponses_acc AM_Param_acc AM_Speed_acc AM_UOI_acc] = LoadData;

%% Initial selection of best units
% to be run once in theory, repeated in case params need changing
FR_thres = 0.1;
AM_UOI_acc = UnitsSelection(ProjectData_acc,AM_UnitResponses_acc,FR_thres);
% save([ DataFolder filesep 'UnitsResponse.mat'], 'AM_UnitResponses', 'AM_Speed', 'AM_Param', 'AM_UOI','-v7.3')
% ATM, units with mean FR>0.1Hz and comparable FR in first and last third
% of recording are kept

%% smooth data with Gaussian filter
% 3rd argument is width of Gaussian filter
% 4th OPTIONAL argument is width of std - sigma=1/3 width by defeault;
Gauss_width = 0.3; % in seconds
AM_UnitResponses_smooth_acc = Smooth_AM_FR(ProjectData_acc,AM_UnitResponses_acc,Gauss_width);

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

%% Selection of best units of main stimulus applied to acceleration protocol


%%%Cond_Sel = {} ; % {optional('run'), optional(grating dir=[0:45:315])}
Units_Sel = AM_UOI(:,1) ; % use best units selected from main stimulus 
% Units_Sel = logical(ones(size(AM_UOI(:,1),1),1));
Run_TH = 3; % running speed threshold is 1 cm/s to be used like this: ..., 'run', Run_TH)

% Cond = 0; % %% all units together, distinguishing for accel. values
% Cond = 1; % %% all units together, distringuishing for accel. values and starting TF
 Cond = 2; % %% plot single units, distringuishing for accel. values
% Cond = 3; % %% plot single units, distringuishing for accel. values and starting TF

SelectedResponses = AccelResponsesSelection(ProjectData_acc,AM_UnitResponses_smooth_acc,AM_Param_acc,AM_Speed_acc(:,:,1:end-1),Units_Sel,Cond); %,'run',Run_TH);
% Units_Sel= true(449,1);
% Units_Sel(UnitPerRec(7)+1:UnitPerRec(8)) = 0;
% Units_Sel(UnitPerRec(9)+1:UnitPerRec(10)) = 0;

%%
% AM_param_acc
% conditions = 1 --> nr of recording
% conditions = 2 --> TF at stimulus onset
% conditions = 3 --> TF at simutlus offset
% conditions = 4 --> dTF shown
% conditions = 5 --> index of when Tf starts changing
% conditions = 6 --> index of when Tf finishes changing
AM_UOI = Units_Sel;

% find all the TF's at the start of the trials
StarTF_2dmask = squeeze(AM_Param_acc(:,:,2));
start_TF_values = unique(StarTF_2dmask(~isnan(StarTF_2dmask(:))));
% separate the responses for the different acceleration values
dTF_2Dmask = squeeze(AM_Param_acc(:,:,4));
dTF_values = unique(dTF_2Dmask(~isnan(dTF_2Dmask(:))));
% find trials in which dTF lasts 0.5s and 1s
dTF_dur_2D = (squeeze(AM_Param_acc(:,:,6))-squeeze(AM_Param_acc(:,:,5)))/60; % 60 is the sampling frequency
dTF_long_idx = dTF_dur_2D>0.8;
% find three trials for each condition for which there is dTF>0; dTF=0;
% dTF<0
clear TF_10_i TF_05_i I_10 I_05
for i = 1:length(dTF_values)
    TF_10_i(i) = min(find((AM_Param_acc(:,:,4).*dTF_long_idx)==dTF_values(i) & (AM_Param_acc(:,:,2).*dTF_long_idx)==6));
    TF_05_i(i) = min(find((AM_Param_acc(:,:,4).*(~dTF_long_idx))==dTF_values(i) & (AM_Param_acc(:,:,2).*(~dTF_long_idx))==6));
end
%convert to 2D indexes
ReferSize = [size(AM_Param_acc,1),size(AM_Param_acc,2)];
[I_10(:,1),I_10(:,2)] = ind2sub(ReferSize,TF_10_i); % from linear to 2D indexes
[I_05(:,1),I_05(:,2)] = ind2sub(ReferSize,TF_05_i); % from linear to 2D indexes

    
%% Process
% tuning curves for acceleration
% tuning curves for various speeds: measure them after the acceleration
% period : 
AccelValues = -3:1:3;
Accel_Dur = [0.5 1];
SpeedValues_05 = start_TF_values+(AccelValues*Accel_Dur(1));
SpeedValues_10 = start_TF_values+(AccelValues*Accel_Dur(2));

% AM_param_acc
% conditions = 1 --> nr of recording
% conditions = 2 --> TF at stimulus onset
% conditions = 3 --> TF at simutlus offset
% conditions = 4 --> dTF shown
% conditions = 5 --> index of when TF starts changing
% conditions = 6 --> index of when TF finishes changing

GoodUnits = find(AM_UOI);
U_speed = unique([SpeedValues_05 SpeedValues_10]);
AM_UnitResponses_acc_s = AM_UnitResponses_smooth_acc(:,AM_UOI,:);
AM_Param_acc_s = AM_Param_acc(:,AM_UOI,:);
clear Resp_mean Resp_sem
for i = 1:length(U_speed)
    for j = 1:size(AM_UnitResponses_acc_s,2)
        Indexes = AM_Param_acc_s(:,j,3) == U_speed(i);
        Resp_mean(j,i) = mean(nanmean(squeeze(AM_UnitResponses_acc_s(Indexes,j,end-trialSide_samples-60:end-trialSide_samples-1))));
        Resp_sem(j,i) = std(nanmean(squeeze(AM_UnitResponses_acc_s(Indexes,j,end-trialSide_samples-60:end-trialSide_samples-1))))/sqrt(60);
    end
end

UnitsPos_idx = find(PertRespUnits_pos==1);
figure
for i = 1:length(UnitsPos_idx)
    subplot(13,5,i)
    errorbar(U_speed,Resp_mean(UnitsPos_idx(i),:),Resp_sem(UnitsPos_idx(i),:))
    set(gca,'box','off','TickDir','out','xTick',0:1:9,'xTickLabel',[])
    xlim([0 9])
    if i>60
        set(gca,'box','off','TickDir','out','xTick',0:1:9,'xTickLabel',0:1:9)
        xlabel('TF @ 0^o');
    end 
    if sum(i == 1:5:length(UnitsPos_idx))>0
        ylabel('spikes/s');
    end
end

%% plot heat map of all 65 speed tuning curves
% positively modulated
UnitsPos_idx = find(PertRespUnits_pos==1);
Resp2Plot = Resp_mean(UnitsPos_idx,3:end);
Resp2Plot = Resp2Plot./max(Resp2Plot,[],2);
% fill in missing values
for i = 1:size(Resp2Plot,1)
    i_nan = find(isnan(Resp2Plot(i,:)));
    i_nnan = find(~isnan(Resp2Plot(i,:)));
    fillv = interp1([-flintmax i_nnan],[0 Resp2Plot(i,i_nnan)],i_nan,'next');
    Resp2Plot(i,i_nan) = fillv;
end     
[B I_inc] = sortrows(Resp2Plot,-([1:1:size(Resp2Plot,2)]));
figure
subplot(1,2,1)
imagesc(U_speed(3:end-3),1:size(Resp2Plot,1), Resp2Plot(:,1:end-3))
ylabel('Units'); xlabel('TF @ 0^o');
set(gca,'box','off','TickDir','out','xTick',0:1:max(U_speed(3:end-3)),'xTickLabel',0:1:max(U_speed(3:end-3)))
% negatively modulated
UnitsPos_idx = find(PertRespUnits_neg==1);
Resp2Plot = Resp_mean(UnitsPos_idx,3:end);
Resp2Plot = Resp2Plot./max(Resp2Plot,[],2);
% fill in missing values
for i = 1:size(Resp2Plot,1)
    i_nan = find(isnan(Resp2Plot(i,:)));
    i_nnan = find(~isnan(Resp2Plot(i,:)));
    fillv = interp1([-flintmax i_nnan],[0 Resp2Plot(i,i_nnan)],i_nan,'next');
    Resp2Plot(i,i_nan) = fillv;
end     
[B I_inc] = sortrows(Resp2Plot,-([1:1:size(Resp2Plot,2)]));
subplot(1,2,2)
imagesc(U_speed(3:end-3),1:size(Resp2Plot,1), Resp2Plot(:,1:end-3))
ylabel('Units'); xlabel('TF @ 0^o');
set(gca,'box','off','TickDir','out','xTick',0:1:max(U_speed(3:end-3)),'xTickLabel',0:1:max(U_speed(3:end-3)))
colormap parula
colorbar


%%
% evaluate mean firing rate at stimulus onset for the 3 diff starting
% speeds
StartSpeeds = unique(AM_Param_acc_s(:,1,2));
clear Resp_mean_1 Reqsp_sem_1
for i = 1:length(StartSpeeds)
    for j = 1:size(AM_UnitResponses_acc_s,2)   
        Indexes = AM_Param_acc_s(:,j,2) == StartSpeeds(i);
        Resp_mean_1(j,i) = mean(nanmean(squeeze(AM_UnitResponses_acc_s(Indexes,j,trialSide_samples+15:nanmedian(AM_Param_acc_s(Indexes,j,5)))))); 
        Resp_sem_1(j,i) = std(nanmean(squeeze(AM_UnitResponses_acc_s(Indexes,j,trialSide_samples+15:nanmedian(AM_Param_acc_s(Indexes,j,5))))))/sqrt(30);
    end
end
    

%% compute acceleration tuning curves
clear RespAcc_mean RespAcc_sem BL_acc1 BL_acc
for i = 1:length(dTF_values)
    for j = 1:size(AM_UnitResponses_acc_s,2)
        Indexes = AM_Param_acc_s(:,j,4) == dTF_values(i);
        RespAcc_mean(j,i) = mean(nanmean(squeeze(AM_UnitResponses_acc_s(Indexes,j, nanmedian(AM_Param_acc_s(Indexes,j,5)):nanmedian(AM_Param_acc_s(Indexes,j,6)) ))));
        RespAcc_sem(j,i) = std(nanmean(squeeze(AM_UnitResponses_acc_s(Indexes,j, nanmedian(AM_Param_acc_s(Indexes,j,5)):nanmedian(AM_Param_acc_s(Indexes,j,6)) ))))/sqrt(60);
        BL_acc1(j,1) = mean(squeeze(nanmean(AM_UnitResponses_acc_s(:,j, trialSide_samples+15:nanmedian(AM_Param_acc_s(Indexes,j,5)) ))));
        BL_acc1(j,2) = mean(squeeze(nanstd(AM_UnitResponses_acc_s(:,j, trialSide_samples+15:nanmedian(AM_Param_acc_s(Indexes,j,5)) )))/sqrt(size(AM_UnitResponses_acc_s,1)));
        BL_acc(j,1) = mean(squeeze(nanmean(AM_UnitResponses_acc_s(:,j,[1:trialSide_samples size(AM_UnitResponses_acc_s,3)-trialSide_samples+15:size(AM_UnitResponses_acc_s,3)]))));
        BL_acc(j,2) = mean(squeeze(nanstd(AM_UnitResponses_acc_s(:,j,[1:trialSide_samples size(AM_UnitResponses_acc_s,3)-trialSide_samples+15:size(AM_UnitResponses_acc_s,3)])))/sqrt(size(AM_UnitResponses_acc_s,1)));
    end
end

figure
for i = 1:length(UnitsPos_idx)
    subplot(13,5,i)
    errorbar(dTF_values,RespAcc_mean(UnitsPos_idx(i),:),RespAcc_sem(UnitsPos_idx(i),:))
end

save('AccelerationInfoCTRL.mat','U_speed','Resp_mean','Resp_sem','dTF_values','RespAcc_mean','RespAcc_sem','BL_acc','BL_acc1')

%% produce some plots
figure
title('example trial with starting TF=6Hz & dTF=-3:1:3; dTF duration 1s');
for i = 1:size(I_05,1)
    plot(ProjectData_acc.Trial_data{AM_Param(I_05(i,1),I_05(i,2),1),1}.time_in_trial{I_05(i,1),1}, ...
        ProjectData_acc.Trial_data{AM_Param(I_05(i,1),I_05(i,2),1),1}.TF_in_trial{I_05(i,1),1},'Color',[0.5 0.5 0.5])
    hold on
end
plot([0.5 0.5],[0 9.5],'r--','LineWidth',2); plot([1 1],[0 9.5],'r--','LineWidth',2) % dTF period
plot([0 0],[0 9.5],'b:','LineWidth',4); plot([3 3],[0 9.5],'b:','LineWidth',4); % trial period
set(gca,'TickDir','out','XTickLabel',[0 0.5 1 3])
box off
ylabel('TF'); xlabel('Seconds')
ylim([0 9.5]); xlim([-0.5 3.5]);
%text(0.25,1,'example trial with starting TF=6Hz & dTF=-3:1:3')
h = text(-0.1,2, 'trial start'); set(h,'Rotation',90); 
h = text(3.1,2, 'trial end'); set(h,'Rotation',90); 
h = text(0.45,2, 'dTF start'); set(h,'Rotation',90); 
h = text(1.55,2, 'dTF stop'); set(h,'Rotation',90); 
title('example trial with starting TF=6Hz & dTF=-3:1:3; dTF duration 0.5 s');
set(gca,'FontSize',13)


%% plot the 18 units that are perturbation responsive in the CTRL animals

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AM_Param :
% conditions = 1 --> nr of recording
% conditions = 2 --> grating direction [0:45:315]
% conditions = 3 --> 1 pert ON, 0 pert OFF
% conditions = 4 --> pert onset
% conditions = 5 --> pert offset

% find angle direction of each trial
Param = AM_Param(:,AM_UOI,2);
% smooth responses to the main stimulus
Gauss_width = 0.3;
AM_UnitResponses_smooth__ = Smooth_AM_FR(ProjectData,AM_UnitResponses,Gauss_width);
UnitResponses = AM_UnitResponses_smooth__(:,AM_UOI,:)*60;
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
%UnitResponses = AM_UnitResponses_smooth__(:,AM_UOI,:)*60;
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% plot them all
UnitsPos_idx = find(PertRespUnits_pos==1);

for  i = 1:length(UnitsPos_idx)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    UI_idx = UnitsPos_idx(i);
%%%%%%%%%%%%%%%%%%%
if DM(UI_idx)
    param_sel = AM_Param(:,SelectedCells,:);
    Pert_lims =[param_sel(min(find(param_sel(:,UI_idx,3)==1)),UI_idx,4) param_sel(min(find(param_sel(:,UI_idx,3)==1)),UI_idx,5)];
    TimeLine = linspace(-1,8.33,size(AM_UnitResponses_smooth__,3));
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
    
    % % % % % % %% % % % % % %% % % % % % %% % % % % % %
    % acceleration tuning
    subplot(14, 14, [1 2 3 4 15 16 17 18 29 30 31 32])
    errorbar(dTF_values,RespAcc_mean(UI_idx,:)-BL_acc1(UI_idx,1),RespAcc_sem(UI_idx,:),'Color','k')
    set(gca,'box','off','TickDir','out')
    xlabel('Temporal dTF (cycles/s/s)'); ylabel('Firing rate (spikes/s)')
    % % % % % % %% % % % % % %% % % % % % %% % % % % % %
    % speed tuning
    subplot(14, 14, [1 2 3 4 15 16 17 18 29 30 31 32]+10)
    %errorbar(U_speed, Resp_mean(UI_idx,:)-BL_acc(UI_idx,1), Resp_sem(UI_idx,:),'Color','k')
    errorbar(U_speed, Resp_mean(UI_idx,:), Resp_sem(UI_idx,:),'Color','k')
    set(gca,'box','off','TickDir','out')
    xlabel('Temporal TF (cycles/s)'); ylabel('Firing rate (spikes/s)')
    xlim([0 max(U_speed)])
    
    saveas(gcf,['E:\Dropbox\UCL\Tomaso-Aman\VisPerturbation\Fig_3\Acceleration\PertUnits_PertResponses_new' filesep 'Unit_' num2str(UI_idx) '_noBL.pdf']);
    %   saveas(gcf,['E:\Dropbox\UCL\Tomaso-Aman\VisPerturbation\Fig_3\panels\Naive\gratingDirExamples' filesep 'Unit_' num2str(UI_idx) '.pdf']);
    %
    %         close all
    %     saveas(gcf,['X:\DATA\PROJECTS\VisPerturbation\Figures\Fig_2\panels\ExampleUnits' filesep FileSuffix '_' num2str(UI) '.fig']);
    %
    %input('');
    close all
end


end



%% compute where the max response to acceleration is recorded, all units and perturbation units
% count nr of trials of each acceleration presentation
for i = 1:length(dTF_values)
    for j = 1:size(AM_UnitResponses_smooth_acc,2)
        IndexesCount_accel(i,j) = sum(AM_Param_acc(:,j,4) == dTF_values(i));
    end
end
%%%%
for i = 1:length(RespAcc_mean)
    [i_accmax(i,1) i_accmax(i,2)] = max(RespAcc_mean(i,:));
end

%AccResp2plot = i_accmax(AM_UOI,:);
AccResp2plot = i_accmax;
figure
histogram(AccResp2plot(~PertRespUnits_pos,2),'FaceColor',[0.5 0.5 0.5],'Normalization','probability','EdgeColor','none')
hold on
histogram(AccResp2plot(PertRespUnits_pos,2),'FaceColor',[1 0 0],'EdgeColor','none','Normalization','probability')
hold on
set(gca,'XTickLabel',dTF_values)
xlabel('Preferred dTF'); ylabel('Probability')
box off
set(gca,'FontSize',13,'TickDir','out');

dist=KLDiv(AccResp2plot(~PertRespUnits_pos,2),AccResp2plot(PertRespUnits_pos,2))
%  dist = KLDiv(P,Q) Kullback-Leibler divergence of two discrete probability
%  distributions
%  P and Q  are automatically normalised to have the sum of one on rows
% have the length of one at each 
% P =  n x nbins
% Q =  1 x nbins or n x nbins(one to one)
% dist = n x 1

%% compute where the max response to speed is recorded, all units and perturbation units.
% count nr of trials of each speed presentation
for i = 1:length(U_speed)
    for j = 1:size(AM_UnitResponses_smooth_acc,2)
        IndexesCount_speed(i,j) = sum(AM_Param_acc(:,j,3) == U_speed(i));
    end
end
IndexesCount_speed = IndexesCount_speed(3:end,:);
%%%
clear i_speedmax
for i = 1:length(Resp_mean)
    [i_speedmax(i,1) i_speedmax(i,2)] = nanmax(Resp_mean(i,3:end)-BL_acc(i,1));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SpeedResp2plot = i_speedmax(AM_UOI,:);
SpeedResp2plot = i_speedmax;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
histogram(SpeedResp2plot(~PertRespUnits_pos,2),'FaceColor',[0.5 0.5 0.5],'EdgeColor','none')
hold on
histogram(SpeedResp2plot(PertRespUnits_pos,2),'FaceColor',[1 0 0],'EdgeColor','none')
hold on
set(gca,'XTick',[1:length(U_speed(3:end))],'XTickLabel',U_speed(3:end))
xlabel('Preferred TF'); ylabel('Probability')
box off
set(gca,'FontSize',13,'TickDir','out')%,'YTick',[1 2 3 4]);

%% compare max response to perturbation VS max response to drifting grating
SpeedResp2plot = i_speedmax;%(AM_UOI,:);
figure
histogram((GratingDirResp(~PertRespUnits_pos,1)-SpeedResp2plot(~PertRespUnits_pos,1))./ ...
          (GratingDirResp(~PertRespUnits_pos,1)+SpeedResp2plot(~PertRespUnits_pos,1)) ...
    ,'FaceColor',[0.5 0.5 0.5],'EdgeColor','none','Normalization','probability')
hold on
histogram((GratingDirResp(PertRespUnits_pos,1)-SpeedResp2plot(PertRespUnits_pos,1))./ ...
        (GratingDirResp(PertRespUnits_pos,1)+SpeedResp2plot(PertRespUnits_pos,1)) ...
          ,'FaceColor',[1 0 0],'EdgeColor','none','Normalization','probability')
xlabel('Perturbation response @ 0^o - Preferred TF (spikes/s)');
ylabel('Probability')
box off
set(gca,'FontSize',13,'TickDir','out')


figure
plot(GratingDirResp(~PertRespUnits_pos,1),SpeedResp2plot(~PertRespUnits_pos,1),'.k')
hold on
plot(GratingDirResp(PertRespUnits_pos,1),SpeedResp2plot(PertRespUnits_pos,1),'or')
plot(GratingDirResp(PertRespUnits_pos,1),SpeedResp2plot(PertRespUnits_pos,1),'.r')
xlabel('Perturbation response @ 0^o (spikes/s)');
ylabel('Preferred TF (spikes/s)')
box off
set(gca,'FontSize',13,'TickDir','out')
hold on
plot([0 16],[0 16],'k')
xlim([-10 40]);ylim([-10 40])

%% verify where the max grating response is to check the negative acceleration...

TuningPropsV.P_dir(~PertRespUnits_pos)
% select units with singificant direction tuning
Dir__idx = TuningPropsV.p_tt<=0.01;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% perturbation responsive units
DirResp_idx = find(PertRespUnits_pos & Dir__idx);
Angle_v = TuningPropsV.P_dir(PertRespUnits_pos & Dir__idx); % angles at which pert. resp. units respond significantly to grating
Accel_v = AccResp2plot(PertRespUnits_pos & Dir__idx,2); % index of dTF values where max response was recorded
clear Acc_Value_idx
Acc_Value_idx = false(length(Angle_v),3); % three as grouping of the directions
Acc_Value_idx(:,1) = Angle_v>=315 | Angle_v<45;
Acc_Value_idx(:,2) = Angle_v>=135 & Angle_v<225;
Acc_Value_idx(:,3) = (Angle_v>=45 & Angle_v<135) | (Angle_v>=225 & Angle_v<315);

% now I have the nr of units with their max responses. Plot as stacked bars
% the # of units that make each accel value!!!
clear SumAll
for i = 1:length(dTF_values)    
        SumAll(i,:) = sum(Acc_Value_idx(dTF_values(Accel_v) == dTF_values(i),:),1);
end
figure
bar(SumAll,'grouped')
set(gca,'XTickLabel',dTF_values)
xlabel('Preferred dTF'); ylabel('Units')
box off
set(gca,'FontSize',13,'TickDir','out');
legend({'Temporal','Nasal','Vertical'})
title(['Pert.Resp. & Dir. tuned Units, n=' num2str(length(Acc_Value_idx))])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rest of the units
DirResp_idx = find(~PertRespUnits_pos & Dir__idx);
Angle_v = TuningPropsV.P_dir(~PertRespUnits_pos & Dir__idx); % angles at which pert. resp. units respond significantly to grating
Accel_v = AccResp2plot(~PertRespUnits_pos & Dir__idx,2); % index of dTF values where max response was recorded
clear Acc_Value_idx
Acc_Value_idx = false(length(Angle_v),3); % three as grouping of the directions
Acc_Value_idx(:,1) = Angle_v>=315 | Angle_v<45;
Acc_Value_idx(:,2) = Angle_v>=135 & Angle_v<225;
Acc_Value_idx(:,3) = (Angle_v>=45 & Angle_v<135) | (Angle_v>=225 & Angle_v<315);
% now I have the nr of units with their max responses. Plot as stacked bars
% the # of units that make each accel value!!!
clear SumAll
for i = 1:length(dTF_values)    
        SumAll(i,:) = sum(Acc_Value_idx(dTF_values(Accel_v) == dTF_values(i),:),1);
end
figure
bar(SumAll,'grouped')
set(gca,'XTickLabel',dTF_values)
xlabel('Preferred dTF'); ylabel('Units')
box off
set(gca,'FontSize',13,'TickDir','out');
legend({'Temporal','Nasal','Vertical'})
title(['Not Pert.Resp. & Dir. tuned Units, n=' num2str(length(Acc_Value_idx))])

%% compare max response to perturbation in pref direction VS max response to drifting grating

TuningPropsV.P_dir(PertRespUnits_neg)
% select units with singificant direction tuning
Dir__idx = TuningPropsV.p_tt<0.01;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Negatively perturbation responsive units
DirResp_idx = find(PertRespUnits_pos & Dir__idx);
Angle_v = TuningPropsV.P_dir(PertRespUnits_neg & Dir__idx); % angles at which pert. resp. units respond significantly to grating
Speed_v = SpeedResp2plot(PertRespUnits_neg & Dir__idx,2); % index of TF values where max response was recorded
clear Speed_Value_idx
Speed_Value_idx = false(length(Angle_v),3); % three as grouping of the directions
Speed_Value_idx(:,1) = Angle_v>=315 | Angle_v<45;
Speed_Value_idx(:,2) = Angle_v>=135 & Angle_v<225;
Speed_Value_idx(:,3) = (Angle_v>=45 & Angle_v<135) | (Angle_v>=225 & Angle_v<315);
% now I have the nr of units with their max responses. Plot as stacked bars
% the # of units that make each speed value!!!
clear SumAll
U_speed_v = U_speed(3:end);
for i = 1:length(U_speed_v)    
        SumAll(i,:) = sum(Speed_Value_idx(U_speed_v(Speed_v) == U_speed_v(i),:),1);
end

% group together the speed values 0&0.5, 1&1.5, etc
A = SumAll'
B = reshape(A,3,2,9)
C = squeeze(sum(B,2))'

figure
bar(C,'grouped')
set(gca,'XTick',1:length(U_speed_v)/2)
xlabel('Preferred TF'); ylabel('Units')
box off
set(gca,'FontSize',13,'TickDir','out');
legend({'Temporal','Nasal','Vertical'})
title(['Pert.Resp. & Dir. tuned Units, n=' num2str(C(:))])
title(['Only Dir. tuned Units, n=' num2str(sum(C(:)))])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% perturbation responsive units
Angle_v = TuningPropsV.P_dir(PertRespUnits_pos & Dir__idx); % angles at which pert. resp. units respond significantly to grating
Speed_v = SpeedResp2plot(PertRespUnits_pos & Dir__idx,2); % index of TF values where max response was recorded
clear Speed_Value_idx
Speed_Value_idx = false(length(Angle_v),3); % three as grouping of the directions
Speed_Value_idx(:,1) = Angle_v>=315 | Angle_v<45;
Speed_Value_idx(:,2) = Angle_v>=135 & Angle_v<225;
Speed_Value_idx(:,3) = (Angle_v>=45 & Angle_v<135) | (Angle_v>=225 & Angle_v<315);
% now I have the nr of units with their max responses. Plot as stacked bars
% the # of units that make each speed value!!!
clear SumAll
U_speed_v = U_speed(3:end);
for i = 1:length(U_speed_v)    
        SumAll(i,:) = sum(Speed_Value_idx(U_speed_v(Speed_v) == U_speed_v(i),:),1);
end

% group together the speed values 0&0.5, 1&1.5, etc
A = SumAll';
B = reshape(A,3,2,9)
C1 = squeeze(sum(B,2))'

figure
bar(C1,'grouped')
set(gca,'XTick',1:length(U_speed_v)/2)
xlabel('Preferred TF'); ylabel('Units')
box off
set(gca,'FontSize',13,'TickDir','out');
legend({'Temporal','Nasal','Vertical'})
title(['Pert.Resp. & Dir. tuned Units, n=' num2str(sum(C(:)))])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SumAll_PN = C;%/sum(C(:));
SumAll_PP = C1;%/sum(C1(:));
clear stackData
stackData(:,:,1) = [SumAll_PP(:,1) SumAll_PN(:,1)]; % (Group, Stack, StackElement)
stackData(:,:,2) = [SumAll_PP(:,2) SumAll_PN(:,3)];
stackData(:,:,3) = [SumAll_PP(:,2) SumAll_PN(:,3)];
groupLabels = {'PP', 'PN'};
plotBarStackGroups(stackData, groupLabels)
set(gca,'XTick',1:length(U_speed_v)/2)
xlabel('Preferred TF'); ylabel('Probability')
box off
set(gca,'FontSize',13,'TickDir','out');
legend({'Temporal','Nasal','Vertical'})
title(['Pert.Resp. & Dir. tuned Units, n=' num2str(sum(SumAll(:)))])
title(['Only Dir. tuned Units, n=' num2str(sum(SumAll(:)))])


%% plot MI vs preferred TF at 0 deg for only dir. tuned units

% additional selection criteria for reliable units
% apply 0.1 Hz threshold also for FR during acceleration stimulus
I_main = find(AM_UOI);
I_accel = find(AM_UOI_acc);
for i = 1:length(I_main)
    for j = 1:length(I_accel)
        if I_main(i) == I_accel(j)
            I_bi(i) = true;
        end
    end
end
% try to increase it to see what you get!
% manually discard those units?

% select all units responding to perturbation 
p_pert_th_01 = prctile(Sh_responses(:),99);
PertResp_units_01 = (AUC_shuffled(:,1)>p_pert_th_01);
p_pert_th_05 = prctile(Sh_responses(:),95);
PertResp_units_05 = (AUC_shuffled(:,1)>p_pert_th_05);

I_bi = true(1,length(Resp_mean));
% select preferred TF frequency
Resp_mean_bi = Resp_mean(I_bi,:);
[Pref_TF TF_i] = max(Resp_mean_bi,[],2);

DM_bi = DM(I_bi);
PertResp_units_01 = PertResp_units_01(I_bi);
PertResp_units_05 = PertResp_units_05(I_bi);

Dir__idx = TuningPropsV.p_tt<0.01;
Dir__idx__ = Dir__idx(I_bi);

figure
% p<0.01
subplot(2,2,1)
%DM_2plot = DM_bi(PertResp_units_01 & Dir__idx__); U_speed2plot = max(U_speed(TF_i(PertResp_units_01 & Dir__idx__)),0);
DM_2plot = DM_bi(PertResp_units_01); U_speed2plot = max(U_speed(TF_i(PertResp_units_01)),0);
plot(U_speed2plot(DM_2plot>0),DM_2plot(DM_2plot>0), '.r','MarkerSize',15)
hold on
plot(U_speed2plot(DM_2plot<0),DM_2plot(DM_2plot<0), '.k','MarkerSize',15)
ylabel('Perturbation MI_F_R')
ylim([-1 1])
set(gca,'box','off','TickDir','out','XTickLabel',[],'YTick',[-1 0 1])
grid on;
legend({['MI>0 n=' num2str(length(DM_2plot(DM_2plot>0)))],['MI<0 n=' num2str(length(DM_2plot(DM_2plot<0)))]})
subplot(2,2,1)
title('All pert. units (AUC - p<0.01)')
subplot(2,2,3)
histogram(U_speed2plot(DM_2plot>0),[0:0.5:9],'FaceColor',[1 0 0],'EdgeColor','none')
hold on
histogram(U_speed2plot(DM_2plot<0),[0:0.5:9],'FaceColor',[0 0 0],'EdgeColor','none')
xlabel('Preferred TF @ 0^o');
ylabel('Units');
set(gca,'box','off','TickDir','out')
subplot(2,2,3)
p = ranksum(U_speed2plot(DM_2plot>0),U_speed2plot(DM_2plot<0))
title(['ranksum test p=' num2str(p,2)]);
% p<0.05
subplot(2,2,2)
%DM_2plot = DM_bi(PertResp_units_05 & Dir__idx__); U_speed2plot = max(U_speed(TF_i(PertResp_units_05 & Dir__idx__)),0);
DM_2plot = DM_bi(PertResp_units_05); U_speed2plot = max(U_speed(TF_i(PertResp_units_05)),0);
plot(U_speed2plot(DM_2plot>0),DM_2plot(DM_2plot>0), '.r','MarkerSize',15)
hold on
plot(U_speed2plot(DM_2plot<0),DM_2plot(DM_2plot<0), '.k','MarkerSize',15)
ylabel('Perturbation MI_F_R')
ylim([-1 1])
set(gca,'box','off','TickDir','out','XTickLabel',[],'YTick',[-1 0 1])
grid on;
legend({['MI>0 n=' num2str(length(DM_2plot(DM_2plot>0)))],['MI<0 n=' num2str(length(DM_2plot(DM_2plot<0)))]})
subplot(2,2,2)
title('All pert. units (AUC - p<0.05)')
subplot(2,2,4)
histogram(U_speed2plot(DM_2plot>0),[0:0.5:9],'FaceColor',[1 0 0],'EdgeColor','none')
hold on
histogram(U_speed2plot(DM_2plot<0),[0:0.5:9],'FaceColor',[0 0 0],'EdgeColor','none')
xlabel('Preferred TF @ 0^o');
ylabel('Units');
set(gca,'box','off','TickDir','out')
subplot(2,2,4)
p = ranksum(U_speed2plot(DM_2plot>0),U_speed2plot(DM_2plot<0))
title(['ranksum test p=' num2str(p,2)]);





