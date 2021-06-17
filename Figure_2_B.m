%% Tomaso Muzzu - UCL - 06/03/2020 - DIgure 2B

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

%% FIGURE 2B
%%%% with components congruent with running direction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AngleCombo = [0 ; 180]; 
AngleCombo = [45 ; 225];
AngleCombo = [90 ; 270];
AngleCombo = [135 ; 315];

AngleCombo = [0 180; 90 270];
AngleCombo = [0 45 315; 135 180 225];

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



%PlotResponses(ProjectData,{SelectedResponses{1,1}(UnitsPertPos,:) SelectedResponses{2,1}(UnitsPertPos,:)},AM_Param,AM_Speed,UnitsPertPos,0); 
%PlotResponses(ProjectData,{SelectedResponses{1,1}(UnitsPertNeg,:) SelectedResponses{2,1}(UnitsPertNeg,:)},AM_Param,AM_Speed,UnitsPertNeg,0); 

%PlotMeanCurveAngle(ProjectData,{SelectedResponses{1,1}(UnitsPertPos,:) SelectedResponses{2,1}(UnitsPertPos,:)},AM_Param,AM_Speed,UnitsPertPos,AngleCombo)
%PlotMeanCurveAngle(ProjectData,{SelectedResponses{1,1}(UnitsPertNeg,:) SelectedResponses{2,1}(UnitsPertNeg,:)},AM_Param,AM_Speed,UnitsPertNeg,AngleCombo)

%% PLOT CURVES ON SESSION BY SESSION BASIS. ONLY IF NECESSARY
% show different angle pop. responses for stimulus onset
% show different angle pop. responses for perturbation onset
% all units and only pert. responsive units
param_sel = AM_Param(:,AM_UOI,:);
RecordingsNR = unique(param_sel(1,:,1));
for i = 1:length(RecordingsNR)
    RecordingsNR_i(i) = min(find(param_sel(1,:,1)==RecordingsNR(i))); % find first index of new recording in 3D Matrix
end
GratingDir = 0:45:315;
AM_UnitResponses_sm = AM_UnitResponses_smooth(:,SelectedCells,:)*BonvisionFR; % responses of all units in spikes/s
clear AngResponse AngResponsePert AngResponsePertOFF
for Rec_OI = 1:length(RecordingsNR_i)
    % responses of all units session by session
    clear SingleSession
    if Rec_OI == length(RecordingsNR_i)
        SingleSession = AM_UnitResponses_sm(:,RecordingsNR_i(Rec_OI):end,:);
    else
        SingleSession = AM_UnitResponses_sm(:,RecordingsNR_i(Rec_OI):RecordingsNR_i(Rec_OI+1)-1,:); 
    end
    for Angle_OI = 1:length(GratingDir) % scan through the directions
        % save the mean responses for each angle of all units in a struct
        AngResponse{Rec_OI,Angle_OI} = SingleSession(param_sel(:,RecordingsNR_i(Rec_OI),2)==GratingDir(Angle_OI),:,:); % all trials
        AngResponsePert{Rec_OI,Angle_OI} = SingleSession(param_sel(:,RecordingsNR_i(Rec_OI),2)==GratingDir(Angle_OI) & ...
                                                     param_sel(:,RecordingsNR_i(Rec_OI),3)==1,:,:); % only perturbation trials
        AngResponsePertOFF{Rec_OI,Angle_OI} = SingleSession(param_sel(:,RecordingsNR_i(Rec_OI),2)==GratingDir(Angle_OI) & ...
                                                     param_sel(:,RecordingsNR_i(Rec_OI),3)==0,:,:); % only non perturbation trials
    end
end


%% Directional responses of ensembles of pert. responsive units
% define variables of the stimulus
PertTrial_Ex(1) = 1;
PertTrial_Ex(2) = min(find(ProjectData.Trial_data{1,1}.PerturbationON==1))
TimeLine = linspace(min(ProjectData.Trial_data{PertTrial_Ex(1),1}.time_in_trial{PertTrial_Ex(2)}), ...
                    max(ProjectData.Trial_data{PertTrial_Ex(1),1}.time_in_trial{PertTrial_Ex(2)}),...
                    size(AM_UnitResponses_sm,3)) -trialSide_seconds; % -1 seconds
Contrast = ProjectData.Trial_data{PertTrial_Ex(1),1}.contrast_in_trial{PertTrial_Ex(2),1}(1:end-1);
TF = ProjectData.Trial_data{PertTrial_Ex(1),1}.TF_in_trial{PertTrial_Ex(2),1}(1:end-1);
TrialEnd = TimeLine(min(find(Contrast(length(Contrast)/2:end)==0))+round(length(Contrast)/2));

PertRespUnits_OI = PertRespUnits_pos;
Sessions_OI = unique(param_sel(1,PertRespUnits_OI,1)); % only sessions with positively modulated units
% Sessions_OI = unique(param_sel(1,PertRespUnits_neg,1)); % only sessions with positively modulated units

% array with 2 columns for ID-ing recording of all perturbation resp. unis
clear Sessions_Units
Sessions_Units(:,2) = param_sel(1,:,1);
Sessions_Units(:,1) = PertRespUnits_OI;
clear YValues1 YValues2 EnsembleSize YValues2_sem AllUnitsAnglePertResp
for i = 1:length(Sessions_OI) % scans through all sessions of interest
    clear PopResponses
    Units2plot = find(Sessions_Units(find(Sessions_Units(:,2) == Sessions_OI(i)),1)==1); % find all pert. resp. units of the current recording
    for j = 1:size(AngResponsePert,2) % scans through all directions
        Response = AngResponsePert{Sessions_OI(i),j}; % population response of current direction for current recording
        Response = squeeze(nanmean(Response(:,Units2plot,:),1)); % compute the mean of the responses for each pert. resp. unit
        Response = Response./max(Response,[],2); 
        if size(Response,2)>1
            PopResponses(j,:) = mean(Response(~isnan(sum(Response,2)),:),1);
            AllUnitsAnglePertResp{i,:}(j,:) = mean(Response(~isnan(sum(Response,2)),:),2);
        else
            PopResponses(j,:) = Response;
        end
    end
        % stimulus onset
        YValues1(i,:) = [mean(PopResponses(:,trialSide_samples:trialSide_samples+4*BonvisionFR),2)- ...
                    mean(PopResponses(:,1:trialSide_samples),2)]';
        % perturbation onset
        %YValues2(i,:) = [mean(PopResponses(:,min(find(TF==0)):max(find(TF==0))),2)- ...
        %            mean(PopResponses(:,min(find(TF==0))-60:min(find(TF==0))),2)]';
        % use inter trial interval as baseline for perturbation responses
        YValues2(i,:) = [mean(PopResponses(:,min(find(TF==0))+1:min(find(TF==0))+trialSide_samples)- ...
                    PopResponses(:,1:trialSide_samples),2)]';
        YValues2_sem(i,:) = [std((PopResponses(:,min(find(TF==0))+1:min(find(TF==0))+trialSide_samples) - ...
                    PopResponses(:,1:trialSide_samples))')/sqrt(length(size(PopResponses,1)))];                
        EnsembleSize(i) = min(size(Response));
        
end

% figure
% boxplot(YValues2(find(EnsembleSize~=1),1)-YValues2(find(EnsembleSize~=1),:),'Notch','on','Labels',{'0','45','90','135','180','225','270','315'},'symbol','')
figure
Data2Scatter = YValues2(find(EnsembleSize~=1),1)-YValues2(find(EnsembleSize~=1),:);
UnivarScatter(Data2Scatter(:,2:end),'Label',{'45','90','135','180','225','270','315'},'MarkerFaceColor',[ 0.5 0.5 0.5])
xlabel('Perturbation direction');
ylabel('Differential response');
set(gca,'TickDir','out','box','off')
title('Direction')
ylim([-max(max(abs(YValues2(find(EnsembleSize~=1),1)-YValues2(find(EnsembleSize~=1),:)))) ...
       max(max(abs(YValues2(find(EnsembleSize~=1),1)-YValues2(find(EnsembleSize~=1),:))))]);

%% plot heat map of population tuning curves for all recordings   
% sort them by preference for 0 degree
PopResponses2plot = YValues2(find(EnsembleSize~=1),:) - YValues2(find(EnsembleSize~=1),1);
[B I_inc] = sortrows(PopResponses2plot,-([2 8 3 7 4 6 5]));
figure
Range = [4:8; 1:5];
imagesc([ PopResponses2plot(I_inc,Range(1,:)) PopResponses2plot(I_inc,Range(2,:))])
ylabel('Recordings'); xlabel('Direction');
set(gca,'box','off','TickDir','out','xTick',1:length(Range(:)), 'xTickLabel',[GratingDir(Range(1,:)) GratingDir(Range(2,:))])
% extend x axis and place the 0 degree in the middle
% indicate nr of units and mouse id
EnsembleSize2Plot = EnsembleSize(find(EnsembleSize~=1));
EnsembleSize2Plot = EnsembleSize2Plot(I_inc);
for i = 1:size(PopResponses2plot,1)
    text(1,i,num2str(EnsembleSize2Plot(i)));
end
title('Population response minus response at 0^o')
hold on
plot([length(Range(:))/2 length(Range(:))/2] - 0.5, [1 size(PopResponses2plot,1)],'Color',[0.5 0.5 0.5])
plot([length(Range(:))/2 length(Range(:))/2] + 2.5, [1 size(PopResponses2plot,1)],'Color',[0.5 0.5 0.5])
redblue_WhiteIsZero(PopResponses2plot)

figure
PopResponses2plot_m = YValues2(find(EnsembleSize~=1),:);
PopResponses2plot_s = YValues2_sem(find(EnsembleSize~=1),:);
AllUnitsAnglePertResp2plot = AllUnitsAnglePertResp(find(EnsembleSize~=1));
for i = 1:length(PopResponses2plot_m)
    subplot(6,7,i)    
    hold on
%     for j = 1:size(AllUnitsAnglePertResp2plot{i,:},2)
%         plot(0:45:315,AllUnitsAnglePertResp2plot{i,:}(:,j),'Color',[0.5 0.5 0.5])
%     end
    errorbar(0:45:315, PopResponses2plot_m(i,:),PopResponses2plot_s(i,:),'k')
    set(gca,'box','off','TickDir','out','xTick',0:90:315)
    if i>34
        xlabel('Direction');
    end
    if mod(i,7)==1 
        ylabel('Norm. response')
    end
end


figure
YValuesOri =   [[YValues2(find(EnsembleSize~=1),1)+YValues2(find(EnsembleSize~=1),5)]/2, ...
                [YValues2(find(EnsembleSize~=1),2)+YValues2(find(EnsembleSize~=1),6)]/2, ...
                [YValues2(find(EnsembleSize~=1),3)+YValues2(find(EnsembleSize~=1),7)]/2, ...
                [YValues2(find(EnsembleSize~=1),4)+YValues2(find(EnsembleSize~=1),8)]/2];
% boxplot(YValuesOri(:,1) - YValuesOri,'Notch','on','Labels',{'0','45','90','135'},'symbol','')
Data2Scatter = YValuesOri(:,1) - YValuesOri;
UnivarScatter(Data2Scatter(:,2:end),'Label',{'45','90','135'},'MarkerFaceColor',[ 0.5 0.5 0.5])
% hold on
% for i = 1:size(YValuesOri,2)
%     plot(i-0.1+0.2*rand(size(YValuesOri,1),1), YValuesOri(:,1)-YValuesOri(:,i), '.', 'Color',[0.5 0.5 0.5]);
% end
title('Orientation')
xlabel('Perturbation orientation');
ylabel('Differential response');
set(gca,'TickDir','out','box','off')
hold on
ylim([-max(max(abs(YValuesOri(:,1) - YValuesOri))) max(max(abs(YValuesOri(:,1) - YValuesOri)))])

[B I_inc] = sortrows(YValuesOri,-([1 2 3 4]));
figure
Range = [2:4; 1:3];
imagesc([ YValuesOri(I_inc,Range(1,:)) YValuesOri(I_inc,Range(2,:))]-YValuesOri(I_inc,1))
ylabel('Recordings'); xlabel('Orientation');
set(gca,'box','off','TickDir','out','xTick',1:length(Range(:)), 'xTickLabel',[GratingDir(Range(1,:)) GratingDir(Range(2,:))])
% extend x axis and place the 0 degree in the middle
% indicate nr of units and mouse id
EnsembleSize2Plot = EnsembleSize(find(EnsembleSize~=1));
EnsembleSize2Plot = EnsembleSize2Plot(I_inc);
for i = 1:size(YValuesOri,1)
    text(1,i,num2str(EnsembleSize2Plot(i)));
end
title('Population response minus response at 0^o')
redblue_WhiteIsZero(PopResponses2plot)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% select only trials while the mouse runs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trialSide_samples = 60;
trialSide_seconds = 1;

PertOnsets = AM_Param(:,:,4); % 2D matrix
PertOffsets = AM_Param(:,:,5); % 2D matrix
[v min_i(1)] = max(PertOnsets(:)); % find the latest moment when pert onset happens
[v min_i(2)] = max(PertOffsets(:)); % find the latest moment when pert offset happens

ReferSize = [size(AM_Param,1),size(AM_Param,2)];
[trial_el(1), trial_el(2)] = ind2sub(ReferSize,min_i(1)); % 2D index of minimum pert onset time

Rec = AM_Param(trial_el(1), trial_el(2),1);
% define timeline of example recording
TrialStart = ProjectData.Session_data{Rec,1}.ACInfo{1,1}.trialStartsEnds(trial_el(1),1);
TrialEnd = ProjectData.Session_data{Rec,1}.ACInfo{1,1}.trialStartsEnds(trial_el(1),2)-TrialStart;
TrialStart = 0;
TimeLine = linspace(TrialStart-trialSide_seconds, ...
                    TrialEnd+trialSide_seconds,...
                    size(AM_UnitResponses_smooth,3)); % -1 seconds
[v po_i(1)] = min(abs(TimeLine-min(PertOnsets(:))));
[v po_i(2)] = min(abs(TimeLine-min(PertOffsets(:))));
SelectedCells = AM_UOI;
UnitResponses_smooth = AM_UnitResponses_smooth(:,SelectedCells,:) ;

%% select responses for running and still trials
% select trials of interest and control trials as well
Trials_PertON = AM_Param(:,:,3)==1; % find indexes where perturbation is on
Trials_PertOFF= AM_Param(:,:,3)==0; % find indexes where perturbation is off
SelectedCells = AM_UOI;
[AM_SpeedResponse AM_SpeedResponseControl] = ComputeSpeedResponse(AM_Param, AM_Speed,SelectedCells,Trials_PertON,Trials_PertOFF,ProjectData);
% AM_SpeedResponse(:,:,1) = mean speed -0.5s to +2s at visual stimulus onset
% AM_SpeedResponse(:,:,2) = mean speed -0.5s to end of perturbation period
% AM_SpeedResponse(:,:,3) = mean speed -0.5s to +1s at perturbation offset
% AM_SpeedResponse(:,:,4) = mean speed across entire trial
Run_TH = 2;
minTrials = 1; % min nr of trials for each condition
[SelectedTrials_OI SelectedTrials_OI_control TotTrials] = SpeedFiltering(AM_Speed, AM_SpeedResponse, AM_SpeedResponseControl, Trials_PertON,Trials_PertOFF,...
                                                                        Run_TH, minTrials);
                                                                    
Param_run = SelectedTrials_OI{2,1}(:,SelectedCells);
Param_still = SelectedTrials_OI{2,2}(:,SelectedCells);

%% Directional responses of ensembles of pert. responsive units
param_sel = AM_Param(:,AM_UOI,:);
RecordingsNR = unique(param_sel(1,:,1));
for i = 1:length(RecordingsNR)
    RecordingsNR_i(i) = min(find(param_sel(1,:,1)==RecordingsNR(i)));
end
GratingDir = 0:45:315;
clear AngResponse AngResponsePert AngResponsePertOFF
for Rec_OI = 1:length(RecordingsNR_i)
    % responses of all units 
    clear SingleSession
    if Rec_OI == length(RecordingsNR_i)
        SingleSession = AM_UnitResponses_sm(:,RecordingsNR_i(Rec_OI):end,:);
    else
        SingleSession = AM_UnitResponses_sm(:,RecordingsNR_i(Rec_OI):RecordingsNR_i(Rec_OI+1)-1,:);
    end
    for Angle_OI = 1:length(GratingDir) 
        % save the mean responses for each angle of all units in a struct
        AngResponse{Rec_OI,Angle_OI} = SingleSession(param_sel(:,RecordingsNR_i(Rec_OI),2)==GratingDir(Angle_OI) & ...
                                                      Param_run(:,RecordingsNR_i(Rec_OI)),:,:);
        AngResponsePert{Rec_OI,Angle_OI} = SingleSession(param_sel(:,RecordingsNR_i(Rec_OI),2)==GratingDir(Angle_OI) & ...
                                                     param_sel(:,RecordingsNR_i(Rec_OI),3)==1 & ...
                                                     Param_run(:,RecordingsNR_i(Rec_OI)),:,:);
        AngResponsePertOFF{Rec_OI,Angle_OI} = SingleSession(param_sel(:,RecordingsNR_i(Rec_OI),2)==GratingDir(Angle_OI) & ...
                                                     param_sel(:,RecordingsNR_i(Rec_OI),3)==0 & ...
                                                     Param_run(:,RecordingsNR_i(Rec_OI)),:,:);
    end
end

PertTrial_Ex(1) = 1;
PertTrial_Ex(2) = min(find(ProjectData.Trial_data{1,1}.PerturbationON==1))
TimeLine = linspace(min(ProjectData.Trial_data{PertTrial_Ex(1),1}.time_in_trial{PertTrial_Ex(2)}), ...
                    max(ProjectData.Trial_data{PertTrial_Ex(1),1}.time_in_trial{PertTrial_Ex(2)}),...
                    size(AM_UnitResponses_sm,3)) -trialSide_seconds; % -1 seconds
Contrast = ProjectData.Trial_data{PertTrial_Ex(1),1}.contrast_in_trial{PertTrial_Ex(2),1}(1:end-1);
TF = ProjectData.Trial_data{PertTrial_Ex(1),1}.TF_in_trial{PertTrial_Ex(2),1}(1:end-1);
TrialEnd = TimeLine(min(find(Contrast(length(Contrast)/2:end)==0))+round(length(Contrast)/2));

Sessions_OI = unique(param_sel(1,PertRespUnits_pos,1));

clear Sessions_Units
Sessions_Units(:,2) = param_sel(1,:,1);
Sessions_Units(:,1) = PertRespUnits_pos;
clear YValues1 YValues2 
for i = 1:length(Sessions_OI) 
    clear PopResponses
    Units2plot = find(Sessions_Units(find(Sessions_Units(:,2) == Sessions_OI(i)),1)==1);
    for j = 1:size(AngResponse,2) 
        Response = AngResponse{Sessions_OI(i),j};
        Response = squeeze(nanmean(Response(:,Units2plot,:),1));
        if size(Response,2)>1
            PopResponses(j,:) = mean(Response(~isnan(sum(Response,2)),:),1);
        else
            PopResponses(j,:) = Response;
        end
    end
        % stimulus onset
        YValues1(i,:) = [mean(PopResponses(:,trialSide_samples:trialSide_samples+4*BonvisionFR),2)- ...
                    mean(PopResponses(:,1:trialSide_samples),2)]';
        % perturbation onset
        %YValues2(i,:) = [mean(PopResponses(:,min(find(TF==0)):max(find(TF==0))),2)- ...
        %            mean(PopResponses(:,min(find(TF==0))-60:min(find(TF==0))),2)]';
        % use inter trial interval as baseline for perturbation responses
        YValues2(i,:) = [mean(PopResponses(:,min(find(TF==0)):max(find(TF==0))),2)- ...
                    mean(PopResponses(:,1:trialSide_samples),2)]';
        EnsembleSize(i) = size(Response,2);
end
                                                                    
figure
Data2Scatter = YValues2(find(EnsembleSize~=1),1)-YValues2(find(EnsembleSize~=1),:);
UnivarScatter(Data2Scatter(:,2:end),'Label',{'45','90','135','180','225','270','315'},'MarkerFaceColor',[ 0.5 0.5 0.5])
xlabel('Perturbation direction');
ylabel('Differential response');
set(gca,'TickDir','out','box','off')
title('Direction - running only')
ylim([-max(max(abs(YValues2(find(EnsembleSize~=1),1)-YValues2(find(EnsembleSize~=1),:)))) ...
       max(max(abs(YValues2(find(EnsembleSize~=1),1)-YValues2(find(EnsembleSize~=1),:))))]);

figure
YValuesOri =   [[YValues2(find(EnsembleSize~=1),1); YValues2(find(EnsembleSize~=1),5)], ...
                [YValues2(find(EnsembleSize~=1),2); YValues2(find(EnsembleSize~=1),6)], ...
                [YValues2(find(EnsembleSize~=1),3); YValues2(find(EnsembleSize~=1),7)], ...
                [YValues2(find(EnsembleSize~=1),4); YValues2(find(EnsembleSize~=1),8)]];
% boxplot(YValuesOri(:,1) - YValuesOri,'Notch','on','Labels',{'0','45','90','135'},'symbol','')
Data2Scatter = YValuesOri(:,1) - YValuesOri;
UnivarScatter(Data2Scatter(:,2:end),'Label',{'45','90','135'},'MarkerFaceColor',[ 0.5 0.5 0.5])

% hold on
% for i = 1:size(YValuesOri,2)
%     plot(i-0.1+0.2*rand(size(YValuesOri,1),1), YValuesOri(:,1)-YValuesOri(:,i), '.', 'Color',[0.5 0.5 0.5]);
% end
title('Orientation - running only')
xlabel('Perturbation orientation');
ylabel('Differential response');
set(gca,'TickDir','out','box','off')
hold on
ylim([-max(max(abs(YValuesOri(:,1) - YValuesOri))) max(max(abs(YValuesOri(:,1) - YValuesOri)))])
                                                                    
                                                                    
                                                                    
                                                                    
                                                                    
                                                                    
                                                                    
                                                                    
                                                                    
                                                                    
                                                                    
                                                                    
                                                                    
                                                                    
                                                                    
                                                                    
                                                                    
                                                                    
                                                                    
