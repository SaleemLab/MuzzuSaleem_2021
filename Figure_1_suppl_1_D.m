%% Tomaso Muzzu - 09/03/2020 - UCL - Figure suppl. running/still classifier


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
thres = 95;
% first 7 animals
if  size(ProjectData,1)>10
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
    
else
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
end


%% fit linear regression model to see if behaviour can predict trial type
% normalise speed
AM_SpeedSel = AM_Speed(:,AM_UOI,1:559);
UnitResponses = AM_UnitResponses_smooth(:,AM_UOI,:);
clear UnitResponses_n
for t = 1:size(UnitResponses,2)
    OneUnitResponse = squeeze(UnitResponses(:,t,:));
    UnitResponses_n(:,t,:) = OneUnitResponse./nanmax(OneUnitResponse,[],2);
end
UnitResponses = UnitResponses_n;
% normalise speed
SpeedResponse_n = AM_SpeedSel;
for t = 1:length(SpeedResponse_n)
    OneSpeedResponse_n = squeeze(SpeedResponse_n(:,t,:));
    SpeedResponse_n(:,t,:) = OneSpeedResponse_n./nanmax(OneSpeedResponse_n,[],2);
end

% look at the speed profiles of each recording
AM_ParamSel = AM_Param(:,AM_UOI,:);
for i = 1:length(unique(AM_ParamSel(1,:,1)))
    AllRec_idx(i) = min(find(AM_ParamSel(1,:,1)==i));
end
SpeedResponse_nSel = SpeedResponse_n(:,AllRec_idx,:);
% find beginning and end of perturbation trials
% scan every recording and look at the mean speed responses during
% perturbation
Rec_OI = unique(AM_ParamSel(1,PertRespUnits_pos,1)); % recordings of interest

for i = 1:length(Rec_OI)
   RecIndex(i) = min(find((AM_ParamSel(1,:,1)==Rec_OI(i))))
end
clear MeanRecSpeed SemRecSpeed
for i = 1:length(RecIndex)
   % get perturbation trials
   PertTrials_1 = AM_ParamSel(:,RecIndex(i),3)==1; 
   PertTrials_0 = AM_ParamSel(:,RecIndex(i),3)==0;
    
   MeanRecSpeed(1,:,i) = nanmean(AM_SpeedSel(PertTrials_1,RecIndex(i),:));
   MeanRecSpeed(2,:,i) = nanmean(AM_SpeedSel(PertTrials_0,RecIndex(i),:));
   SemRecSpeed(1,:,i) = nanstd(AM_SpeedSel(PertTrials_0,RecIndex(i),:))/sqrt(sum(PertTrials_1));
   SemRecSpeed(2,:,i) = nanstd(AM_SpeedSel(PertTrials_0,RecIndex(i),:))/sqrt(sum(PertTrials_0));   
end

BonvisionFR = 60; %Hz  
trialSide_samples = 60; trialSide_seconds = 1; 
PertStartEnds = [round(nanmin(AM_ParamSel(:,RecIndex,4),[],1)*BonvisionFR+trialSide_samples) ;...
                round(nanmax(AM_ParamSel(:,RecIndex,5),[],1)*BonvisionFR+trialSide_samples)];
% extract predictors from speed signals
po_i = [min(PertStartEnds(1,:)) max(PertStartEnds(2,:))];
clear DiffResponsePercent_Pert
for i = 1:size(SpeedResponse_nSel,2)
    
    UnitResponse = squeeze(SpeedResponse_nSel(:,i,:));
    
    DiffResponsePercent_Pert(:,i,1) = (nanmean(UnitResponse(:,po_i(1):po_i(2)),2)./nanmean(UnitResponse(:,po_i(1)-trialSide_samples:po_i(1)),2)); % mean(pert)/mean(prepert)
    DiffResponsePercent_Pert(:,i,2) = (nanmean(UnitResponse(:,po_i(1):po_i(2)),2)-nanmean(UnitResponse(:,po_i(1)-trialSide_samples:po_i(1)),2)); % mean(pert)-mean(prepert)
    DiffResponsePercent_Pert(:,i,3) = (nansum(UnitResponse(:,po_i(1):po_i(2)),2)); % integral of mean(pert)
    DiffResponsePercent_Pert(:,i,4) = (nansum(UnitResponse(:,po_i(1)-trialSide_samples:po_i(2)),2)); % integral of mean(pert) and preceding 1s
    DiffResponsePercent_Pert(:,i,5) = (nanmean(UnitResponse(:,po_i(1):po_i(2)-trialSide_samples/2),2)-nanmean(UnitResponse(:,po_i(1)-trialSide_samples:po_i(1)),2))./...
                                      (nanmean(UnitResponse(:,po_i(1):po_i(2)-trialSide_samples/2),2)+nanmean(UnitResponse(:,po_i(1)-trialSide_samples:po_i(1)),2)); % mean(pert)-mean(prepert)/mean(pert+max(prepert) a.k.a depth of modulation
    DiffResponsePercent_Pert(:,i,6) = nanmean(UnitResponse(:,po_i(1):po_i(2)-trialSide_samples/2),2)./...
                                      (nanmean(UnitResponse(:,po_i(1):po_i(2)-trialSide_samples/2),2)+nanmean(UnitResponse(:,po_i(1)-trialSide_samples:po_i(1)),2)); % mean(pert)/mean(pert)+mean(prepert) a.k.a modulation index    
end

%% step 2: use these values with the actual response to fit a logistic model
Param = squeeze(AM_Param(:,AllRec_idx,3)); % save trial type 
clear SingleUnitResponse AUC_speed
plotted=0;
for i = 1:size(DiffResponsePercent_Pert,2) 
    % perturbation period
    SingleUnitResponse = squeeze(DiffResponsePercent_Pert(:,i,:));
    % select trial responses
    Responses = Param(:,i);
    % scrub responses
    Responses = logical(Responses(~isnan(Responses)));
    % scrub predictor
    SingleUnitResponse = SingleUnitResponse(1:length(Responses),:);
    SingleUnitResponse(SingleUnitResponse==Inf,:) = 0;
    SingleUnitResponse(isnan(SingleUnitResponse(:))) = 0;
    % apply logistic model
    mdl = fitglm(SingleUnitResponse,Responses,'Distribution','binomial','Link','logit');
    scores = mdl.Fitted.Probability;
    % step 3: compute the ROC and evaluate the AUC for each unit
    [X,Y,T,AUC1] = perfcurve(Responses,scores,1);
    AUC_speed(i,1) = AUC1;
end
tic
clear AUC_shuffled_Sp
parfor i = 1:size(DiffResponsePercent_Pert,2)
    SingleUnitResponse = squeeze(DiffResponsePercent_Pert(:,i,:));
    % select trial responses
    Responses = Param(:,i);
    % scrub responses
    Responses = logical(Responses(~isnan(Responses)));
    % scrub predictor
    SingleUnitResponse = SingleUnitResponse(1:length(Responses),:);
    SingleUnitResponse(SingleUnitResponse==Inf,:) = 0;
    SingleUnitResponse(isnan(SingleUnitResponse(:))) = 0;
    
    for sh_i = 1:1000
        % shuffle perturbation trals
        Responses_sh = Responses(randperm(length(Responses)));
        
        % apply logistic model
        mdl = fitglm(SingleUnitResponse,Responses_sh,'Distribution','binomial','Link','logit');
        scores = mdl.Fitted.Probability;
        % step 3: compute the ROC and evaluate the AUC for each unit
        [X,Y,T,AUC1] = perfcurve(Responses,scores,1);
        AUC_shuffled_Sp(i,sh_i) = AUC1;
    end
    i
end
toc
% AUC_shuffled_Sp = cat(2,AUC_speed(:,1),AUC_shuffled_Sp);
save('AUC_shuffled_Speed.mat','AUC_shuffled_Sp')


%%
if ~exist('AUC_shuffled_Sp','var')
    load('AUC_shuffled_Speed.mat')
    AUC_speed = AUC_shuffled_Sp(:,1);
    AUC_shuffled_Sp = AUC_shuffled_Sp(:,2:end);
end

p_sh = prctile(AUC_shuffled_Sp(:),thres);
find(AUC_speed>p_sh)
AM_ParamSel = AM_Param(:,AM_UOI,:);
RecIndex = unique(AM_ParamSel(1,PertRespUnits_pos,1));
AUC_2_plot = AUC_shuffled_Sp(RecIndex,:);
AUC_lin = AUC_2_plot(:,2:end);
p_sh = prctile(AUC_lin(:),thres);
%find(AUC_speed(AM_ParamSel(1,RecIndex,1))>p_sh)

%% show AUC of shuffled and original data
RecIndex = unique(AM_ParamSel(1,PertRespUnits_pos,1));
AUC_2_plot = AUC_shuffled_Sp(RecIndex,:);

figure
h = histogram(AUC_2_plot(:,2:end),[0:0.02:1],'FaceColor',[0.5 0.5 0.5],'Normalization','probability','EdgeColor','none')
hold on
h1 = histogram(AUC_speed,[0:0.02:1],'FaceColor',[1 0 0],'Normalization','probability','EdgeColor','none')
%plot([p_sh(1) p_sh(1)], [0 max(h.BinCounts/length(AUC_shuffled_Sp(:)))*1.1],'k') 
plot([p_sh p_sh], [0 max(h1.BinCounts/length(AUC_speed(AM_ParamSel(1,RecIndex,1))))*1.1],'k--','LineWidth',3)
xlabel('AUC of shuffled data'); ylabel('probability')
ll = legend({['AUC of shuffled data n=' num2str(size(AUC_2_plot,1)) 'k'],['AUC of sessions with resp units n=' num2str(size(AUC_2_plot,1))], '99% shuffled data'})
ll.Color = 'none'; ll.EdgeColor = 'none';
box off
set(gca,'FontSize',13,'TickDir','out');

sum(AUC_speed(RecIndex)>p_sh)

% how many units come from these 6 recordings?
R_psh = RecIndex((AUC_speed(RecIndex)>p_sh));
for i =1:length(R_psh)
    nrUnitsXrec(i) = sum(AM_ParamSel(1,PertRespUnits_pos,1)==R_psh(i))
end


%% compute the speed score (Pearson's correlation) between units and speed
% smoothing parameters
Width = round(0.3*BonvisionFR);
Sigma = Width/3; % standard deviation in number of samples (converted from time in seconds)
x_g = linspace(-Width/2, Width/2, Width);
gaussFilter = exp(-x_g.^2/(2*Sigma^2));
gaussFilter_ = gaussFilter / sum (gaussFilter); % normalize
k = 1;
clear FR_session CorrCoeffSpeed CorrCoeffSpeed_sh
for i = 1:size(ProjectData,1)
    SessionTime = ProjectData.Session_data{i,1}.Time{1,1}; 
    MouseSpeed = conv(ProjectData.Session_data{i,1}.MouseSpeed{1,1},gaussFilter_,'same');
    for j = 1:size(ProjectData.Units_Info{i,1},1)
        FR_session{i}(j,:) = conv(histcounts(ProjectData.Units_Info{i,1}.Spiketimes{j,1},SessionTime),gaussFilter_,'same');
        Rho = corrcoef(FR_session{i}(j,:),MouseSpeed);
        CorrCoeffSpeed(k) = Rho(1,2);
%         for ll = 1:1000
%             Lag = rand*20;
%             Sh_SC = ProjectData.Units_Info{i,1}.Spiketimes{j,1}+Lag;
%             Sh_SC(Sh_SC>max(SessionTime)) = Sh_SC(Sh_SC>max(SessionTime))-max(SessionTime);
%             Sh_SC = sort(Sh_SC);
%             FR_session_sh = conv(histcounts(Sh_SC,SessionTime),gaussFilter_,'same');
%             Rho = corrcoef(FR_session_sh,MouseSpeed);
%             CorrCoeffSpeed_sh(k,ll) = Rho(1,2);
%         end
        k = k + 1;
    end
end
save('CorrCoeffSpeed_sh.mat','CorrCoeffSpeed_sh');
% shuffle the firing rate and recompute the corr coeff. 1000 times if
% select only the units worth analysing
CorrCoeffSpeed_sh = CorrCoeffSpeed_sh(AM_UOI,:);

thres_1_99 = [prctile(CorrCoeffSpeed_sh(:),1) prctile(CorrCoeffSpeed_sh(:),99)];

for u = 1:length(CorrCoeffSpeed)
    p_speed(u) = sum(CorrCoeffSpeed(u)>thres_99(u))
end
CorrCoeffSpeed = CorrCoeffSpeed(AM_UOI);
figure
h = histogram(CorrCoeffSpeed_sh,[-1:0.02:1],'FaceColor',[0.5 0.5 0.5],'Normalization','probability','EdgeColor','none')
hold on
plot([thres_1_99(1) thres_1_99(1)],[0 max(h.BinCounts)/length(CorrCoeffSpeed_sh(:))], 'k:');
plot([thres_1_99(2) thres_1_99(2)],[0 max(h.BinCounts)/length(CorrCoeffSpeed_sh(:))], 'k:');
h1 = histogram(CorrCoeffSpeed,[-1:0.02:1],'FaceColor',[1 0 0],'Normalization','probability','EdgeColor','none')
xlabel('Speed Neuro Corr.'); ylabel('probability')
box off
set(gca,'TickDir','out')
title('Speed score')
% count how many units are modulated by running
sum(CorrCoeffSpeed>thres_1_99(2))
sum(CorrCoeffSpeed<thres_1_99(1))
Units_run = [find(CorrCoeffSpeed>thres_1_99(2)) find(CorrCoeffSpeed<thres_1_99(1))];
gg = 1; clear idx
for i = 1:length(Units_run)
    if ~isempty(find(Units_run(i)==RespUnits_idx))
        idx(gg) = find(Units_run(i)==RespUnits_idx);
        gg = gg + 1;
    end
end
% recordings that contain units modulated by speed
AM_ParamSel(1,RespUnits_idx(idx),1)

AM_ParamSel(1,RecIndex,1)


