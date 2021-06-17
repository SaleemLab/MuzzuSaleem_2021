%% Tomaso Muzzu - UCL - 08/06/2020
% slope of MI's of firing rate computed at each trial

%% Functions to plot Figure 1 suppl 1 - running modulation
if ~exist('ProjectData','var')
    [ProjectData AM_UnitResponses AM_Param AM_Speed AM_UOI SelectedResponses AM_UnitResponses_smooth] = LoadDataALL;
end

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


%% step 1: compute the % difference in response during the perturbation
% period
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
minTrials = 4; % min nr of trials for each condition
[SelectedTrials_OI SelectedTrials_OI_control TotTrials] = SpeedFiltering(AM_Speed, AM_SpeedResponse, AM_SpeedResponseControl, Trials_PertON,Trials_PertOFF,...
                                                                        Run_TH, minTrials);
% SelectedTrials_OI{4,2} = 4 rows as above, column 1 = pert ON & run; column 2 = pert ON & static
% SelectedTrials_OI_control{4,2} = 4 rows as above, column 1 = pert OFF & run; column 2 = pert OFF & static

%% plot modulation indexes
% consider only running during perturbation period
Param = squeeze(AM_Param(:,SelectedCells,3)); % all units
Param_run = SelectedTrials_OI{2,1}(:,SelectedCells); 
Param_still = SelectedTrials_OI{2,2}(:,SelectedCells); 

% Param = squeeze(AM_Param(:,RespUnits(:,1),3)); % only units responsive to pert stim 
% Param_run = SelectedTrials_OI{2,1}(:,RespUnits(:,1)); 
% Param_still = SelectedTrials_OI{2,2}(:,RespUnits(:,1)); 
% Depth of modulation
r = 1; clear DM_rs DepthMod Responses ResponseON 
for i = 1:size(Param,2)
    % select trial responses
    Responses_ON = Param(:,i);
    % scrub responses
    Responses_ON = logical(Responses_ON(~isnan(Responses_ON)));
    
    % running trials 
    Responses = Responses_ON & Param_run(1:length(Responses_ON),i);
    UnitResponse = nanmean(squeeze(UnitResponses_smooth(Responses==1,i,:)))*60; 
    Dir = sign(sum(UnitResponse(po_i(1):po_i(1)+trialSide_samples))-sum(UnitResponse(po_i(1)-trialSide_samples:po_i(1))));
     DM_rs(i,1) = (nanmean(UnitResponse(po_i(1):po_i(1)+trialSide_samples))-nanmean(UnitResponse(po_i(1)-trialSide_samples:po_i(1))))/...
                (nanmean(UnitResponse(po_i(1)-trialSide_samples:po_i(1))));
            index_i(i) = i;
%     switch Dir
%         case 1
%             DM_rs(i,1) = (nanmax(UnitResponse(po_i(1):po_i(1)+trialSide_samples))-nanmax(UnitResponse(po_i(1)-trialSide_samples:po_i(1))))/...
%                 (nanmax(UnitResponse(po_i(1):po_i(1)+trialSide_samples))+nanmax(UnitResponse(po_i(1)-trialSide_samples:po_i(1))));
%             index_i(i) = i;
%         case -1
%             DM_rs(i,1) = (min(UnitResponse(po_i(1):po_i(1)+trialSide_samples))-min(UnitResponse(po_i(1)-trialSide_samples:po_i(1))))/...
%                 (min(UnitResponse(po_i(1):po_i(1)+trialSide_samples))+min(UnitResponse(po_i(1)-trialSide_samples:po_i(1))));
%             index_i(i) = i;
%     end 
     
    % still trials
    Responses = Responses_ON & Param_still(1:length(Responses_ON),i);
    UnitResponse = mean(squeeze(UnitResponses_smooth(Responses==1,i,:)))*60;
    Dir = sign(sum(UnitResponse(po_i(1):po_i(1)+trialSide_samples))-sum(UnitResponse(po_i(1)-trialSide_samples:po_i(1))));
    DM_rs(i,2) = (nanmean(UnitResponse(po_i(1):po_i(1)+trialSide_samples))-nanmean(UnitResponse(po_i(1)-trialSide_samples:po_i(1))))/...
                (nanmean(UnitResponse(po_i(1)-trialSide_samples:po_i(1))));
            index_i(i) = i;
%     switch Dir
%         case 1
%             DM_rs(i,2) = (nanmax(UnitResponse(po_i(1):po_i(1)+trialSide_samples))-nanmax(UnitResponse(po_i(1)-trialSide_samples:po_i(1))))/...
%                 (nanmax(UnitResponse(po_i(1):po_i(1)+trialSide_samples))+nanmax(UnitResponse(po_i(1)-trialSide_samples:po_i(1))));
%             index_i(i) = i;
%         case -1
%             DM_rs(i,2) = (min(UnitResponse(po_i(1):po_i(1)+trialSide_samples))-min(UnitResponse(po_i(1)-trialSide_samples:po_i(1))))/...
%                 (min(UnitResponse(po_i(1):po_i(1)+trialSide_samples))+min(UnitResponse(po_i(1)-trialSide_samples:po_i(1))));
%             index_i(i) = i;
%     end
     
end
DM_rs(size(DM_rs,1)+1:size(Param_run,2),:) = zeros(size(Param_run,2)-size(DM_rs,1),2);

%% look at the correlation between firing rate and speed
% normalise firing rate
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

%% look at trial by trial response to see if there is a correlation whatsoever
RespUnits_idx = find(PertRespUnits_pos); % indexes of pert. resp. units
AM_ParamSel = AM_Param(:,AM_UOI,:); % select parameters for reliable units
AM_SpeedSel = AM_Speed(:,AM_UOI,1:559); % select speeds for reliable units
Rec_OI = unique(AM_ParamSel(1,PertRespUnits_pos,1)); % recordings of interest (there were pert. resp. units)
RecIndex = unique(AM_ParamSel(1,PertRespUnits_pos,1)); % indexes of recordings of interest
UnitsPertPos = PertRespUnits_pos;
BonvisionFR = 60; %Hz  
trialSide_samples = 60; trialSide_seconds = 1; 
PertStartEnds = [round(nanmin(AM_ParamSel(:,RecIndex,4),[],1)*BonvisionFR+trialSide_samples) ;...
                round(nanmax(AM_ParamSel(:,RecIndex,5),[],1)*BonvisionFR+trialSide_samples)];
k = 1;
clear MI_neuro MI_corr Speed_trial rho MI_speed
for i = 1:length(RecIndex)
    Unit_RoI = find(AM_ParamSel(1,UnitsPertPos,1)==RecIndex(i)); %AM_ParamSel(1,RecIndex(i),1));
    Unit_idx = RespUnits_idx(Unit_RoI);
    PertTrials_1 = AM_ParamSel(:,RecIndex(i),3)==1;
    PertTrials_0 = AM_ParamSel(:,RecIndex(i),3)==0;
    for j = 1:length(Unit_idx)
        clear UnitResponseData SpeedResponseData
        UnitResponseData = squeeze(UnitResponses_n(PertTrials_1,Unit_idx(j),PertStartEnds(1,i)-60:PertStartEnds(2,i)));
        SpeedResponseData = squeeze(SpeedResponse_n(PertTrials_1,Unit_idx(j),PertStartEnds(1,i)-60:PertStartEnds(2,i)));
        SpeedResponseData_real = squeeze(AM_SpeedSel(PertTrials_1,Unit_idx(j),PertStartEnds(1,i)-60:PertStartEnds(2,i)));
        rho_ = diag(corr(UnitResponseData',SpeedResponseData'));
        rho{k} = rho_(~isnan(rho_));
        % MI for speed and neural response
        Speed_trial{k} = mean(SpeedResponseData_real(:,1:60),2);
        MI_neuro = (nanmean(UnitResponseData(:,60:end),2)-nanmean(UnitResponseData(:,1:60),2))./...
                   (nanmean(UnitResponseData(:,1:60),2));
        MI_speed = (nanmean(SpeedResponseData(:,60:end),2)-nanmean(SpeedResponseData(:,1:60),2))./...
                   (nanmean(SpeedResponseData(:,1:60),2));
        MI_corr{k} = [MI_neuro, MI_speed];      
        k = k + 1;
    end
end

figure
suptitle('Correlation between speed and FR')
for k = 1:length(rho)
    subplot(9,10,k)
    histogram(abs(rho{k}),0:0.05:1)
end

figure
for k = 1:length(rho)
    CorrSpeedPert(k) = median(abs(rho{k}));
end
histogram(CorrSpeedPert,0:0.05:1)
ylabel('distribution of median R')
title(['Median correlation with pert period speed, n=' num2str(sum(UnitsPertPos))])
set(gca,'TickDir','out','box','off')
median(CorrSpeedPert)
std(CorrSpeedPert)
% look at correlation between mean firing rate response and speed trace...

%%

figure
%suptitle('MIneuro VS MIspeed')
clear p_MI_slope
for k = 1:90 %length(MI_corr)
    subplot(9,10,k)
    %plot(MI_corr{k}(:,2),MI_corr{k}(:,1),'.');
    hold on
    x = MI_corr{k}(:,2);
    y = MI_corr{k}(:,1);
    idx = (~isnan(x) & ~isinf(x) & ~isnan(y) & ~isinf(y));
    x = x(idx);    y = y(idx);
    fitobject = fit(x ,y ,'poly1');
    plot(fitobject,x,y)
    mdl = fitlm(x,y);
    p_MI_slope(k,:) = [mdl.Coefficients.Estimate(2) mdl.Coefficients.pValue(2)];   
    legend off
end

%% show examples
k = 10
figure
x = MI_corr{k}(:,2);
y = MI_corr{k}(:,1);
idx = (~isnan(x) & ~isinf(x) & ~isnan(y) & ~isinf(y));
x = x(idx);    y = y(idx);
fitobject = fit(x(~isnan(y)),y(~isnan(y)),'poly1');
plot(fitobject,x,y)
mdl = fitlm(x(~isnan(y)),y(~isnan(y)));
xlim([-1 1]); ylim([-1 1]);
set(gca,'TickDir','out','box','off')
xlabel('MI speed'); ylabel('MI firing rate')
title(['unit #' num2str(k)])
ll = legend({'trials' ,['y=' num2str(mdl.Coefficients.Estimate(1),1) '+' num2str(mdl.Coefficients.Estimate(2),1) 'x']},'fontsize',10);
ll.Color = 'none'; ll.EdgeColor = 'none';
text(-0.8, 0.8, ['p=' num2str(p_MI_slope(k,2)),2],'fontsize',10)
   
sum(p_MI_slope(:,2)<0.01)
find(p_MI_slope(:,2)<0.01)

figure
histogram(p_MI_slope(:,1),-1:0.1:1,'FaceColor',[0.5 0.5 0.5],'EdgeColor','none') 
hold on
histogram(p_MI_slope(p_MI_slope(:,2)<0.01,1),-1:0.1:1,'FaceColor','r','EdgeColor','none') 
ylabel('units');
xlabel('slope of linear fit');
set(gca,'TickDir','out','box','off')
ll = legend({['Pert. resp. units n=' num2str(length(MI_corr)-sum(p_MI_slope(:,2)<0.01))], ...
             ['Pert. resp. units with slope ~=0, n=' num2str(sum(p_MI_slope(:,2)<0.01))]},'fontsize',10);
ll.Color = 'none'; ll.EdgeColor = 'none';

%% check modulation wrt speed pre perturbation

figure
%suptitle('MIneuro VS speed pre-pert')
clear p_slope
for k = 1:size(Speed_trial,2)
   % subplot(12,10,k)
    % fit a linear plot
    x = Speed_trial{k};
    y = MI_corr{k}(:,1);
    %fitobject = fit(x(~isnan(y)),y(~isnan(y)),'poly1');
	%hold on
    %plot(fitobject,x,y)
    %legend off
    mdl = fitlm(x(~isnan(y)),y(~isnan(y)));
    p_slope(k,:) =  [mdl.Coefficients.Estimate(2) mdl.Coefficients.pValue(2)];
end

figure
plot(p_slope,'.')
find(p_slope<0.01)

%% show examples
figure
x = Speed_trial{k};
y = MI_corr{k}(:,1);
fitobject = fit(x(~isnan(y)),y(~isnan(y)),'poly1');
plot(fitobject,x,y)
mdl = fitlm(x(~isnan(y)),y(~isnan(y)));
%xlim([0 max(x)]); ylim([-1 1]);
set(gca,'TickDir','out','box','off')
xlabel('Speed pre-pert (cm/s)'); ylabel('MI firing rate')
title(['unit #' num2str(k)])
ll = legend({'trials' ,['y=' num2str(mdl.Coefficients.Estimate(1),1) '+' num2str(mdl.Coefficients.Estimate(2),1) 'x']},'fontsize',10);
ll.Color = 'none'; ll.EdgeColor = 'none';
text(4,0.3, ['p=' num2str(mdl.Coefficients.pValue(2),2)],'fontsize',10)
   
sum(p_slope<0.01)
find(p_slope<0.01)

figure
histogram(p_slope(:,1),-0.1:0.01:0.1,'FaceColor',[0.5 0.5 0.5],'EdgeColor','none') 
hold on
histogram(p_slope(p_slope(:,2)<0.01,1),-0.1:0.01:0.1,'FaceColor','r','EdgeColor','none') 
ylabel('units');
xlabel('1st coeff. of fit');
set(gca,'TickDir','out','box','off')




%% rate of change correlation 
RespUnits_idx=find(UnitsPertPos); k = 1; clear rho MI_corr MI_speed MI_neuro Speed_trial
for i = 1:length(RecIndex)
    Unit_RoI = find(AM_ParamSel(1,UnitsPertPos,1)==AM_ParamSel(1,RecIndex(i),1));
    Unit_idx = RespUnits_idx(Unit_RoI);
    PertTrials_1 = AM_ParamSel(:,RecIndex(i),3)==1;
    PertTrials_0 = AM_ParamSel(:,RecIndex(i),3)==0;
    for j = 1:length(Unit_idx)
        clear UnitResponseData SpeedResponseData
        UnitResponseData = diff(squeeze(UnitResponses_n(PertTrials_1,Unit_idx(j),PertStartEnds(1,i)-60:PertStartEnds(2,i))),5);
        SpeedResponseData = diff(squeeze(SpeedResponse_n(PertTrials_1,Unit_idx(j),PertStartEnds(1,i)-60:PertStartEnds(2,i))),5);
        SpeedResponseData_real = diff(squeeze(AM_SpeedSel(PertTrials_1,Unit_idx(j),PertStartEnds(1,i)-60:PertStartEnds(2,i))),5);
        % Pearson's correlation coefficient
        rho_ = diag(corr(UnitResponseData',SpeedResponseData'));
        rho{k} = rho_(~isnan(rho_));
        % MI for speed and neural response
        Speed_trial{k} = mean(SpeedResponseData_real(:,1:60),2);
        MI_neuro = (mean(UnitResponseData(:,60:end),2)-mean(UnitResponseData(:,1:60),2))./...
                   (mean(UnitResponseData(:,60:end),2)+mean(UnitResponseData(:,1:60),2));
        MI_speed = (mean(SpeedResponseData(:,60:end),2)-mean(SpeedResponseData(:,1:60),2))./...
                   (mean(SpeedResponseData(:,60:end),2)+mean(SpeedResponseData(:,1:60),2));
        MI_corr{k} = [MI_neuro, MI_speed];      
        k = k + 1;
    end
end

figure
suptitle('Correlation between changes of speed and FR')
for k = 1:length(rho)
    subplot(12,10,k)
    histogram(abs(rho{k}),0:0.05:1)
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
















