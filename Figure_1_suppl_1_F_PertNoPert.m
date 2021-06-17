%% Tomaso Muzzu - UCL - 29/04/2020 - running modulation effects on grating stimuli

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
Param = squeeze(AM_Param(:,SelectedCells,3)); % pert ON/OFF for all units
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


%%
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
[SelectedTrials_OI SelectedTrials_OI_control TotTrials] = SpeedFiltering(AM_Speed, AM_SpeedResponse, AM_SpeedResponseControl, Trials_PertON,Trials_PertOFF,Run_TH, minTrials);
Condition = 2;
switch Condition
    case 1
        Param_run = SelectedTrials_OI{Condition,1}(:,SelectedCells);
        Param_still = SelectedTrials_OI{Condition,2}(:,SelectedCells);
    case 2
        Param_run = SelectedTrials_OI{Condition,1}(:,SelectedCells);
        Param_still = SelectedTrials_OI{Condition,2}(:,SelectedCells);
        Param_run_OFF = SelectedTrials_OI_control{Condition,1}(:,SelectedCells);
        Param_still_OFF = SelectedTrials_OI_control{Condition,2}(:,SelectedCells);
    case 3
        Param_run = SelectedTrials_OI{Condition,1}(:,SelectedCells);
        Param_still = SelectedTrials_OI{Condition,2}(:,SelectedCells);
    case 4
        Param_run = SelectedTrials_OI{Condition,1}(:,SelectedCells);
        Param_still = SelectedTrials_OI{Condition,2}(:,SelectedCells);
end

% go through unit by unit
% adjust the baseline so the two conditions are comparable
% compensate for the different number of trials, compare samples of
% similar sizes
% ranksum distributions of time bin-by-bin running vs still curves

r = 1; clear DM_rs DepthMod Responses ResponseON Peak_rs PopRes_run PopRes_still ranksum_p ranksum_h  AreaUC
clear ranksum_p_mean ranksum_h_mean ranksum_p_mean UnitResponseNP_g UnitResponsePT_g ttest_h_mean ttest_p_mean DM_RS_
clear UnitResponse_r_mean_ON UnitResponse_r_mean_OFF UnitResponse_s_mean_ON UnitResponse_s_mean_OFF UnitResponse_r_mean UnitResponse_s_mean
for i = 1:size(Param,2)
    % select perturbation trial responses
    Responses_ON = Param(:,i);
    % scrub responses
    Responses_ON = logical(Responses_ON(~isnan(Responses_ON)));
    %% save all trials for the two conditions
    clear UnitResponse_r UnitResponse_s
    % running trials PERTURBATION ON
    Responses = Responses_ON & Param_run(1:length(Responses_ON),i);
    UnitResponse_run_ON = squeeze(UnitResponses_smooth(Responses==1,i,:))*60; 
    % running trials PERTURBATION OFF
    Responses = ~Responses_ON & Param_run_OFF(1:length(Responses_ON),i);
    UnitResponse_run_OFF = squeeze(UnitResponses_smooth(Responses==1,i,:))*60; 
    
    % still trials PERTURBATION ON
    Responses = Responses_ON & Param_still(1:length(Responses_ON),i);
    UnitResponse_still_ON = squeeze(UnitResponses_smooth(Responses==1,i,:))*60;
    % still trials PERTURBATION OFF
    Responses = ~Responses_ON & Param_still_OFF(1:length(Responses_ON),i);
    UnitResponse_still_OFF = squeeze(UnitResponses_smooth(Responses==1,i,:))*60;
    
    % All perturbation trials
    Responses = Responses_ON;
    UnitResponsePT_g(i,:) = nanmean(squeeze(UnitResponses_smooth(Responses==1,i,:))*60);
    % normal trials
    Responses_OFF = Param(:,i)==0;
    % scrub responses
    Responses_OFF = Responses_OFF(1:length(Responses_ON));
    Responses = Responses_OFF;
    UnitResponseNP_g(i,:) = nanmean(squeeze(UnitResponses_smooth(Responses==1,i,:))*60);
    
    % find the smaller sample size
    [SampleSize SampleId] = min([size(UnitResponse_run_ON,1) size(UnitResponse_still_ON,1)]);
    [SampleSize1 SampleId] = min([size(UnitResponse_run_ON,1) size(UnitResponse_run_OFF,1)]);
    [SampleSize2 SampleId] = min([size(UnitResponse_still_ON,1) size(UnitResponse_still_OFF,1)]);
    
    %% select equal number of trials for each condition
    for k = 1:1000
        %% re-sample trials with replacement - compute mean curves with equal number or trials for run_ON and still_ON
        if size(UnitResponse_run_ON,1)==0 |  size(UnitResponse_still_ON,1)==0
            trials_run = 0; trials_still = 0;
            UnitResponse_r_mean_t = nanmean(UnitResponse_run_ON(1:SampleSize,po_i(1)-trialSide_samples:po_i(1)+trialSide_samples));
            UnitResponse_s_mean_t = nanmean(UnitResponse_still_ON(1:SampleSize,po_i(1)-trialSide_samples:po_i(1)+trialSide_samples));
        else
            trials_run = randi([1 size(UnitResponse_run_ON,1)],1,20);
            trials_still = randi([1 size(UnitResponse_still_ON,1)],1,20);
            UnitResponse_r_mean_t = nanmean(UnitResponse_run_ON(trials_run,po_i(1)-trialSide_samples:po_i(1)+trialSide_samples));
            UnitResponse_s_mean_t = nanmean(UnitResponse_still_ON(trials_still,po_i(1)-trialSide_samples:po_i(1)+trialSide_samples));
        end
        %% re-sample trials with replacement - compute mean curves with equal number or trials for run_ON and run_OFF
        if size(UnitResponse_run_ON,1)==0 |  size(UnitResponse_run_OFF,1)==0
            trials_run_ON = 0; trials_run_OFF = 0;
            UnitResponse_r_mean_t_ON = nanmean(UnitResponse_run_ON(1:SampleSize1,po_i(1)-trialSide_samples:po_i(1)+trialSide_samples));
            UnitResponse_r_mean_t_OFF = nanmean(UnitResponse_run_OFF(1:SampleSize1,po_i(1)-trialSide_samples:po_i(1)+trialSide_samples));
        else
            trials_run_ON = randi([1 size(UnitResponse_run_ON,1)],1,20);
            trials_run_OFF = randi([1 size(UnitResponse_run_OFF,1)],1,20);
            UnitResponse_r_mean_t_ON = nanmean(UnitResponse_run_ON(trials_run_ON,po_i(1)-trialSide_samples:po_i(1)+trialSide_samples));
            UnitResponse_r_mean_t_OFF = nanmean(UnitResponse_run_OFF(trials_run_OFF,po_i(1)-trialSide_samples:po_i(1)+trialSide_samples));
        end 
        %% re-sample trials with replacement - compute mean curves with equal number or trials for still_ON and still_OFF
        if size(UnitResponse_still_ON,1)==0 |  size(UnitResponse_still_OFF,1)==0
            trials_still_ON = 0; trials_still_OFF = 0;
            UnitResponse_s_mean_t_ON = nanmean(UnitResponse_still_ON(1:SampleSize2,po_i(1)-trialSide_samples:po_i(1)+trialSide_samples));
            UnitResponse_s_mean_t_OFF = nanmean(UnitResponse_still_OFF(1:SampleSize2,po_i(1)-trialSide_samples:po_i(1)+trialSide_samples));
        else
            trials_still_ON = randi([1 size(UnitResponse_still_ON,1)],1,20);
            trials_still_OFF = randi([1 size(UnitResponse_still_OFF,1)],1,20);
            UnitResponse_s_mean_t_ON = nanmean(UnitResponse_still_ON(trials_still_ON,po_i(1)-trialSide_samples:po_i(1)+trialSide_samples));
            UnitResponse_s_mean_t_OFF = nanmean(UnitResponse_still_OFF(trials_still_OFF,po_i(1)-trialSide_samples:po_i(1)+trialSide_samples));
        end 
        %% re-sample trials without replacement
%         UnitResponse_r = UnitResponse_run(randperm(size(UnitResponse_run,1)),:);
%         UnitResponse_s = UnitResponse_still(randperm(size(UnitResponse_still,1)),:);
%         % compute the mean and subtract the baseline
%         UnitResponse_r_mean = nanmean(UnitResponse_r(1:SampleSize,po_i(1)-trialSide_samples:po_i(1)+trialSide_samples));
%         UnitResponse_s_mean = nanmean(UnitResponse_s(1:SampleSize,po_i(1)-trialSide_samples:po_i(1)+trialSide_samples));
        %%
        %UnitResponse_r_mean(k,:) = UnitResponse_r_mean-mean(UnitResponse_r_mean(1:60));
        %UnitResponse_s_mean(k,:) = UnitResponse_s_mean-mean(UnitResponse_s_mean(1:60));
        UnitResponse_r_mean(k,:) = UnitResponse_r_mean_t;%-mean(UnitResponse_r_mean_t(1:60));
        UnitResponse_s_mean(k,:) = UnitResponse_s_mean_t;%-mean(UnitResponse_s_mean_t(1:60));
        
        UnitResponse_r_mean_ON(k,:) = UnitResponse_r_mean_t_ON;
        UnitResponse_r_mean_OFF(k,:) = UnitResponse_r_mean_t_OFF;
        
        UnitResponse_s_mean_ON(k,:) = UnitResponse_s_mean_t_ON;
        UnitResponse_s_mean_OFF(k,:) = UnitResponse_s_mean_t_OFF;
        
        %% take 50ms bins and compare (ranksum) the two mean curves for the entire pert period
        try
            [ranksum_p(i,k) ranksum_h(i,k)] = ranksum(movmean(UnitResponse_r_mean(k, 61:end),3), movmean(UnitResponse_s_mean(k, 61:end),3))  ;
        catch
            ranksum_p(i,k) = nan;
            ranksum_h(i,k) = 2;
        end
        %% compute modulation index for curve
        DM_RS_(i,k,1) = (nanmean(UnitResponse_r_mean(k,61:end))-nanmean(UnitResponse_r_mean(k,1:60)))/...
                        nanmean(UnitResponse_r_mean(k,1:60));
        DM_RS_(i,k,2) = (nanmean(UnitResponse_s_mean(k,61:end))-nanmean(UnitResponse_s_mean(k,1:60)))/...
                        nanmean(UnitResponse_s_mean(k,1:60));
                    
        AreaUC(i,k,1) = sum(UnitResponse_r_mean(k,61:end)); %-mean(UnitResponse_r_mean_t(1:60))); % running pert ON
        AreaUC(i,k,2) = sum(UnitResponse_s_mean(k,61:end)); %-mean(UnitResponse_s_mean_t(1:60))); % still pert ON
        
        AreaUC(i,k,3) = sum(UnitResponse_r_mean_ON(k,61:end)); %-mean(UnitResponse_r_mean_t_ON(1:60))); % running pert ON
        AreaUC(i,k,4) = sum(UnitResponse_r_mean_OFF(k,61:end)); %-mean(UnitResponse_r_mean_t_OFF(1:60))); % running pert OFF
        
        AreaUC(i,k,5) = sum(UnitResponse_s_mean_ON(k,61:end)); %-mean(UnitResponse_s_mean_t_ON(1:60))); % still pert ON
        AreaUC(i,k,6) = sum(UnitResponse_s_mean_OFF(k,61:end)); %-mean(UnitResponse_s_mean_t_OFF(1:60))); % still pert OFF
    end
    try
        [ranksum_p_mean(i) ranksum_h_mean(i)] = ranksum(movmean(mean(UnitResponse_r_mean(:, 61:end)),3), movmean(mean(UnitResponse_s_mean(:, 61:end)),3))  ;
        [ttest_h_mean(i) ttest_p_mean(i)] = ttest2(movmean(mean(UnitResponse_r_mean(:, 61:end)),3), movmean(mean(UnitResponse_s_mean(:, 61:end)),3))  ;
    catch
        ranksum_p_mean(i) = 2;
        ranksum_h_mean(i) = 2;
        ttest_h_mean(i) = 2 ;
        ttest_p_mean(i) = 2;
    end
    
    % running trials 
    Responses = Responses_ON & Param_run(1:length(Responses_ON),i);
    UnitResponse = nanmean(squeeze(UnitResponses_smooth(Responses==1,i,:)))*60;
    % UnitResponse = UnitResponse - nanmean(UnitResponse(po_i(1)-trialSide_samples:po_i(1)));
    DM_rs(i,1) = (nanmean(UnitResponse(po_i(1):po_i(1)+trialSide_samples))-nanmean(UnitResponse(po_i(1)-trialSide_samples:po_i(1))))/...
               nanmean(UnitResponse(po_i(1)-trialSide_samples:po_i(1)));
    Peak_rs(i,1) = nanmax(UnitResponse(po_i(1):po_i(1)+trialSide_samples));            
    PopRes_run(i,:) = UnitResponse;
      
    % still trials
    Responses = Responses_ON & Param_still(1:length(Responses_ON),i);
    UnitResponse = mean(squeeze(UnitResponses_smooth(Responses==1,i,:)))*60;
    % UnitResponse = UnitResponse - nanmean(UnitResponse(po_i(1)-trialSide_samples:po_i(1)));
    DM_rs(i,2) = (nanmean(UnitResponse(po_i(1):po_i(1)+trialSide_samples))-nanmean(UnitResponse(po_i(1)-trialSide_samples:po_i(1))))/...
                nanmean(UnitResponse(po_i(1)-trialSide_samples:po_i(1)));
    Peak_rs(i,2) = nanmax(UnitResponse(po_i(1):po_i(1)+trialSide_samples));
    PopRes_still(i,:) = UnitResponse;

end


%% do the neurons respond to perturbation 
clear DM_rs_i DM_rs_Vi
DM_rs_i = DM_rs(PertRespUnits_pos,:);
Excl_i = ~isnan(DM_rs_i(:,1)) & ~isnan(DM_rs_i(:,2)) & (DM_rs_i(:,1)~=0 & DM_rs_i(:,2)~=0);
AreaUC_i = AreaUC(PertRespUnits_pos,:,:);
AreaUC_i = AreaUC_i(Excl_i,:,:);


% make the three plots indicating the units sign. modulated by running
% show the p for the population subset 
% adjust the plot for the one siggested by aman

PopRes_run_i = PopRes_run(PertRespUnits_pos,:);
PopRes_still_i = PopRes_still(PertRespUnits_pos,:);
PopRes_run_i = PopRes_run_i(Excl_i,:);
PopRes_still_i = PopRes_still_i(Excl_i,:);
jj = 1; clear index_of_i
for j = 1:size(PopRes_run_i,1)
    Peak(j,1) = max( PopRes_run_i(j,po_i(1):po_i(1)+trialSide_samples) - mean(PopRes_run_i(j,po_i(1)-trialSide_samples:po_i(1))) );
    Peak(j,2) = max( PopRes_still_i(j,po_i(1):po_i(1)+trialSide_samples) - mean(PopRes_still_i(j,po_i(1)-trialSide_samples:po_i(1))) );
    
    if ( (sum(AreaUC_i(j,:,1)>AreaUC_i(j,:,2))>950) & Peak(j,1)>Peak(j,2) )
        plot(linspace(-1,1,length(PopRes_run_i(j,po_i(1)-trialSide_samples:po_i(1)+trialSide_samples))),...
            PopRes_run_i(j,po_i(1)-trialSide_samples:po_i(1)+trialSide_samples) - mean(PopRes_run_i(j,po_i(1)-trialSide_samples:po_i(1))) ...
            ,'r','LineWidth',3)
        index_of_i(jj) = j;
        jj = jj + 1;
    end
    
end

figure
subplot(1,3,1)
plot(mean(AreaUC_i(:,:,2),2),mean(AreaUC_i(:,:,1),2),'.r','MarkerSize',8)
hold on
plot(mean(AreaUC_i(index_of_i,:,2),2),mean(AreaUC_i(index_of_i,:,1),2),'or','MarkerSize',5)
set(gca,'TickDir','out','FontSize',12,'FontName','Arial')
grid on
box off
hold on
plot([-500 500], [-500 500], ':k')
xlabel('AUC_s_t_i_l_l_-_O_N'); ylabel('AUC_r_u_n_-_O_N');
title(['AUC_r_u_n_-_O_N VS AUC_s_t_i_l_l_-_O_N n = ' num2str(size(AreaUC_i,1))])
[p h] = signrank(mean(AreaUC_i(:,:,2),2)-mean(AreaUC_i(:,:,1),2))
text(500,30,['p=' num2str(p)])
% [p h] = signrank(AreaUC_i(index_of_i,1,2)-AreaUC_i(index_of_i,1,1))
% text(500,-300,['p=' num2str(p)])
set(gca,'xTick',[-1000:500:1000],'yTick',[-1000:500:1000])
xlim([-1000 1000]); ylim([-1000 1000])

subplot(1,3,2)
plot(mean(AreaUC_i(:,:,4),2),mean(AreaUC_i(:,1,3),2),'.r','MarkerSize',8)
hold on
%plot(mean(AreaUC_i(index_of_i,:,4),2),mean(AreaUC_i(index_of_i,:,3),2),'*r','MarkerSize',5)
set(gca,'TickDir','out','FontSize',12,'FontName','Arial')
grid on
box off
hold on
plot([-500 500], [-500 500], ':k')
xlabel('AUC_r_u_n_-_O_F_F'); ylabel('AUC_r_u_n_-_O_N');
title(['AUC_r_u_n_-_O_N VS AUC_r_u_n_-_O_F_F n = ' num2str(size(AreaUC_i,1))])
[p h] = signrank(mean(AreaUC_i(:,:,3),2)-mean(AreaUC_i(:,:,4),2))
text(250,30,['p=' num2str(p)])
set(gca,'xTick',[-1000:500:1000],'yTick',[-1000:500:1000])
xlim([-1000 1000]); ylim([-1000 1000])

subplot(1,3,3)
plot(mean(AreaUC_i(:,:,6),2),mean(AreaUC_i(:,:,5),2),'.r','MarkerSize',8)
hold on
%plot(mean(AreaUC_i(index_of_i,:,6),2),mean(AreaUC_i(index_of_i,:,5),2),'*r','MarkerSize',5)
set(gca,'TickDir','out','FontSize',12,'FontName','Arial')
grid on
box off
hold on
plot([-500 500], [-500 500], ':k')
xlabel('AUC_s_t_i_l_l_-_O_F_F'); ylabel('AUC_s_t_i_l_l_-_O_N');
title(['AUC_s_t_i_l_l_-_O_N VS AUC_s_t_i_l_l_-_O_F_F n = ' num2str(size(AreaUC_i,1))])
[p h] = signrank(mean(AreaUC_i(:,:,6),2)-mean(AreaUC_i(:,:,5),2))
text(250,30,['p=' num2str(p)])
set(gca,'xTick',[-1000:500:1000],'yTick',[-1000:500:1000])
xlim([-1000 1000]); ylim([-1000 1000])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
plot(mean(AreaUC_i(:,:,2),2)-mean(AreaUC_i(:,:,6),2),mean(AreaUC_i(:,:,1),2)-mean(AreaUC_i(:,:,4),2),'.r','MarkerSize',8)
hold on
plot(mean(AreaUC_i(index_of_i,:,2),2)-mean(AreaUC_i(index_of_i,:,6),2),mean(AreaUC_i(index_of_i,:,1),2)-mean(AreaUC_i(index_of_i,:,4),2),'or','MarkerSize',5)
set(gca,'xTick',[-1000:500:1000],'yTick',[-1000:500:1000])
xlim([-500 1000]); ylim([-500 1000])
set(gca,'TickDir','out','FontSize',12,'FontName','Arial')
grid on
box off
hold on
plot([-1000 1000], [-1000 1000], ':k')
plot([-500 500], [-1000 1000], ':k')
[p h] = signrank((mean(AreaUC_i(index_of_i,:,2),2)-mean(AreaUC_i(index_of_i,:,6),2))-(mean(AreaUC_i(index_of_i,:,1),2)-mean(AreaUC_i(index_of_i,:,4),2)))
text(250,30,['p=' num2str(p)])
xlabel('AUC_s_t_i_l_l_-_O_N - AUC_s_t_i_l_l_-_O_F_F'); ylabel('AUC_r_u_n_-_O_N - AUC_r_u_n_-_O_F_F');
title(['n = ' num2str(size(AreaUC_i(index_of_i,:,:),1))])

[p h] = signrank(mean(AreaUC_i(index_of_i,:,2),2)-mean(AreaUC_i(index_of_i,:,6),2),mean(AreaUC_i(index_of_i,:,1),2)-mean(AreaUC_i(index_of_i,:,4),2))

length(mean(AreaUC_i(index_of_i,:,2),2)-mean(AreaUC_i(index_of_i,:,6),2))

sum((mean(AreaUC_i(:,:,1),2)-mean(AreaUC_i(:,:,4),2)) > 2*(mean(AreaUC_i(:,:,2),2)-mean(AreaUC_i(:,:,6),2)))

sum((mean(AreaUC_i(index_of_i,:,1),2)-mean(AreaUC_i(index_of_i,:,4),2)) > 2*(mean(AreaUC_i(index_of_i,:,2),2)-mean(AreaUC_i(index_of_i,:,6),2)))

