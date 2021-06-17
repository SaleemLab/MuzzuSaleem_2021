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


%% modulation index of perturbation stimulus responses
% consider only running during visual stimulus onset period (-0.5s +2s)
Param = squeeze(AM_Param(:,SelectedCells,3)); % pert ON/OFF for all units
Condition = 2;
switch Condition
    case 1
        Param_run = SelectedTrials_OI{Condition,1}(:,SelectedCells);
        Param_still = SelectedTrials_OI{Condition,2}(:,SelectedCells);
    case 2
        Param_run = SelectedTrials_OI{Condition,1}(:,SelectedCells);
        Param_still = SelectedTrials_OI{Condition,2}(:,SelectedCells);
    case 3
        Param_run = SelectedTrials_OI{Condition,1}(:,SelectedCells);
        Param_still = SelectedTrials_OI{Condition,2}(:,SelectedCells);
    case 4
        Param_run = SelectedTrials_OI{Condition,1}(:,SelectedCells);
        Param_still = SelectedTrials_OI{Condition,2}(:,SelectedCells);
end
% Modulation Index - MI - perturbation responses
r = 1; clear DM_rs DepthMod Responses ResponseON Peak_rs PopRes_run PopRes_still
for i = 1:size(Param,2)
    % select trial responses
    Responses_ON = Param(:,i);
    % scrub responses
    Responses_ON = logical(Responses_ON(~isnan(Responses_ON)));
    
    % running trials 
    Responses = Responses_ON & Param_run(1:length(Responses_ON),i);
    UnitResponse = nanmean(squeeze(UnitResponses_smooth(Responses==1,i,:)))*60;
    %UnitResponse = UnitResponse - mean(UnitResponse(po_i(1)-trialSide_samples:po_i(1)));
    Dir = sign(sum(UnitResponse(po_i(1):po_i(1)+trialSide_samples))-sum(UnitResponse(po_i(1)-trialSide_samples:po_i(1))));
    
    DM_rs(i,1) = (nanmax(UnitResponse(po_i(1):po_i(1)+trialSide_samples))-nanmax(UnitResponse(po_i(1)-trialSide_samples:po_i(1))))/...
                 nanmax(UnitResponse(po_i(1)-trialSide_samples:po_i(1)));
%     DM_rs(i,1) = (nanmean(UnitResponse(po_i(1):po_i(1)+trialSide_samples))-nanmean(UnitResponse(po_i(1)-trialSide_samples:po_i(1))))/...
%                  nanmean(UnitResponse(po_i(1)-trialSide_samples:po_i(1)));
    
    Peak_rs(i,1) = nanmax(UnitResponse(po_i(1):po_i(1)+trialSide_samples));            
    PopRes_run(i,:) = UnitResponse;
%     switch Dir
%         case 1
%             DM_rs(i,1) = (max(UnitResponse(po_i(1):po_i(1)+trialSide_samples))-max(UnitResponse(po_i(1)-trialSide_samples:po_i(1))))/...
%                 (max(UnitResponse(po_i(1):po_i(1)+trialSide_samples))+max(UnitResponse(po_i(1)-trialSide_samples:po_i(1))));
%             index_i(i) = i;
%         case -1
%             DM_rs(i,1) = (min(UnitResponse(po_i(1):po_i(1)+trialSide_samples))-min(UnitResponse(po_i(1)-trialSide_samples:po_i(1))))/...
%                 (min(UnitResponse(po_i(1):po_i(1)+trialSide_samples))+min(UnitResponse(po_i(1)-trialSide_samples:po_i(1))));
%             index_i(i) = i;
%     end 
     
    % still trials
    Responses = Responses_ON & Param_still(1:length(Responses_ON),i);
    UnitResponse = mean(squeeze(UnitResponses_smooth(Responses==1,i,:)))*60;
    %UnitResponse = UnitResponse - mean(UnitResponse(po_i(1)-trialSide_samples:po_i(1)));
    Dir = sign(sum(UnitResponse(po_i(1):po_i(1)+trialSide_samples))-sum(UnitResponse(po_i(1)-trialSide_samples:po_i(1))));
    
    DM_rs(i,2) = (nanmax(UnitResponse(po_i(1):po_i(1)+trialSide_samples))-nanmax(UnitResponse(po_i(1)-trialSide_samples:po_i(1))))/...
                    nanmax(UnitResponse(po_i(1)-trialSide_samples:po_i(1)));
%     DM_rs(i,2) = (nanmean(UnitResponse(po_i(1):po_i(1)+trialSide_samples))-nanmean(UnitResponse(po_i(1)-trialSide_samples:po_i(1))))/...
%                     nanmean(UnitResponse(po_i(1)-trialSide_samples:po_i(1)));
    
    Peak_rs(i,2) = nanmax(UnitResponse(po_i(1):po_i(1)+trialSide_samples));
    PopRes_still(i,:) = UnitResponse;
%     switch Dir
%         case 1
%             DM_rs(i,2) = (max(UnitResponse(po_i(1):po_i(1)+trialSide_samples))-max(UnitResponse(po_i(1)-trialSide_samples:po_i(1))))/...
%                 (max(UnitResponse(po_i(1):po_i(1)+trialSide_samples))+max(UnitResponse(po_i(1)-trialSide_samples:po_i(1))));
%             index_i(i) = i;
%         case -1
%             DM_rs(i,2) = (min(UnitResponse(po_i(1):po_i(1)+trialSide_samples))-min(UnitResponse(po_i(1)-trialSide_samples:po_i(1))))/...
%                 (min(UnitResponse(po_i(1):po_i(1)+trialSide_samples))+min(UnitResponse(po_i(1)-trialSide_samples:po_i(1))));
%             index_i(i) = i;
%     end
     
end
% DM_rs(size(DM_rs,1)+1:size(Param_run,2),:) = zeros(size(Param_run,2)-size(DM_rs,1),2);
% Peak_rs(size(Peak_rs,1)+1:size(Param_run,2),:) = zeros(size(Param_run,2)-size(DM_rs,1),2);
% PopRes_run(size(DM_rs,1)+1:size(Param_run,2),:) = zeros(size(Param_run,2)-size(DM_rs,1),2);
% PopRes_still(size(DM_rs,1)+1:size(Param_run,2),:) = zeros(size(Param_run,2)-size(DM_rs,1),2);
%% modulation index of grating stimulus responses
Condition = 1;
switch Condition
    case 1
        Param_run = SelectedTrials_OI{Condition,1}(:,SelectedCells);
        Param_still = SelectedTrials_OI{Condition,2}(:,SelectedCells);
    case 2
        Param_run = SelectedTrials_OI{Condition,1}(:,SelectedCells);
        Param_still = SelectedTrials_OI{Condition,2}(:,SelectedCells);
    case 3
        Param_run = SelectedTrials_OI{Condition,1}(:,SelectedCells);
        Param_still = SelectedTrials_OI{Condition,2}(:,SelectedCells);
    case 4
        Param_run = SelectedTrials_OI{Condition,1}(:,SelectedCells);
        Param_still = SelectedTrials_OI{Condition,2}(:,SelectedCells);
end
% Modulation Index - MI - grating stimulus
r = 1; clear DM_rs_V DepthMod Responses ResponseON 
for i = 1:size(Param,2)
    % select trial responses
    Responses_ON = Param(:,i);
    % scrub responses
    Responses_ON = logical(Responses_ON(~isnan(Responses_ON)));
    
    % running trials 
    Responses = Responses_ON & Param_run(1:length(Responses_ON),i);
    UnitResponse = nanmean(squeeze(UnitResponses_smooth(Responses==1,i,:)))*60; 
    Dir = sign(sum(UnitResponse(trialSide_samples+1:trialSide_samples+240))-sum(UnitResponse(1:trialSide_samples)));
    DM_rs_V(i,1) = (nanmax(UnitResponse(trialSide_samples+1:trialSide_samples+240))-nanmax(UnitResponse(1:trialSide_samples)))/...
                       nanmax(UnitResponse(1:trialSide_samples));
%     switch Dir
%         case 1
%             DM_rs_V(i,1) = (max(UnitResponse(trialSide_samples+1:trialSide_samples+240))-max(UnitResponse(1:trialSide_samples)))/...
%                 (max(UnitResponse(trialSide_samples+1:trialSide_samples+240))+max(UnitResponse(1:trialSide_samples)));
%             index_i(i) = i;
%         case -1
%             DM_rs_V(i,1) = (min(UnitResponse(trialSide_samples+1:trialSide_samples+240))-min(UnitResponse(1:trialSide_samples)))/...
%                 (min(UnitResponse(trialSide_samples+1:trialSide_samples+240))+min(UnitResponse(1:trialSide_samples)));
%             index_i(i) = i;
%     end 
     
    % still trials
    Responses = Responses_ON & Param_still(1:length(Responses_ON),i);
    UnitResponse = mean(squeeze(UnitResponses_smooth(Responses==1,i,:)))*60;
    Dir = sign(sum(UnitResponse(trialSide_samples+1:trialSide_samples+240))-sum(UnitResponse(1:trialSide_samples)));
    DM_rs_V(i,2) = (nanmax(UnitResponse(trialSide_samples+1:trialSide_samples+240))-nanmax(UnitResponse(1:trialSide_samples)))/...
                    nanmax(UnitResponse(1:trialSide_samples));
%     switch Dir
%         case 1
%             DM_rs_V(i,2) = (max(UnitResponse(trialSide_samples+1:trialSide_samples+240))-max(UnitResponse(1:trialSide_samples)))/...
%                 (max(UnitResponse(trialSide_samples+1:trialSide_samples+240))+max(UnitResponse(1:trialSide_samples)));
%             index_i(i) = i;
%         case -1
%             DM_rs_V(i,2) = (min(UnitResponse(trialSide_samples+1:trialSide_samples+240))-min(UnitResponse(1:trialSide_samples)))/...
%                 (min(UnitResponse(trialSide_samples+1:trialSide_samples+240))+min(UnitResponse(1:trialSide_samples)));
%             index_i(i) = i;
%     end
     
end
DM_rs_V(size(DM_rs_V,1)+1:size(Param_run,2),:) = zeros(size(Param_run,2)-size(DM_rs_V,1),2);


%% sort the cells by their response to perturbation (reference for heat map)
Units_Sel = AM_UOI;
[B,I,PertTrial_Ex] = PlotSortingOrder(ProjectData,SelectedResponses,AM_Param,AM_Speed,Units_Sel,0);
% I = descending order taking into account:[pert response; 2s vis stim respo; post-pert reponse]
I = I(~isnan(B(:,1)));
B = reshape(B(~isnan(B(:))),length(B(~isnan(B(:))))/3,3);

clear DM_rs_i DM_rs_Vi
DM_rs_i = DM_rs(PertRespUnits_pos,:);
Excl_i = ~isnan(DM_rs_i(:,1)) & ~isnan(DM_rs_i(:,2)) & (DM_rs_i(:,1)~=0 & DM_rs_i(:,2)~=0);
DM_rs_i = DM_rs_i(Excl_i,:);

DM_rs_Vi = DM_rs_V(PertRespUnits_pos,:);
Excl_iV = ~isnan(DM_rs_Vi(:,1)) & ~isnan(DM_rs_Vi(:,2)) & (DM_rs_Vi(:,1)~=0 & DM_rs_Vi(:,2)~=0);


figure
plot(DM_rs_i(:,2),DM_rs_i(:,1),'.r','MarkerSize',8)
xlabel('MI still'); ylabel('MI run');
set(gca,'TickDir','out','FontSize',12,'FontName','Arial')
xlim([-1 1]); ylim([-1 1]);
grid on
box off
hold on
plot([-1 1], [-1 1], ':k')
title(['MI_r_u_n VS MI_s_t_i_l_l n = ' num2str(length(DM_rs_i))])
[p h] = signrank(DM_rs_i(:,1)-DM_rs_i(:,2))
text(-0.9,-0.3,['p = ' num2str(p,2)])
figure
histogram(DM_rs_i(:,1)-DM_rs_i(:,2),[-1.5:0.1:1.5],'FaceColor','r','Normalization','probability','EdgeColor','none')
xlabel(['MI_r_u_n - MI_s_t_i_l_l n = ' num2str(length(DM_rs_i))]);
set(gca,'TickDir','out','FontSize',12,'FontName','Arial')
ylabel('Probability');
box off
hold on
title(['MI_r_u_n - MI_s_t_i_l_l n = ' num2str(length(DM_rs_i))])

% plot the peaks responses
Peak_rs_i = Peak_rs(PertRespUnits_pos,:);
Peak_rs_i = Peak_rs_i(Excl_i,:);
Peak_rs_i = Peak_rs_i(~isnan(Peak_rs_i(:,1)) | ~isnan(Peak_rs_i(:,2)),:);
% Peak_rs_i_n = Peak_rs_i./MaxValues; % from 15 lines below
figure
plot(Peak_rs_i(:,2),Peak_rs_i(:,1),'.r','MarkerSize',8)
%plot(Peak_rs_i_n(:,2),Peak_rs_i_n(:,1),'.r','MarkerSize',8)
set(gca,'TickDir','out','FontSize',12,'FontName','Arial')
grid on
box off
hold on
plot([0 max(Peak_rs_i(:,2))], [0 max(Peak_rs_i(:,2))], ':k')
%plot([min(Peak_rs_i_n(:,2)) max(Peak_rs_i_n(:,2))], [min(Peak_rs_i_n(:,1)) max(Peak_rs_i_n(:,1))], ':k')
xlabel('FR Peak still'); ylabel('FR Peak run');
title(['Peak_r_u_n VS Peak_s_t_i_l_l n = ' num2str(length(DM_rs_i))])
text(10,30,['p=' num2str(p)])
[p h] = signrank(Peak_rs_i(:,1)-Peak_rs_i(:,2))
[p h] = signrank(Peak_rs_i_n(:,1)-Peak_rs_i_n(:,2))


PopRes_run_i = PopRes_run(PertRespUnits_pos,:);

PopRes_still_i = PopRes_still(PertRespUnits_pos,:); 
MaxValues = nanmax(nanmax(PopRes_run_i,[],2),nanmax(PopRes_still_i,[],2));
figure
plot(nanmean(PopRes_run_i./MaxValues),'r')
hold on
plot(nanmean(PopRes_still_i./MaxValues),'k')

%% re-compute MI with data from SelectedResponses
figure
plot(nanmean(SelectedResponses{1,1}(PertRespUnits_pos,:)),'r')
hold on
plot(nanmean(SelectedResponses{1,2}(PertRespUnits_pos,:)),'k')

figure
for i = 1:length(DM_rs_i)
   plot([1 2], [DM_rs_i(:,1) DM_rs_i(:,2)],'Color',[0.5 0.5 0.5])
   hold on
end
boxplot(DM_rs_i,'Notch','on','labels',{'MI_run', 'MI_still'})
ylabel('Modulation Index');
set(gca,'TickDir','out','box','off','fontSize',13)
title('Effect of running')
set(0, 'DefaultFigureRenderer', 'painters');

%% plot modulation index of grating stimuli runnning 
DM_rs_Vi = DM_rs_Vi(Excl_iV,:);

figure
plot(DM_rs_Vi(:,2),DM_rs_Vi(:,1),'.r','MarkerSize',8)
xlabel('MI still'); ylabel('MI run');
set(gca,'TickDir','out','FontSize',12,'FontName','Arial')
xlim([-1 1]); ylim([-1 1]);
grid on
box off
hold on
plot([-1 1], [-1 1], ':k')
title(['MI_r_u_n VS MI_s_t_i_l_l n = ' num2str(length(DM_rs_Vi))])

figure
histogram(DM_rs_Vi(:,1)-DM_rs_Vi(:,2),[-1.5:0.1:1.5],'FaceColor','r','Normalization','probability','EdgeColor','none')
xlabel(['MI_r_u_n - MI_s_t_i_l_l n = ' num2str(length(DM_rs_Vi))]);
set(gca,'TickDir','out','FontSize',12,'FontName','Arial')
ylabel('Probability');
box off
hold on
title(['MI_r_u_n - MI_s_t_i_l_l n = ' num2str(length(DM_rs_Vi))])
[h p] = signrank(DM_rs_Vi(:,1)-DM_rs_Vi(:,2))


figure
for i = 1:length(DM_rs_Vi)
   plot([1 2], [DM_rs_Vi(:,1) DM_rs_Vi(:,2)],'Color',[0.5 0.5 0.5])
   hold on
end
boxplot(DM_rs_Vi,'Notch','on','labels',{'MI_run', 'MI_still'})
ylabel('Modulation Index');
set(gca,'TickDir','out','box','off','fontSize',13)
title('Effect of running')
set(0, 'DefaultFigureRenderer', 'painters');


%% look at running modulation of single units
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
Param = squeeze(AM_Param(:,SelectedCells,3)); % pert ON/OFF for all units

r = 1; clear DM_rs DepthMod Responses ResponseON Peak_rs PopRes_run PopRes_still ranksum_p ranksum_h  AreaUC
clear ranksum_p_mean ranksum_h_mean ranksum_p_mean UnitResponseNP_g UnitResponsePT_g ttest_h_mean ttest_p_mean DM_RS_
for i = 1:size(Param,2)
    % select perturbation trial responses
    Responses_ON = Param(:,i);
    % scrub responses
    Responses_ON = logical(Responses_ON(~isnan(Responses_ON)));
    %% save all trials for the two conditions
    clear UnitResponse_r UnitResponse_s
    % running trials 
    Responses = Responses_ON & Param_run(1:length(Responses_ON),i);
    UnitResponse_run = squeeze(UnitResponses_smooth(Responses==1,i,:))*60; 
    % still trials
    Responses = Responses_ON & Param_still(1:length(Responses_ON),i);
    UnitResponse_still = squeeze(UnitResponses_smooth(Responses==1,i,:))*60;
    
    % All perturbation trials
    Responses = Responses_ON;
    UnitResponsePT_g(i,:) = nanmean(squeeze(UnitResponses_smooth(Responses==1,i,:))*60);
    % normal trials
    Responses_OFF = Param(:,i)==0;
    % scrub responses
    Responses_OFF = Responses_OFF(1:length(Responses_ON));
    Responses = Responses_OFF;
    UnitResponseNP_g(i,:) = nanmean(squeeze(UnitResponses_smooth(Responses==1,i,:))*60);
    % find the smalle sample size
    [SampleSize SampleId] = min([size(UnitResponse_run,1) size(UnitResponse_still,1)]);
    %% select equal number of trials for each condition
    for k = 1:1000
        %% re-sample trials with replacement
        if size(UnitResponse_run,1)==0 | size(UnitResponse_still,1)==0
            trials_run = 0; trials_still = 0;
            UnitResponse_r_mean_t = nanmean(UnitResponse_run(1:SampleSize,po_i(1)-trialSide_samples:po_i(1)+trialSide_samples));
            UnitResponse_s_mean_t = nanmean(UnitResponse_still(1:SampleSize,po_i(1)-trialSide_samples:po_i(1)+trialSide_samples));
        else
            trials_run = randi([1 size(UnitResponse_run,1)],1,20);
            trials_still = randi([1 size(UnitResponse_still,1)],1,20);
            UnitResponse_r_mean_t = nanmean(UnitResponse_run(trials_run,po_i(1)-trialSide_samples:po_i(1)+trialSide_samples));
            UnitResponse_s_mean_t = nanmean(UnitResponse_still(trials_still,po_i(1)-trialSide_samples:po_i(1)+trialSide_samples));
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
        AreaUC(i,k,1) = sum(UnitResponse_r_mean(k,61:end)-mean(UnitResponse_r_mean_t(1:60))); 
        AreaUC(i,k,2) = sum(UnitResponse_s_mean(k,61:end)-mean(UnitResponse_s_mean_t(1:60))); 
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

clear DM_rs_i DM_rs_Vi
DM_rs_i = DM_rs(PertRespUnits_pos,:);
Excl_i = ~isnan(DM_rs_i(:,1)) & ~isnan(DM_rs_i(:,2)) & (DM_rs_i(:,1)~=0 & DM_rs_i(:,2)~=0);
DM_rs_i = DM_rs_i(Excl_i,:);

DM_rs_Vi = DM_rs_V(PertRespUnits_pos,:);
Excl_iV = ~isnan(DM_rs_Vi(:,1)) & ~isnan(DM_rs_Vi(:,2)) & (DM_rs_Vi(:,1)~=0 & DM_rs_Vi(:,2)~=0);

PopRes_run_i = PopRes_run(PertRespUnits_pos,:);
PopRes_still_i = PopRes_still(PertRespUnits_pos,:);
PopRes_run_i = PopRes_run_i(Excl_i,:);
PopRes_still_i = PopRes_still_i(Excl_i,:);

UnitResponsePT_g_i = UnitResponsePT_g(PertRespUnits_pos,:);
UnitResponsePT_g_i = UnitResponsePT_g_i(Excl_i,:);
UnitResponseNP_g_i = UnitResponseNP_g(PertRespUnits_pos,:);
UnitResponseNP_g_i = UnitResponseNP_g_i(Excl_i,:);

ranksum_p_mean_i = ranksum_p_mean(PertRespUnits_pos);
ranksum_p_mean_i = ranksum_p_mean_i(Excl_i);
ttest_p_mean_i = ttest_p_mean(PertRespUnits_pos);
ttest_p_mean_i = ttest_p_mean_i(Excl_i);

ranksum_p_i = ranksum_p(PertRespUnits_pos,:);
ranksum_p_i = ranksum_p_i(Excl_i,:);
for j = 1:size(ranksum_p_i,1)
    ind_p_01(j) = sum(ranksum_p_i(j,:)<0.01);
end

DM_RS_i = DM_RS_(PertRespUnits_pos,:,:);
DM_RS_i = DM_RS_i(Excl_i,:,:);

AreaUC_i = AreaUC(PertRespUnits_pos,:,:);
AreaUC_i = AreaUC_i(Excl_i,:,:);
jj = 1; clear index_of_i
figure
for j = 1:size(PopRes_run_i,1)
    Peak(j,1) = max( PopRes_run_i(j,po_i(1):po_i(1)+trialSide_samples) - mean(PopRes_run_i(j,po_i(1)-trialSide_samples:po_i(1))) );
    Peak(j,2) = max( PopRes_still_i(j,po_i(1):po_i(1)+trialSide_samples) - mean(PopRes_still_i(j,po_i(1)-trialSide_samples:po_i(1))) );
    
    if j==91
        figure
    end
    if j>90
        subplot(9,10,j-90)
    else
        subplot(9,10,j)
    end
    plot(linspace(-1,1,length(PopRes_run_i(j,po_i(1)-trialSide_samples:po_i(1)+trialSide_samples))),...
        PopRes_run_i(j,po_i(1)-trialSide_samples:po_i(1)+trialSide_samples) - mean(PopRes_run_i(j,po_i(1)-trialSide_samples:po_i(1))) ...
        ,'r')
    hold on
    plot(linspace(-1,1,length(PopRes_run_i(j,po_i(1)-trialSide_samples:po_i(1)+trialSide_samples))),...
        PopRes_still_i(j,po_i(1)-trialSide_samples:po_i(1)+trialSide_samples) - mean(PopRes_still_i(j,po_i(1)-trialSide_samples:po_i(1))) ...
        ,'k')
    plot(linspace(-1,1,length(PopRes_run_i(j,po_i(1)-trialSide_samples:po_i(1)+trialSide_samples))),...
        UnitResponsePT_g_i(j,po_i(1)-trialSide_samples:po_i(1)+trialSide_samples) - mean(UnitResponsePT_g_i(j,po_i(1)-trialSide_samples:po_i(1))) ...
        ,'-','Color',[255 192 203]/255)
    plot(linspace(-1,1,length(PopRes_run_i(j,po_i(1)-trialSide_samples:po_i(1)+trialSide_samples))),...
        UnitResponseNP_g_i(j,po_i(1)-trialSide_samples:po_i(1)+trialSide_samples) - mean(UnitResponseNP_g_i(j,po_i(1)-trialSide_samples:po_i(1))) ...
        ,'-','Color',[0.5 0.5 0.5])
    % if ((1000-ind_p_01(j))/1000)<=0.05 & Peak(j,1)>Peak(j,2)
    %if (ttest_p_mean_i(j)<=0.05 & Peak(j,1)>Peak(j,2))
    if ( (sum(AreaUC_i(j,:,1)>AreaUC_i(j,:,2))>950) & Peak(j,1)>Peak(j,2) )
        plot(linspace(-1,1,length(PopRes_run_i(j,po_i(1)-trialSide_samples:po_i(1)+trialSide_samples))),...
            PopRes_run_i(j,po_i(1)-trialSide_samples:po_i(1)+trialSide_samples) - mean(PopRes_run_i(j,po_i(1)-trialSide_samples:po_i(1))) ...
            ,'r','LineWidth',3)
        index_of_i(jj) = j; 
        jj = jj + 1;
    end
    % title(['MI_r=' num2str(DM_rs_i(j,1),2) ' MI_s' num2str(DM_rs_i(j,2),2) ' p=' num2str((1000-ind_p_01(j))/1000,2)])
    %title(['MI_r=' num2str(DM_rs_i(j,1),2) ' MI_s' num2str(DM_rs_i(j,2),2) ' p=' num2str(ttest_p_mean_i(j),2)])
    title(['MI_r=' num2str(DM_rs_i(j,1),2) ' MI_s=' num2str(DM_rs_i(j,2),2) ' p=' num2str((1000-sum(AreaUC_i(j,:,1)>AreaUC_i(j,:,2)))/1000,2)])
    %title(['Ratio=' num2str(mean(AreaUC_i(j,:,1))/mean(AreaUC_i(j,:,2)),2)  ' p=' num2str((1000-sum(AreaUC_i(j,:,1)>AreaUC_i(j,:,2)))/1000,2)])
    box off
    grid on
end

%% highlights units that show MI change due to modulation
figure
plot(DM_rs_i(:,2),DM_rs_i(:,1),'.r','MarkerSize',8)
hold on
plot(DM_rs_i(index_of_i,2),DM_rs_i(index_of_i,1),'+r','MarkerSize',12)
xlabel('MI still'); ylabel('MI run');
set(gca,'TickDir','out','FontSize',12,'FontName','Arial')
xlim([-1 1]); ylim([-1 1]);
grid on
box off
hold on
plot([-1 1], [-1 1], ':k')
title(['MI_r_u_n VS MI_s_t_i_l_l n = ' num2str(length(DM_rs_i))])
[p h] = signrank(DM_rs_i(:,1)-DM_rs_i(:,2))
[p h] = signrank(DM_rs_i(index_of_i,1)-DM_rs_i(index_of_i,2))


figure
plot(Peak_rs_i(:,2),Peak_rs_i(:,1),'.r','MarkerSize',8)
%plot(Peak_rs_i_n(:,2),Peak_rs_i_n(:,1),'.r','MarkerSize',8)
set(gca,'TickDir','out','FontSize',12,'FontName','Arial')
grid on
box off
hold on
plot([0 max(Peak_rs_i(:,2))], [0 max(Peak_rs_i(:,2))], ':k')
%plot([min(Peak_rs_i_n(:,2)) max(Peak_rs_i_n(:,2))], [min(Peak_rs_i_n(:,1)) max(Peak_rs_i_n(:,1))], ':k')
xlabel('FR Peak still'); ylabel('FR Peak run');
title(['Peak_r_u_n VS Peak_s_t_i_l_l n = ' num2str(length(DM_rs_i))])
text(10,30,['p=' num2str(p)])
[p h] = signrank(Peak_rs_i(:,1)-Peak_rs_i(:,2))
[p h] = signrank(Peak_rs_i_n(:,1)-Peak_rs_i_n(:,2))

figure
plot(mean(AreaUC_i(:,:,2),2),mean(AreaUC_i(:,:,1),2),'or','MarkerSize',8)
hold on
plot(mean(AreaUC_i(index_of_i,:,2),2),mean(AreaUC_i(index_of_i,:,1),2),'.r','MarkerSize',16)
set(gca,'TickDir','out','FontSize',12,'FontName','Arial')
grid on
box off
hold on
plot([-500 500], [-500 500], ':k')
xlabel('AUC still'); ylabel('AUC run');
title(['AUC_r_u_n VS AUC_s_t_i_l_l n = ' num2str(length(DM_rs_i))])
[p h] = signrank(mean(AreaUC_i(:,:,2),2)-mean(AreaUC_i(:,:,1),2))
text(500,30,['p=' num2str(p)])
[p h] = signrank(AreaUC_i(index_of_i,1,2)-AreaUC_i(index_of_i,1,1))
text(500,-300,['p=' num2str(p)])
xlim([-500 1000])

figure
boxplot(AreaUC_i(:,1,2)./AreaUC_i(:,1,1),'Notch','on')
hold on
boxplot(AreaUC_i(index_of_i,1,2)./AreaUC_i(index_of_i,1,1),'Notch','on',)


%%
clear Peak
for j = 1:size(PopRes_run_i,1)
  
    Peak(j,1) = max( PopRes_run_i(j,po_i(1)-trialSide_samples:po_i(1)+trialSide_samples) - mean(PopRes_run_i(j,po_i(1)-trialSide_samples:po_i(1))) );
    Peak(j,2) = max( PopRes_still_i(j,po_i(1)-trialSide_samples:po_i(1)+trialSide_samples) - mean(PopRes_still_i(j,po_i(1)-trialSide_samples:po_i(1))) );
end

find( ((1000-ind_p_01)/1000<=0.05)' & Peak(:,1)>Peak(:,2))



% count how many units are recorded from these sessions

PopRes_run_i = PopRes_run(PertRespUnits_pos,:);

Rec_OI = squeeze(AM_Param(1,SelectedCells,1)); % recording ID
Rec_OI_i = Rec_OI(PertRespUnits_pos);
Rec_OI_i = Rec_OI_i(Excl_i);
Rec_OI = unique(Rec_OI_i);
AM_Param_ = AM_Param(1,AM_UOI,1);
totUnits = 0;
for i = 1:length(Rec_OI)
    totUnits = sum(AM_Param_ == Rec_OI(i)) + totUnits
end


%% example
figure
for j = 34
    Peak(j,1) = max( PopRes_run_i(j,po_i(1):po_i(1)+trialSide_samples) - mean(PopRes_run_i(j,po_i(1)-trialSide_samples:po_i(1))) );
    Peak(j,2) = max( PopRes_still_i(j,po_i(1):po_i(1)+trialSide_samples) - mean(PopRes_still_i(j,po_i(1)-trialSide_samples:po_i(1))) );
  
    plot(linspace(-1,1,length(PopRes_run_i(j,po_i(1)-trialSide_samples:po_i(1)+trialSide_samples))),...
        PopRes_run_i(j,po_i(1)-trialSide_samples:po_i(1)+trialSide_samples)  ...
        ,'r')
    hold on
    plot(linspace(-1,1,length(PopRes_run_i(j,po_i(1)-trialSide_samples:po_i(1)+trialSide_samples))),...
        PopRes_still_i(j,po_i(1)-trialSide_samples:po_i(1)+trialSide_samples)  ...
        ,'k')
    plot(linspace(-1,1,length(PopRes_run_i(j,po_i(1)-trialSide_samples:po_i(1)+trialSide_samples))),...
        UnitResponsePT_g_i(j,po_i(1)-trialSide_samples:po_i(1)+trialSide_samples)  ...
        ,'-','Color',[255 192 203]/255)
    plot(linspace(-1,1,length(PopRes_run_i(j,po_i(1)-trialSide_samples:po_i(1)+trialSide_samples))),...
        UnitResponseNP_g_i(j,po_i(1)-trialSide_samples:po_i(1)+trialSide_samples) ...
        ,'-','Color',[0.5 0.5 0.5])
    % if ((1000-ind_p_01(j))/1000)<=0.05 & Peak(j,1)>Peak(j,2)
    %if (ttest_p_mean_i(j)<=0.05 & Peak(j,1)>Peak(j,2))
    if ( (sum(AreaUC_i(j,:,1)>AreaUC_i(j,:,2))>950) & Peak(j,1)>Peak(j,2) )
        plot(linspace(-1,1,length(PopRes_run_i(j,po_i(1)-trialSide_samples:po_i(1)+trialSide_samples))),...
            PopRes_run_i(j,po_i(1)-trialSide_samples:po_i(1)+trialSide_samples)  ...
            ,'r','LineWidth',3)
    end
    % title(['MI_r=' num2str(DM_rs_i(j,1),2) ' MI_s' num2str(DM_rs_i(j,2),2) ' p=' num2str((1000-ind_p_01(j))/1000,2)])
    %title(['MI_r=' num2str(DM_rs_i(j,1),2) ' MI_s' num2str(DM_rs_i(j,2),2) ' p=' num2str(ttest_p_mean_i(j),2)])
    title(['MI_r=' num2str(DM_rs_i(j,1),2) ' MI_s=' num2str(DM_rs_i(j,2),2) ' p=' num2str((1000-sum(AreaUC_i(j,:,1)>AreaUC_i(j,:,2)))/1000,2)])
    %title(['Ratio=' num2str(mean(AreaUC_i(j,:,1))/mean(AreaUC_i(j,:,2)),2)  ' p=' num2str((1000-sum(AreaUC_i(j,:,1)>AreaUC_i(j,:,2)))/1000,2)])
    box off
    grid on
end

