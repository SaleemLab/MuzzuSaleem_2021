%% Tomaso Muzzu - UCL - 23/03/2020 - running modulation effects

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
    load('DM_pert_shuffled_CTRL.mat')
    DM_CTRL = DM;
    load('DM_pert_shuffled.mat')
    DM = cat(1,DM_sh(:,1),DM');
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
    switch Dir
        case 1
            DM_rs(i,1) = (max(UnitResponse(po_i(1):po_i(1)+trialSide_samples))-max(UnitResponse(po_i(1)-trialSide_samples:po_i(1))))/...
                (max(UnitResponse(po_i(1):po_i(1)+trialSide_samples))+max(UnitResponse(po_i(1)-trialSide_samples:po_i(1))));
            index_i(i) = i;
        case -1
            DM_rs(i,1) = (min(UnitResponse(po_i(1):po_i(1)+trialSide_samples))-min(UnitResponse(po_i(1)-trialSide_samples:po_i(1))))/...
                (min(UnitResponse(po_i(1):po_i(1)+trialSide_samples))+min(UnitResponse(po_i(1)-trialSide_samples:po_i(1))));
            index_i(i) = i;
    end 
     
    % still trials
    Responses = Responses_ON & Param_still(1:length(Responses_ON),i);
    UnitResponse = mean(squeeze(UnitResponses_smooth(Responses==1,i,:)))*60;
    Dir = sign(sum(UnitResponse(po_i(1):po_i(1)+trialSide_samples))-sum(UnitResponse(po_i(1)-trialSide_samples:po_i(1))));
    switch Dir
        case 1
            DM_rs(i,2) = (max(UnitResponse(po_i(1):po_i(1)+trialSide_samples))-max(UnitResponse(po_i(1)-trialSide_samples:po_i(1))))/...
                (max(UnitResponse(po_i(1):po_i(1)+trialSide_samples))+max(UnitResponse(po_i(1)-trialSide_samples:po_i(1))));
            index_i(i) = i;
        case -1
            DM_rs(i,2) = (min(UnitResponse(po_i(1):po_i(1)+trialSide_samples))-min(UnitResponse(po_i(1)-trialSide_samples:po_i(1))))/...
                (min(UnitResponse(po_i(1):po_i(1)+trialSide_samples))+min(UnitResponse(po_i(1)-trialSide_samples:po_i(1))));
            index_i(i) = i;
    end
     
end
DM_rs(size(DM_rs,1)+1:size(Param_run,2),:) = zeros(size(Param_run,2)-size(DM_rs,1),2);


%% sort the cells by their response to perturbation (reference for heat map)
Units_Sel = AM_UOI;
[B,I,PertTrial_Ex] = PlotSortingOrder(ProjectData,SelectedResponses,AM_Param,AM_Speed,Units_Sel,0);
% I = descending order taking into account:[pert response; 2s vis stim respo; post-pert reponse]
I = I(~isnan(B(:,1)));
B = reshape(B(~isnan(B(:))),length(B(~isnan(B(:))))/3,3);

clear DM_rs_i
DM_rs_i = DM_rs(PertRespUnits_pos,:);
DM_rs_i(isnan(DM_rs_i(:,1)),:)=[];
DM_rs_i(isnan(DM_rs_i(:,2)),:)=[];
DM_rs_i(DM_rs_i(:,1)==0 & DM_rs_i(:,2)==0,:)=[];

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

figure
histogram(DM_rs_i(:,1)-DM_rs_i(:,2),[-1.5:0.1:1.5],'FaceColor','r','Normalization','probability','EdgeColor','none')
xlabel(['MI_r_u_n - MI_s_t_i_l_l n = ' num2str(length(DM_rs_i))]);
set(gca,'TickDir','out','FontSize',12,'FontName','Arial')
ylabel('Probability');
box off
hold on
title(['MI_r_u_n - MI_s_t_i_l_l n = ' num2str(length(DM_rs_i))])
[h p] = kstest(DM_rs_i(:,1)-DM_rs_i(:,2))


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


