%% Tomaso Muzzu - UCL - 17/09/2019

% quantify neural responses to perturbation 

% predictors: 
    % 1 - ratio of FR during pert period over FR of 1s before
    % 2 - difference of FR during pert period mins FR of 1s before
    % 3 - integral of the mean FR during pert period
    % 4 - integral of the mean FR from -1s to end of pert period
% same predictors for post perturbatino period
% responses: the actual visual perturbation


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
minTrials = 5; % min nr of trials for each condition
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
    
%     % running  trials
%     Responses = Responses_ON & Param_run(1:length(Responses_ON),i);
%     UnitResponse = mean(squeeze(UnitResponses_smooth(Responses==1,i,:)))*60;
%     DM_rs(i,1) = (sum(UnitResponse(po_i(1):po_i(1)+trialSide_samples))-sum(UnitResponse(po_i(1)-trialSide_samples:po_i(1))))/...
%          (sum(UnitResponse(po_i(1):po_i(1)+trialSide_samples))+sum(UnitResponse(po_i(1)-trialSide_samples:po_i(1))));
%     % still trials
%     Responses = Responses_ON & Param_still(1:length(Responses_ON),i);
%     UnitResponse = mean(squeeze(UnitResponses_smooth(Responses==1,i,:)))*60;
%     DM_rs(i,2) = (sum(UnitResponse(po_i(1):po_i(1)+trialSide_samples))-sum(UnitResponse(po_i(1)-trialSide_samples:po_i(1))))/...
%          (sum(UnitResponse(po_i(1):po_i(1)+trialSide_samples))+sum(UnitResponse(po_i(1)-trialSide_samples:po_i(1)))); 
    
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
    
%     if AUC_units(i,1)>p_sh(2)
%         DepthMod(r,1:2) = [DM(i) AUC_units(i)];
%         r = r + 1;
%     end
%     parfor sh=1:1000
%         % shuffle perturbation trals
%         Responses_sh = Responses(randperm(length(Responses)));
%         UnitResponse = mean(squeeze(UnitResponses_smooth(Responses_sh==1,i,:)))*60;
%         DM_sh(i,sh) = (sum(UnitResponse(po_i(1):po_i(1)+trialSide_samples))-sum(UnitResponse(po_i(1)-trialSide_samples:po_i(1))))/...
%             (sum(UnitResponse(po_i(1):po_i(1)+trialSide_samples))+sum(UnitResponse(po_i(1)-trialSide_samples:po_i(1))));
%     end
end

DM_rs(isnan(DM_rs(:,1)),:)=[];
DM_rs(isnan(DM_rs(:,2)),:)=[];
DM_rs(DM_rs(:,1)==0 & DM_rs(:,1)==0,:)=[];

figure
plot(DM_rs(:,2),DM_rs(:,1),'.','MarkerSize',8)
xlabel('MI still'); ylabel('MI run');
set(gca,'TickDir','out','FontSize',12,'FontName','Arial')
xlim([-1 1]); ylim([-1 1]);
grid on
box off
hold on
plot([-1 1], [-1 1], ':k')
title(['DMI_r_u_n VS DMI_s_t_i_l_l n = ' num2str(sum(~isnan(DM_rs(:,1))))])


figure
histogram(DM_rs(:,1),-1:0.05:1)
hold on
histogram(DM_rs(:,2),-1:0.05:1)

% plot box plot showing change of magnitude of response in run/still 
figure
for i = 1:length(DM_rs)
    plot([DM_rs(i,1) DM_rs(i,2)],'Color',[0.5 0.5 0.5])
    hold on
end
xlim([0 3]); ylim([-1 1]);
set(gca,'XTick',[1 2],'XTicklabel',{'run','still'},'TickDir','out')
ylabel('Modulation index'); box off;
plot([nanmean(DM_rs(:,1)) nanmean(DM_rs(:,2))],'k','LineWidth',3)
title('DMI_r_u_n VS DMI_s_t_i_l_l (min 10 trials/condition)');

% difference between run and still conditions for each unit
for i = 1:length(DM_rs(:,1))
    Dist(i) = pdist([DM_rs(i,1) 0; DM_rs(i,2) 0]);
end

MeanDelta = median(Dist);
DM_rs(:,1)-DM_rs(:,2);


%% sort the cells by their response to perturbation (reference for heat map)
[B,I,PertTrial_Ex] = PlotSortingOrder(ProjectData,SelectedResponses,AM_Param,AM_Speed,Units_Sel,0);
% I = descending order taking into account:[pert response; 2s vis stim respo; post-pert reponse]
I = I(~isnan(B(:,1)));
B = reshape(B(~isnan(B(:))),length(B(~isnan(B(:))))/3,3);

%% show mean curves for running and stationary of all units
Units_Index = find(DM_rs(:,1)~=0 & DM_rs(:,2)~=0);
Units_Index_still = find(nansum(SelTrials{1,2},1));
SelTrials{1} = SelectedTrials_OI{2,1}(:,SelectedCells);
SelTrials{2} = SelectedTrials_OI{2,2}(:,SelectedCells);
Units_Index = find(nansum(SelTrials{1,1},1));
param_sel = AM_Param(:,SelectedCells,:);
fig_n=1;
for k = 1:size(Units_Index,2)
   if mod(k,16)==0
       subplot_nr = 16;
   elseif mod(k,16)==1
       FIG = figure('units','normalized','outerposition',[0 0 1 1]);
       subplot_nr = mod(k,16);
   else
       subplot_nr = mod(k,16);
   end
   TimeLine = linspace(-1,8.3,length(SelectedResponses{1,1}(Units_Index(k),:)));
   Pert_lims =[param_sel(min(find(param_sel(:,Units_Index(k),3)==1)),Units_Index(k),4) param_sel(min(find(param_sel(:,Units_Index(k),3)==1)),Units_Index(k),5)]; 
   subplot(4,4,subplot_nr)
   plot(TimeLine,SelectedResponses{1,1}(Units_Index(k),:),'r')
   hold on
   plot(TimeLine,SelectedResponses{1,2}(Units_Index(k),:),'k')
   hold on
   plot([Pert_lims(1)-0.1 Pert_lims(1)-0.1], [0 1],'--b','LineWidth',1);
   hold on
   plot([Pert_lims(2) Pert_lims(2)], [0 1],'--b','LineWidth',1);
   plot([0 0],[0 1],'--k','LineWidth',1.5);
   hold on
   plot([7.3 7.3], [0 1],'--k','LineWidth',1.5);
   ylim([0 1])
   xlim([-1 8])
   ll = legend({['pertON run DM= ' num2str(DM_rs(k,1))] , ['pertON still, DM=' num2str(DM_rs(k,2))]},'fontsize',10);
   ll.Color = 'none'; ll.EdgeColor = 'none';
   title(['id=' num2str(find(I==Units_Index(k))) ', RunTrials=' num2str(nansum(SelTrials{1}(:,Units_Index(k)))) ', StillTrials=' num2str(nansum(SelTrials{2}(:,Units_Index(k))))]);    
   if mod(k,16)==0
       saveas(gcf,['X:\DATA\PROJECTS\VisPerturbation\Figures\Fig_1\panels\RunStillAnalysis\MeanResponseRunStill_' num2str(fig_n) '_a.fig'])
       fig_n = fig_n + 1;
       close
   end
end



%%
if ~exist('PertResp_units','var')
    % select only perturbation responsive units
    load('AUC_shuffled.mat')
    Sh_responses = AUC_shuffled(:,2:end);
    p_pert_th = prctile(Sh_responses(:),99);
    PertResp_units = (AUC_shuffled(:,1)>p_pert_th);
    % select only pos. modulated perturbation responsive units
    load('DM_pert_shuffled.mat')
    DM = DM_sh(:,1);
    DM_sign_i(:,1) = DM>0;
    DM_sign_i(:,2) = DM<=0;
    % select only pos. modulated perturbation responsive units
    UnitsPertPos = PertResp_units & DM_sign_i(:,1);
    UnitsPertNeg = PertResp_units & DM_sign_i(:,2);
end

AUC_units = AUC_shuffled(:,1);
p_sh(2) = p_pert_th;

figure
% get this variable from AUC_units.m
plot(AUC_units(:,1),'.')
hold on
plot([0 length(AUC_units(:,1))],[p_sh(2) p_sh(2)],'k')

% rr = Units_Index(find(AUC_units(Units_Index)>p_sh(2)));
ri = find(AUC_units(:,1)>p_sh(2));
% get only positive units
ri = find(UnitsPertPos(1:length(index_i)));

clear DM_rs_i
DM_rs_i = DM_rs(ri,:);
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


%% plot only reliable units
param_sel = AM_Param(:,SelectedCells,:);
% plot only pert. resp. units from sessions with stat&run conditions
ri = find(UnitsPertPos(1:length(index_i)) & ~isnan(DM_rs(:,1)) & ~isnan(DM_rs(:,2)) & DM_rs(:,1)~=0 & DM_rs(:,2)~=0);
fig_n =1;
for k = 1:length(ri)
   if mod(k,16)==0
       subplot_nr = 16;
   elseif mod(k,16)==1
       FIG = figure('units','normalized','outerposition',[0 0 1 1]);
       subplot_nr = mod(k,16);
   else
       subplot_nr = mod(k,16);
   end
   TimeLine = linspace(-1,8.3,length(SelectedResponses{1,1}(ri(k),:)));
   Pert_lims =[param_sel(min(find(param_sel(:,ri(k),3)==1)),ri(k),4) param_sel(min(find(param_sel(:,ri(k),3)==1)),ri(k),5)]; 
   subplot(4,4,subplot_nr)
   nanmean(UnitResponses_smooth,(ri(k),:))
   
   plot(TimeLine,SelectedResponses{1,1}(ri(k),:),'r')
   hold on
   plot(TimeLine,SelectedResponses{1,2}(ri(k),:),'k')
   hold on
   plot([Pert_lims(1)-0.1 Pert_lims(1)-0.1], [0 1],'--b','LineWidth',1);
   hold on
   plot([Pert_lims(2) Pert_lims(2)], [0 1],'--b','LineWidth',1);
   plot([0 0],[0 1],'--k','LineWidth',1.5);
   hold on
   plot([7.3 7.3], [0 1],'--k','LineWidth',1.5);
   ylim([0 1])
   xlim([-1 8])
   ll = legend({['pertON run DM= ' num2str(DM_rs((ri(k)),1))] , ['pertON still, DM=' num2str(DM_rs((ri(k)),2))]},'fontsize',10);
   ll.Color = 'none'; ll.EdgeColor = 'none';
   title(['id=' num2str(find(I==ri(k))) ', RunTrials=' num2str(nansum(SelTrials{1}(:,ri(k)))) ...
           ', StillTrials=' num2str(nansum(SelTrials{2}(:,ri(k)))) ...   
            'AUC=' num2str(AUC_shuffled(ri(k),1),2) ', ' ...
            'MI=' num2str(DM_sh(ri(k),1),2)]);
   box off
   
%    if mod(k,16)==0
%        saveas(gcf,['X:\DATA\PROJECTS\VisPerturbation\Figures\Fig_1\panels\RunStillAnalysis\MeanResponseRunStill_' num2str(fig_n) '_PertResp.fig'])
%        fig_n = fig_n + 1;
% %        close
%    end
end

%% sort the cells by their response to perturbation (reference for heat map)
[B,I,PertTrial_Ex] = PlotSortingOrder(ProjectData,SelectedResponses,AM_Param,AM_Speed,Units_Sel,0);
% I = descending order taking into account:[pert response; 2s vis stim respo; post-pert reponse]
I = I(~isnan(B(:,1)));
B = reshape(B(~isnan(B(:))),length(B(~isnan(B(:))))/3,3);


%% when were responsive units recorded? during run or during stillness?
% save indeces of units sign. responsive
RespUnits = find(AUC_units(:,1)>p_sh(2));
% go throught the binary params for running/stationary info of pert period and all trial duration
% AM_SpeedResponse(:,:,1) = mean speed -0.5s to +2s at visual stimulus onset
% AM_SpeedResponse(:,:,2) = mean speed -0.5s to end of perturbation period
% AM_SpeedResponse(:,:,3) = mean speed -0.5s to +1s at perturbation offset
% AM_SpeedResponse(:,:,4) = mean speed across entire trial
clear SpeedResponse MeanSpeedResponse
SpeedResponse = squeeze(AM_SpeedResponse(:,SelectedCells,2));
SpeedResponse = SpeedResponse(:,UnitsPertPos);
MeanSpeedResponse(:,1) = nansum(SpeedResponse>Run_TH);
MeanSpeedResponse(:,2) = nansum(SpeedResponse<=Run_TH);
clear ActualMeanSpeed 
for g = 1:size(SpeedResponse,2)
    ActualMeanSpeed(g,1) = mean(SpeedResponse(SpeedResponse(:,g)>Run_TH,g));
    ActualMeanSpeed(g,2) = mean(SpeedResponse(SpeedResponse(:,g)<=Run_TH,g));
end

figure
area(MeanSpeedResponse./sum(MeanSpeedResponse,2)*100)
xlabel('Sign.responsive units'); ylabel('% trials with pert ON');
box off
set(gca,'TickDir','out','FontSize',13,'FontName','Arial')
ll = legend({['>' num2str(Run_TH) 'cm/s'],['<=' num2str(Run_TH) 'cm/s']},'fontsize',10);
ll.Color = 'none'; ll.EdgeColor = 'none';
hold on
V_percent = unique(MeanSpeedResponse(:,1))
for h = 1:length(V_percent)
    x_loc(h) = min(find(MeanSpeedResponse(:,1)==V_percent(h)));
end
x_loc = sort(x_loc);
for h = 1:length(V_percent)
    h = text(x_loc(h),mod(h,2)*5+5,num2str(ActualMeanSpeed(x_loc(h),1),'%.1f'))
    set(h,'Rotation',90);
end
title(['DMI_r_u_n / DMI_s_t_i_l_l'])

% total nr of trials running vs still
RunTrialMean = mean(MeanSpeedResponse./sum(MeanSpeedResponse,2)*100);
% nr of recordings in which we found sign. responsive units
RecNrs = sum(diff(MeanSpeedResponse(:,1))~=0)+1;


% sign. responsive units
x = MeanSpeedResponse./sum(MeanSpeedResponse,2)*100;
% not. sign. responsive units
TrialsType = [TotTrials{2,1}(:,SelectedCells);...   % pert ON (running; still)
              TotTrials{2,2}(:,SelectedCells) ]'; % pert OFF (running; still)
ri = find(AUC_units(:,1)<p_sh(2)); % non sign. pert. responsive units
x_n = TrialsType(ri,1:2)./sum(TrialsType(ri,1:2),2)*100;
bins = 0:1:100;
dist = histcounts(x(:,2), bins); cdist = cumsum(dist/sum(dist));
dist_n = histcounts(x_n(:,2), bins); cdist_n = cumsum(dist_n/sum(dist_n));

figure
plot(cdist)
hold on
plot(cdist_n)
box off
set(gca,'TickDir','out')    
ll = legend({'Sign. responsive units', 'Other units'},'location','northwest','fontsize',10)
ll.Color = 'none'; ll.EdgeColor = 'none'; 
set(gca,'FontSize',13);
set(gca,'fontsize',13,'TickDir','out');
xlabel('% stationary trials','fontsize',13); ylabel('Cumulative probability','fontsize',13);

% take sessions with pert. units.
% proportion of units responsive to perturbation
RecsWPert = unique(param_sel(1,UnitsPertPos,1));
clear TotalUnits_rec PertUnits_rec
for i = 1:length(RecsWPert)
   TotalUnits_rec(i) = sum(param_sel(1,:,1) == RecsWPert(i));
   PertUnits_rec(i) = sum(param_sel(1,UnitsPertPos,1) == RecsWPert(i));
end
PertUnitPercent = PertUnits_rec./TotalUnits_rec*100;
[b , ind] = sort(PertUnitPercent);
% take into account the nr of trials running
TrialsPercent = MeanSpeedResponse./sum(MeanSpeedResponse,2)*100;
clear RunPercent
for i = 1:length(RecsWPert)
   RunPercent(i) = TrialsPercent( min(find(param_sel(1,UnitsPertPos,1) == RecsWPert(i),1)) );
end

figure
plot(100-RunPercent,PertUnitPercent,'.k')
hold on
plot(100-RunPercent,PertUnitPercent,'or')
hold on
% bins = 0:1:100;
% dist = histcounts(PertUnits_rec/sum(PertUnits_rec)*100, bins); cdist = cumsum(dist/sum(dist));
[c, ind] = sort(100-RunPercent);
plot([100-RunPercent(ind) 100],[cumsum(PertUnits_rec/sum(PertUnits_rec)*100) 100],'r')
xlim([0 100]); ylim([0 100])
ylabel('% of perturbation responsive units')
xlabel('% of trials @ 2 < cm/s')
box off
set(gca,'TickDir','out')

% sign. responsive units
x = MeanSpeedResponse./sum(MeanSpeedResponse,2)*100;
% not. sign. responsive units
TrialsType = [TotTrials{2,1}(:,SelectedCells);...   % pert ON (running; still)
              TotTrials{2,2}(:,SelectedCells) ]'; % pert OFF (running; still)
% visually responsive units
DM_vis = squeeze(FR_UOI_VIS(3,:,SelectedCells));
DM_vis_sh = [prctile(DM_vis(2:end,:),1); prctile(DM_vis(2:end,:),99)];
ri = DM_vis(1,:)<=DM_vis_sh(1,:) | DM_vis(1,:)>=DM_vis_sh(2,:);
x_n = TrialsType(ri,1:2)./sum(TrialsType(ri,1:2),2)*100;
bins = 0:1:100;
dist = histcounts(x(:,2), bins); cdist = cumsum(dist/sum(dist));
dist_n = histcounts(x_n(:,2), bins); cdist_n = cumsum(dist_n/sum(dist_n));

figure
plot(cdist)
hold on
plot(cdist_n)
box off
set(gca,'TickDir','out')    
ll = legend({'Sign. responsive units', 'Visually responsive units'},'location','northwest','fontsize',10)
ll.Color = 'none'; ll.EdgeColor = 'none'; 
set(gca,'FontSize',13);
set(gca,'fontsize',13,'TickDir','out');
xlabel('% stationary trials','fontsize',13); ylabel('Cumulative probability','fontsize',13);



% plot box plot showing change of magnitude of response in run/still 
figure
for i = 1:length(DM_rs_i)
    plot([DM_rs_i(i,1) DM_rs_i(i,2)],'Color',[0.5 0.5 0.5])
    hold on
end
xlim([0 3]); ylim([-1 1]);
set(gca,'XTick',[1 2],'XTicklabel',{'run','still'},'TickDir','out')
ylabel('Modulation index'); box off;
plot([mean(DM_rs_i(:,1)) mean(DM_rs_i(:,2))],'k','LineWidth',3)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot mean responses during RUN of Pert ON and Pert OFF, and during STILL of Pert ON and Pert OFF

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
UnitResponses_smooth = AM_UnitResponses_smooth(:,SelectedCells,:) ;

%% select responses for running and still trials
% select trials of interest and control trials as well
Trials_PertON = AM_Param(:,:,3)==1; % find indexes where perturbation is on
Trials_PertOFF= AM_Param(:,:,3)==0; % find indexes where perturbation is off
Trials_PertOFF_lim = logical(zeros(size(Trials_PertOFF,1),size(Trials_PertOFF,2)));
SelectedCells = AM_UOI;
[AM_SpeedResponse AM_SpeedResponseControl] = ComputeSpeedResponse(AM_Param, AM_Speed,SelectedCells,Trials_PertON,Trials_PertOFF,ProjectData);
% AM_SpeedResponse(:,:,1) = mean speed -0.5s to +2s at visual stimulus onset
% AM_SpeedResponse(:,:,2) = mean speed -0.5s to end of perturbation period
% AM_SpeedResponse(:,:,3) = mean speed -0.5s to +1s at perturbation offset
% AM_SpeedResponse(:,:,4) = mean speed across entire trial

[SelectedTrials_OI  SelectedTrials_OI_control TotTrials] = SpeedFiltering(AM_Speed, AM_SpeedResponse, AM_SpeedResponseControl, Trials_PertON,Trials_PertOFF, Run_TH);
% SelectedTrials_OI{4,2} = 4 rows as above, column 1 = pert ON & run; column 2 = pert ON & static
% SelectedTrials_OI_control{4,2} = 4 rows as above, column 1 = pert OFF & run; column 2 = pert OFF & static


% plot distribution of trials based on run/still conditions
% 1 consider the perturbation trials
TrialsType = [TotTrials{2,1}(:,SelectedCells);...   % pert ON (running; still)
              TotTrials{2,2}(:,SelectedCells) ]'; % pert OFF (running; still)

ri = find(AUC_units(:,1)>p_sh(2));

figure
subplot(2,1,1)
area(TrialsType(ri,1:2)./sum(TrialsType(ri,1:2),2)*100)
hold on
plot(ri,zeros(length(ri),1),'^k')
set(gca,'XTickLabels',[], 'TickDir','out');
box off; ylabel('Pert ON trials %'); 
xlim([0 size(TrialsType,1)])
ll = legend({['Pert ON >' num2str(Run_TH) 'cm/s'],['Pert ON <=' num2str(Run_TH) 'cm/s'],'Pert. responsive units'},'fontsize',10);
ll.Color = 'none'; ll.EdgeColor = 'none';

subplot(2,1,2)
area(TrialsType(:,3:4)./sum(TrialsType(:,3:4),2)*100)
hold on
plot(ri,zeros(length(ri),1),'vk')
set(gca, 'TickDir','out'); set(gca, 'YDir','reverse')
box off; ylabel('Pert OFF trials %'); xlabel('units');
xlim([0 size(TrialsType,1)]);
ll = legend({['Pert OFF >' num2str(Run_TH) 'cm/s'],['Pert OFF <=' num2str(Run_TH) 'cm/s'],'Pert. responsive units'},'fontsize',10);
ll.Color = 'none'; ll.EdgeColor = 'none';



% generate 3D matrix for including time responses
t1 = double(repmat(SelectedTrials_OI{2,1}, 1, 1, size(AM_UnitResponses_smooth,3))); % repeating the selection across time
t2 =  double(repmat(SelectedTrials_OI{2,2}, 1, 1, size(AM_UnitResponses_smooth,3))); % repeating the selection across time
% nan all the zero values
t1(t1(:)==0) = nan; t2(t2(:)==0) = nan;

% select only the responses during selected trials (pert ON and OFF separately)
TrialsResp_OI = t1.*AM_UnitResponses_smooth;
TrialsRespControl_OI = t2.*AM_UnitResponses_smooth;
% compute the mean responses for each unit during trials with pert ON
% and OFF
TrialsResp_OI_2D = squeeze(nanmean(TrialsResp_OI,1));
TrialsRespControl_OI_2D = squeeze(nanmean(TrialsRespControl_OI,1));
% normalise the responses of interest
TrialsResp_max = max(nanmax(TrialsResp_OI_2D,[],2),nanmax(TrialsRespControl_OI_2D,[],2));
Response_OI = TrialsResp_OI_2D./TrialsResp_max;
ResponseControl_OI = TrialsRespControl_OI_2D./TrialsResp_max;

% generate 3D matrix for including time responses
t3 = double(repmat(SelectedTrials_OI_control{2,1}, 1, 1, size(AM_UnitResponses_smooth,3))); % repeating the selection across time
t4 =  double(repmat(SelectedTrials_OI_control{2,2}, 1, 1, size(AM_UnitResponses_smooth,3))); % repeating the selection across time
% nan all the zero values
t3(t3(:)==0) = nan; t3(t4(:)==0) = nan;

%% plot modulation indexes
% consider only running during perturbation period
Param = squeeze(AM_Param(:,SelectedCells,3)); % all units
Param_run = SelectedTrials_OI{2,1}(:,SelectedCells); 
Param_still = SelectedTrials_OI{2,2}(:,SelectedCells); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Units_Index = find(DM_rs(:,1)~=0 & DM_rs(:,2)~=0);
Units_Index = find(nansum(SelTrials{1,1},1));
Units_Index_still = find(nansum(SelTrials{1,2},1));
SelTrials{1} = SelectedTrials_OI{2,1}(:,SelectedCells);
SelTrials{2} = SelectedTrials_OI{2,2}(:,SelectedCells);
fig_n=1;
for k = 1:size(Units_Index,2)
   if mod(k,16)==0
       subplot_nr = 16;
   elseif mod(k,16)==1
       FIG = figure('units','normalized','outerposition',[0 0 1 1]);
       subplot_nr = mod(k,16);
   else
       subplot_nr = mod(k,16);
   end
   TimeLine = linspace(-1,8.3,length(SelectedResponses{1,1}(Units_Index(k),:)));
   Pert_lims =[param_sel(min(find(param_sel(:,Units_Index(k),3)==1)),Units_Index(k),4) param_sel(min(find(param_sel(:,Units_Index(k),3)==1)),Units_Index(k),5)]; 
   subplot(4,4,subplot_nr)
   plot(TimeLine,SelectedResponses{1,1}(Units_Index(k),:),'r')
   hold on
   plot(TimeLine,SelectedResponses{1,2}(Units_Index(k),:),'k')
   hold on
   plot([Pert_lims(1)-0.1 Pert_lims(1)-0.1], [0 1],'--b','LineWidth',1);
   hold on
   plot([Pert_lims(2) Pert_lims(2)], [0 1],'--b','LineWidth',1);
   plot([0 0],[0 1],'--k','LineWidth',1.5);
   hold on
   plot([7.3 7.3], [0 1],'--k','LineWidth',1.5);
   ylim([0 1])
   xlim([-1 8])
   ll = legend({['pertON run DM= ' num2str(DM_rs(k,1))] , ['pertON still, DM=' num2str(DM_rs(k,2))]},'fontsize',10);
   ll.Color = 'none'; ll.EdgeColor = 'none';
   title(['id=' num2str(find(I==Units_Index(k))) ', RunTrials=' num2str(nansum(SelTrials{1}(:,Units_Index(k)))) ', StillTrials=' num2str(nansum(SelTrials{2}(:,Units_Index(k))))]);    
   if mod(k,16)==0
       saveas(gcf,['X:\DATA\PROJECTS\VisPerturbation\Figures\Fig_1\panels\RunStillAnalysis\MeanResponseRunStill_' num2str(fig_n) '_a.fig'])
       fig_n = fig_n + 1;
       close
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


