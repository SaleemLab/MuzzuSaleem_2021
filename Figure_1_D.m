%% Tomaso Muzzu - UCL - 30/03/2020


%% Plot example units

%% Functions to plot Figure 2 - direction tuning
if ~exist('ProjectData','var')
    [ProjectData AM_UnitResponses AM_Param AM_Speed AM_UOI SelectedResponses AM_UnitResponses_smooth] = LoadDataALL;
end
thres = 95;
% first 7 animals
if  size(ProjectData,1) == 37
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
        %DM = DM_sh(:,1);
        DM_sign_i(:,1) = DM>0;
        DM_sign_i(:,2) = DM<=0;
        % select only pos. modulated perturbation responsive units
        PertRespUnits_pos = PertResp_units & DM_sign_i(:,1);
        PertRespUnits_neg = PertResp_units & DM_sign_i(:,2);
    end
    
elseif size(ProjectData,1) == 10
    CTRL_exp = 1;
    % naive animals
    Animal_1st_idx = [1 4 7];
    if ~exist('PertResp_units','var')
        % select only perturbation responsive units
        load('AUC_shuffled_CTRL_1.mat')
        Sh_responses = AUC_shuffled(:,2:end);
        p_pert_th = prctile(Sh_responses(:),thres);
        PertResp_units = (AUC_shuffled(:,1)>p_pert_th);
        % select only pos. modulated perturbation responsive units
        load('DM_CTRL.mat')
        %DM = DM_sh(:,1);
        DM_sign_i(:,1) = DM>0;
        DM_sign_i(:,2) = DM<=0;
        % select only pos. modulated perturbation responsive units
        PertRespUnits_pos = PertResp_units & DM_sign_i(:,1);
        PertRespUnits_neg = PertResp_units & DM_sign_i(:,2);
    end
    
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
param_sel = AM_Param(:,SelectedCells,:);
UnitResponses_smooth = AM_UnitResponses_smooth(:,SelectedCells,:);


%% load AUC values of the units from control group
% load('AUC_shuffled_CTRL_1.mat')
AUC_units = AUC_shuffled(:,1);
p_sh = prctile(reshape(AUC_shuffled(:,2:end),1,size(AUC_shuffled(:,2:end),1)*size(AUC_shuffled(:,2:end),2)),thres);


%% compute modulation index for all units
param_sel = AM_Param(:,SelectedCells,:);
Param = squeeze(param_sel(:,:,3));
% % Depth of modulation
r= 1; clear DM DepthMod
for i = 1:size(Param,2)
    % select trial responses
    Responses = Param(:,i);
    % scrub responses
    Responses = logical(Responses(~isnan(Responses)));
    UnitResponse = mean(squeeze(UnitResponses_smooth(Responses==1,i,:)))*60;
    %UnitResponse = UnitResponse - mean(UnitResponse(po_i(1)-trialSide_samples:po_i(1)));
%     Dir = sign(sum(UnitResponse(po_i(1):po_i(1)+trialSide_samples))-sum(UnitResponse(po_i(1)-trialSide_samples:po_i(1))));
%     switch Dir
%         case 1
            DM(i) = (nanmean(UnitResponse(po_i(1):po_i(1)+trialSide_samples))-nanmean(UnitResponse(po_i(1)-trialSide_samples:po_i(1))))/...
                nanmean(UnitResponse(po_i(1)-trialSide_samples:po_i(1)));
%             DM(i) = (nanmax(UnitResponse(po_i(1):po_i(1)+trialSide_samples))-nanmax(UnitResponse(po_i(1)-trialSide_samples:po_i(1))))/...
%                 nanmax(UnitResponse(po_i(1)-trialSide_samples:po_i(1)));
           
%             DM(i) = (nanmean(UnitResponse(po_i(1):po_i(1)+trialSide_samples))-nanmean(UnitResponse(po_i(1)-trialSide_samples:po_i(1))))/...
%                 (nanmean(UnitResponse(po_i(1):po_i(1)+trialSide_samples))+nanmean(UnitResponse(po_i(1)-trialSide_samples:po_i(1))));
%         case -1
%             DM(i) = (min(UnitResponse(po_i(1):po_i(1)+trialSide_samples))-min(UnitResponse(po_i(1)-trialSide_samples:po_i(1))))/...
%                 (min(UnitResponse(po_i(1):po_i(1)+trialSide_samples))+min(UnitResponse(po_i(1)-trialSide_samples:po_i(1))));
%     end
    
    if AUC_units(i,1)>p_sh
        DepthMod(r,:) = [DM(i) AUC_units(i) i];
        r = r + 1;
    end
    
%     parfor sh=1:1000
%         % shuffle perturbation trals
%         Responses_sh = Responses(randperm(length(Responses)));
%         UnitResponse = mean(squeeze(UnitResponses_smooth(Responses_sh==1,i,:)))*60;
%         DM_sh(i,sh) = (sum(UnitResponse(po_i(1):po_i(1)+trialSide_samples))-sum(UnitResponse(po_i(1)-trialSide_samples:po_i(1))))/...
%             (sum(UnitResponse(po_i(1):po_i(1)+trialSide_samples))+sum(UnitResponse(po_i(1)-trialSide_samples:po_i(1))));
%     end

end

% verify disttribution of shuffled DM's
% figure
%  histogram(DM_sh,[-1:0.05:1])
% save('DM_ALL.mat','DM')

[unit_nr(:,1) unit_nr(:,2)] = sort(DM);
161 158 185

figure
DM_h = DM; DM_h(DepthMod(:,3))=nan; DM_h(isnan(DM_h))=[];
h = histogram(DM_h,[-1:0.05:1],'FaceColor',[0.5 0.5 0.5],'Normalization','probability','EdgeColor','none')
hold on
h = histogram(DepthMod(:,1),[-1:0.05:1],'FaceColor',[1 0 0],'Normalization','probability','EdgeColor','none')
hold on
xlabel('Depth of Modulation Index'); ylabel('probability')
box off
set(gca,'FontSize',13,'TickDir','out');
text(-0.5,0.1,['n=' num2str(sum(DepthMod(:,1)<=0))],'FontSize',13)
text(0.5,0.1,['n=' num2str(sum(DepthMod(:,1)>0))],'FontSize',13)
ll = legend({['Not sign. responsive units n=' num2str(length(DM)-size(DepthMod,1))],['Sign. responsive units n=' num2str(size(DepthMod,1))]},'fontsize',10);
ll.Color = 'none'; ll.EdgeColor = 'none';

[h,p1]=kstest2(DepthMod(DepthMod(:,1)>0,1),DM_h(DM_h>0))
[h,p2]=kstest2(DepthMod(DepthMod(:,1)<=0,1),DM_h(DM_h<=0))
[h,p2]=kstest2(DepthMod(:,1),DM_h)

figure
plot(-1+0.05:0.05:1,cumsum(histcounts(DM_h,-1:0.05:1))/(sum(histcounts(DM_h,-1:0.05:1)))*100,'r')
hold on
plot(-1+0.05:0.05:1,cumsum(histcounts(DepthMod(:,1),-1:0.05:1))/(sum(histcounts(DepthMod(:,1),-1:0.05:1)))*100,'k')
box off
ylabel('% units')
xlabel('Modulation Index')
set(gca,'FontSize',13,'TickDir','out');
ll = legend({['Not sign. responsive units n=' num2str(length(DM)-size(DepthMod,1))],['Sign. responsive units n=' num2str(size(DepthMod,1))]},'fontsize',10);
ll.Color = 'none'; ll.EdgeColor = 'none';
text(0.5,75,['p=' num2str(p2,1)],'FontSize',13)


%% show shuffling process
figure
h = histogram(AUC_shuffled(:,2:end),[0:0.02:1],'FaceColor',[0.5 0.5 0.5],'Normalization','probability','EdgeColor','none')
hold on
h1 = histogram(AUC_units(:,1),[0:0.02:1],'FaceColor',[1 0 0],'Normalization','probability','EdgeColor','none')
plot([p_sh(1) p_sh(1)], [0 max(h.BinCounts/length(AUC_shuffled(:)))*1.1],'k') 
%plot([p_sh(2) p_sh(2)], [0 max(h1.BinCounts/length(AUC_units(:,1)))*1.1],'k--','LineWidth',3)
xlabel('AUC of shuffled data'); ylabel('probability')
ll = legend({['AUC of shuffled data n=' num2str(sum(AM_UOI)) 'k'],['AUC of units n=' num2str(sum(AM_UOI))], [num2str(thres) '% shuffled data']})
ll.Color = 'none'; ll.EdgeColor = 'none';
box off
set(gca,'FontSize',13,'TickDir','out');

%% EXAMPLE UNITS
% plot units for which you showing the AUC

%% sort the cells by their response to perturbation (reference for heat map)
Units_Sel = AM_UOI;
Diff = 0;
[B,I,PertTrial_Ex] = PlotSortingOrder(ProjectData,SelectedResponses,AM_Param,AM_Speed,Units_Sel,Diff);
% I = descending order taking into account:[pert response; 2s vis stim respo; post-pert reponse]
I = I(~isnan(B(:,1)));
B = reshape(B(~isnan(B(:))),length(B(~isnan(B(:))))/3,3);

%%
figure
param_sel = AM_Param(:,SelectedCells,:); UI = 20;
Pert_lims =[param_sel(min(find(param_sel(:,UI,3)==1)),UI,4) param_sel(min(find(param_sel(:,UI,3)==1)),UI,5)]; 
TimeLine = linspace(-1,8.3,size(UnitResponses_smooth,3));
shadedErrorBar(TimeLine,...
               smooth(mean(squeeze(UnitResponses_smooth(param_sel(:,UI,3)==1,UI,:)))*60,5),...
               std(squeeze(UnitResponses_smooth(param_sel(:,UI,3)==1,UI,:)))/sqrt(size(UnitResponses_smooth(param_sel(:,UI,3)==1,UI,:),1))*60,...
               'lineprops',{'r-','markerfacecolor','r'});
hold on
shadedErrorBar(TimeLine,...
               smooth(mean(squeeze(UnitResponses_smooth(param_sel(:,UI,3)==0,UI,:)))*60,5),...
               std(squeeze(UnitResponses_smooth(param_sel(:,UI,3)==0,UI,:)))/sqrt(size(UnitResponses_smooth(param_sel(:,UI,3)==0,UI,:),1))*60, ...
               'lineprop',{'k-','markerfacecolor','k'});
hold on
plot([0 0],[0 max(mean(squeeze(UnitResponses_smooth(param_sel(:,UI,3)==1,UI,:)))*60)],'--k','LineWidth',1.5);
hold on
plot([7.3 7.3 ], [0 max(mean(squeeze(UnitResponses_smooth(param_sel(:,UI,3)==1,UI,:)))*60)],'--k','LineWidth',1.5);
hold on
plot([Pert_lims(1)-0.1 Pert_lims(1)-0.1], [0 max(mean(squeeze(UnitResponses_smooth(param_sel(:,UI,3)==1,UI,:)))*60)],'--b','LineWidth',1);
hold on
plot([Pert_lims(2) Pert_lims(2)], [0 max(mean(squeeze(UnitResponses_smooth(param_sel(:,UI,3)==1,UI,:)))*60)],'--b','LineWidth',1);
set(gca,'FontSize',13);
ylabel('Hz','fontsize',13); xlabel('seconds')
ylim([0 max(mean(squeeze(UnitResponses_smooth(param_sel(:,UI,3)==1,UI,:)))*60)*1.3])
xlim([-1 8])
box off
set(gca,'TickDir','out');
ll = legend({'perturbation ON', 'perturbation OFF'},'fontsize',10);
ll.Color = 'none'; ll.EdgeColor = 'none';
title(['Unit = ' num2str(find(I==UI)) ', AUC = ' num2str(AUC_units(UI,1)) ',DM = ' num2str(DM(UI))])           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
     
%% percentages of responsive units per animal per group (Experienced VS naive)
% param_sel :
% conditions = 1 --> nr of recording
% conditions = 2 --> grating direction [0:45:315]
% conditions = 3 --> 1 pert ON, 0 pert OFF
% conditions = 4 --> pert onset
% conditions = 5 --> pert offset
% select only perturbation responsive units
Sh_responses_CTRL = AUC_shuffled_CTRL(:,2:end);
p_pert_th_CTRL = prctile(Sh_responses_CTRL(:),thres);
PertResp_units_CTRL = (AUC_shuffled_CTRL(:,1)>p_pert_th_CTRL);
% select only pos. modulated perturbation responsive units
DM_All = DM;
load('DM_CTRL.mat')
%DM = DM_sh(:,1);
DM_sign_i_CTRL(:,1) = DM>0;
DM_sign_i_CTRL(:,2) = DM<=0;
% select only pos. modulated perturbation responsive units
PertRespUnits_pos_CTRL = PertResp_units_CTRL & DM_sign_i_CTRL(:,1);
PertRespUnits_neg_CTRL = PertResp_units_CTRL & DM_sign_i_CTRL(:,2);

% reunite all with CTRL data
PertRespUnits_neg_n = [PertRespUnits_neg(1:830); PertRespUnits_neg_CTRL];
PertRespUnits_pos_n = [PertRespUnits_pos(1:830); PertRespUnits_pos_CTRL];

Rec_OI = param_sel(1,PertResp_units,1);
clear UnitsPerRecording
for i = 1:max(param_sel(1,:,1))
    UnitsPerRecording(i,1) = sum(param_sel(1,:,1)==i);
    UnitsPerRecording(i,2) = sum(param_sel(1,PertRespUnits_neg_n,1)==i);
    UnitsPerRecording(i,3) = sum(param_sel(1,PertRespUnits_pos_n,1)==i);
end
Mice = unique(ProjectData.Mouse_name);
for i = 1:max(param_sel(1,:,1))
    for j = 1:length(Mice)
        if strcmp(char(ProjectData.Mouse_name(i)),char(Mice(j)))
            tempMice(i) = j;
        end
    end
end
UnitsPerRecording(:,4) = tempMice';
clear RespUnitPerRec
for i = 1:length(Mice)
    idx = UnitsPerRecording(:,4) == i;
    RespUnitPerRec(i,1) = sum(sum(UnitsPerRecording(idx,2:3)))/sum(UnitsPerRecording(idx,1))*100; % all
    RespUnitPerRec(i,2) = sum(sum(UnitsPerRecording(idx,2)))/sum(UnitsPerRecording(idx,1))*100; % negative
    RespUnitPerRec(i,3) = sum(sum(UnitsPerRecording(idx,3)))/sum(UnitsPerRecording(idx,1))*100; % positive
end

figure
subplot(2,3,1) % all units
bar(1,mean(RespUnitPerRec(1:7,1)))
hold on
bar(2,mean(RespUnitPerRec(8:10,1)))
Data2Scatter = [RespUnitPerRec(1:7,1) [RespUnitPerRec(8:10,1); nan(4,1)]];
UnivarScatter(Data2Scatter,'Label',{'Experienced','Naive'},'MarkerFaceColor',[ 0.5 0.5 0.5])
title('All units')
ylabel('%'); ylim([0 100])
subplot(2,3,2) % positive units
bar(1,mean(RespUnitPerRec(1:7,3)))
hold on
bar(2,mean(RespUnitPerRec(8:10,3)))
Data2Scatter = [RespUnitPerRec(1:7,3) [RespUnitPerRec(8:10,3); nan(4,1)]];
UnivarScatter(Data2Scatter,'Label',{'Experienced','Naive'},'MarkerFaceColor',[ 0.5 0.5 0.5])
title('MI>0')
ylabel('%'); ylim([0 100])
subplot(2,3,3) % negative units
bar(1,mean(RespUnitPerRec(1:7,2)))
hold on
bar(2,mean(RespUnitPerRec(8:10,2)))
Data2Scatter = [RespUnitPerRec(1:7,2) [RespUnitPerRec(8:10,2); nan(4,1)]];
UnivarScatter(Data2Scatter,'Label',{'Experienced','Naive'},'MarkerFaceColor',[ 0.5 0.5 0.5])
title('MI<0')
% plot recordings only
Percent2plot(:,1) = ((UnitsPerRecording(:,2)+UnitsPerRecording(:,3))./UnitsPerRecording(:,1))*100;
Percent2plot(:,2) = ((UnitsPerRecording(:,2))./UnitsPerRecording(:,1))*100;
Percent2plot(:,3) = ((UnitsPerRecording(:,3))./UnitsPerRecording(:,1))*100;
ylabel('%'); ylim([0 100])
subplot(2,3,4) % all units
bar(1,mean(Percent2plot(1:37,1)))
hold on
bar(2,mean(Percent2plot(38:end,1)))
Data2Scatter = [Percent2plot(1:37,1) [Percent2plot(38:end,1); nan(27,1)]];
UnivarScatter(Data2Scatter,'Label',{'Experienced','Naive'},'MarkerFaceColor',[ 0.5 0.5 0.5])
title('All units')
ylabel('%'); ylim([0 100])
subplot(2,3,5) % positive units
bar(1,mean(Percent2plot(1:37,3)))
hold on
bar(2,mean(Percent2plot(38:end,3)))
Data2Scatter = [Percent2plot(1:37,3) [Percent2plot(38:end,3); nan(27,1)]];
UnivarScatter(Data2Scatter,'Label',{'Experienced','Naive'},'MarkerFaceColor',[ 0.5 0.5 0.5])
title('MI>0')
subplot(2,3,6) % negative units
bar(1,mean(Percent2plot(1:37,2)))
hold on
bar(2,mean(Percent2plot(38:end,2)))
Data2Scatter = [Percent2plot(1:37,2) [Percent2plot(38:end,2); nan(27,1)]];
UnivarScatter(Data2Scatter,'Label',{'Experienced','Naive'},'MarkerFaceColor',[ 0.5 0.5 0.5])
title('MI<0')
ylabel('%'); ylim([0 100])






