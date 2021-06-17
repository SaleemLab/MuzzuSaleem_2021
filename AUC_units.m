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
UnitResponses_smooth = AM_UnitResponses_smooth(:,SelectedCells,:);
clear DiffResponsePercent_Pert DiffResponsePercent_PostPert DiffResponsePercent_Vis
for i = 1:size(UnitResponses_smooth,2)
    
    UnitResponse = squeeze(UnitResponses_smooth(:,i,:))./max(squeeze(UnitResponses_smooth(:,i,:)),[],2);
    
    DiffResponsePercent_Pert(:,i,1) = (nanmean(UnitResponse(:,po_i(1):po_i(2)),2)./nanmean(UnitResponse(:,po_i(1)-trialSide_samples:po_i(1)),2)); % mean(pert)/mean(prepert)
    DiffResponsePercent_Pert(:,i,2) = (nanmean(UnitResponse(:,po_i(1):po_i(2)),2)-nanmean(UnitResponse(:,po_i(1)-trialSide_samples:po_i(1)),2)); % mean(pert)-mean(prepert)
    DiffResponsePercent_Pert(:,i,3) = (nansum(UnitResponse(:,po_i(1):po_i(2)),2)); % integral of mean(pert)
    DiffResponsePercent_Pert(:,i,4) = (nansum(UnitResponse(:,po_i(1)-trialSide_samples:po_i(2)),2)); % integral of mean(pert) and preceding 1s
    DiffResponsePercent_Pert(:,i,5) = (nanmean(UnitResponse(:,po_i(1):po_i(2)-trialSide_samples/2),2)-nanmean(UnitResponse(:,po_i(1)-trialSide_samples:po_i(1)),2))./...
                                      (nanmean(UnitResponse(:,po_i(1):po_i(2)-trialSide_samples/2),2)+nanmean(UnitResponse(:,po_i(1)-trialSide_samples:po_i(1)),2)); % mean(pert)-mean(prepert)/mean(pert+max(prepert) a.k.a depth of modulation
    DiffResponsePercent_Pert(:,i,6) = nanmean(UnitResponse(:,po_i(1):po_i(2)-trialSide_samples/2),2)./...
                                      (nanmean(UnitResponse(:,po_i(1):po_i(2)-trialSide_samples/2),2)+nanmean(UnitResponse(:,po_i(1)-trialSide_samples:po_i(1)),2)); % mean(pert)/mean(pert)+mean(prepert) a.k.a modulation index
    
    DiffResponsePercent_PostPert(:,i,1) = (nanmean(UnitResponse(:,po_i(2):po_i(2)+trialSide_samples),2)./nanmean(UnitResponse(:,po_i(1):po_i(2)),2)); % mean(postpert)/mean(epert)
    DiffResponsePercent_PostPert(:,i,2) = (nanmean(UnitResponse(:,po_i(2):po_i(2)+trialSide_samples),2)-nanmean(UnitResponse(:,po_i(1):po_i(2)),2)); % mean(postpert)-mean(pert)
    DiffResponsePercent_PostPert(:,i,3) = (nansum(UnitResponse(:,po_i(2):po_i(2)+trialSide_samples),2)); % integral of mean(postpert)
    DiffResponsePercent_PostPert(:,i,4) = (nansum(UnitResponse(:,po_i(1):po_i(2)+trialSide_samples),2)); % integral of mean(postpert) and preceding 1s
    DiffResponsePercent_PostPert(:,i,5) = (nanmean(UnitResponse(:,po_i(2):po_i(2)+trialSide_samples/2),2)-nanmean(UnitResponse(:,po_i(2)-trialSide_samples/2:po_i(2)),2))./...
                                          (nanmean(UnitResponse(:,po_i(2):po_i(2)+trialSide_samples/2),2)+nanmean(UnitResponse(:,po_i(2)-trialSide_samples/2:po_i(2)),2)); % mean(pert)-mean(prepert)/mean(pert)+mean(prepert) a.k.a depth of modulation
    DiffResponsePercent_PostPert(:,i,6) =  nanmean(UnitResponse(:,po_i(2):po_i(2)+trialSide_samples/2),2)./...
                                          (nanmean(UnitResponse(:,po_i(2):po_i(2)+trialSide_samples/2),2)+nanmean(UnitResponse(:,po_i(2)-trialSide_samples/2:po_i(2)),2)); % mean(pert)/mean(pert)+mean(prepert) a.k.a modulation index
%     
%   DiffResponsePercent_Vis(:,i) = (nanmean(UnitResponse(:,trialSide_samples:trialSide_samples+2*BonVision_SR),2)./nanmean(UnitResponse(:,1:trialSide_samples),2)).^2;

end

%% step 2: use these values with the actual response to fit a logistic model
Param = squeeze(AM_Param(:,SelectedCells,3));
clear SingleUnitResponse AUC_units
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
    AUC_units(i,1) = AUC1;
%     if AUC1>0.8 && plotted==0
%         hold on
%         plot(X,Y,'Color',[0.5 0.5 0.5],'LineWidth',2)
%         xlabel('False positive rate')
%         ylabel('True positive rate')
%         title('ROC for Classification by Logistic Regression')
%         saved = [AUC1 i];
%         plotted = 1;
%     end
    
    
    %% post perturbation period 
    SingleUnitResponse = squeeze(DiffResponsePercent_PostPert(:,i,:));
    % scrub predictor
    SingleUnitResponse = SingleUnitResponse(1:length(Responses),:);
    SingleUnitResponse(SingleUnitResponse==Inf,:) = 0;
    SingleUnitResponse(isnan(SingleUnitResponse(:))) = 0;
    % apply logistic model
    mdl = fitglm(SingleUnitResponse,Responses,'Distribution','binomial','Link','logit');
    scores = mdl.Fitted.Probability;
    % step 3: compute the ROC and evaluate the AUC for each unit
    [X,Y,T,AUC2] = perfcurve(Responses,scores,1);
    AUC_units(i,2) = AUC2;
    
end


%% step 4: find cutoff value through random shuffling of Pert trials id
tic
clear AUC_shuffled_pp AUC_shuffled
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
    
    % post perturbation period
    SingleUnitResponse_pp = squeeze(DiffResponsePercent_PostPert(:,i,:));
    % scrub predictor
    SingleUnitResponse_pp = SingleUnitResponse_pp(1:length(Responses),:);
    SingleUnitResponse_pp(SingleUnitResponse_pp==Inf,:) = 0;
    SingleUnitResponse_pp(isnan(SingleUnitResponse_pp(:))) = 0;
    for sh_i = 1:1000
        % shuffle perturbation trals
        Responses_sh = Responses(randperm(length(Responses)));
        
        % apply logistic model
        mdl = fitglm(SingleUnitResponse,Responses_sh,'Distribution','binomial','Link','logit');
        scores = mdl.Fitted.Probability;
        % step 3: compute the ROC and evaluate the AUC for each unit
        [X,Y,T,AUC1] = perfcurve(Responses,scores,1);
        AUC_shuffled(i,sh_i) = AUC1;
        
        % post perturbation period
        mdl = fitglm(SingleUnitResponse_pp,Responses_sh,'Distribution','binomial','Link','logit');
        scores = mdl.Fitted.Probability;
        % step 3: compute the ROC and evaluate the AUC for each unit
        [X,Y,T,AUC2] = perfcurve(Responses,scores,1);
        AUC_shuffled_pp(i,sh_i) = AUC2;
    end
    i
end
toc

AUC_shuffled = cat(2,AUC_units(:,1),AUC_shuffled);
AUC_shuffled_pp = cat(2,AUC_units(:,2),AUC_shuffled_pp);

save('AUC_shuffled_CTRL.mat','AUC_shuffled','AUC_shuffled_pp')

figure
plot(AUC_shuffled(i,:),'.')

% prctile(AUC_shuffled(:),95) = 0.5922  on 18/09/19 p_sh(1) = 0.5922;
% prctile(AUC_shuffled(:),99) = 0.6754  on 18/09/19 p_sh(2) = 0.6754;
p_sh = [prctile(AUC_shuffled(:),95) prctile(AUC_shuffled(:),99)];
% prctile(AUC_shuffled_pp(:),95) = 0.6057  on 18/09/19 pp_sh(1) = 0.6057
% prctile(AUC_shuffled_pp(:),99) = 0.7022  on 18/09/19 pp_sh(2) = 0.7022
pp_sh = [prctile(AUC_shuffled_pp(:),95) prctile(AUC_shuffled_pp(:),99)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% values for data from control experiments
% prctile(AUC_shuffled(:),95) = 0.6411  on 18/09/19 p_sh = [0.6411 0.7484 ]
% prctile(AUC_shuffled(:),99) =  0.7484  on 18/09/19
% prctile(AUC_shuffled_pp(:),95) = 0.6470  on 18/09/19 pp_sh = [0.6470 0.7414]
% prctile(AUC_shuffled_pp(:),99) = 0.7414  on 18/09/19
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% population reference
% select only units sign. modulated by vis. perturbation
BestUnits = [find(AM_UOI) find(AM_UOI)];
idx_1 = BestUnits(AUC_units(:,1)>p_sh(2),1);
idx_2 = BestUnits(AUC_units(:,2)>pp_sh(2),2);
RespUnits = false(length(AM_UOI),2);
RespUnits(idx_1,1) = true;
RespUnits(idx_2,2) = true;


%% self reference for each unit
for i = 1:length(AUC_units)
    AUC_p(i) = sum(AUC_units(i)>AUC_shuffled(i,:));
end

% show AUC for perturbation period
figure
plot(AUC_units(:,1),'r.')
hold on
plot(find(AUC_p>=size(AUC_shuffled,2)*0.99),AUC_units(find(AUC_p>=size(AUC_shuffled,2)*0.99),1),'ro')
plot([0 sum(SelectedCells)], [p_sh(1) p_sh(1)], 'k')
plot([0 sum(SelectedCells)], [p_sh(2) p_sh(2)], 'k--')
xlabel('units'); ylabel('AUC')
ll = legend({'AUC perturbation period','sign. AUC value (s.r.)','95% shuffled trial ID', '99% shuffled trial ID'})
ll = legend({'AUC perturbation period','95% shuffled trial ID', '99% shuffled trial ID'})
ll.Color = 'none'; ll.EdgeColor = 'none';
xlim([0 sum(SelectedCells)+1]); 
% show AUC for post perturbation period
figure
plot(AUC_units(:,2),'b.')
hold on
plot([0 sum(SelectedCells)], [pp_sh(1) pp_sh(1)], 'k')
plot([0 sum(SelectedCells)], [pp_sh(2) pp_sh(2)], 'k--')
xlabel('units'); ylabel('AUC')
legend({'AUC post-pert period','95% shuffle', '99% shuffle'})
xlim([0 sum(SelectedCells)+1]); 
% show AUC1 vs AUC2
figure
scatterDiagHist(AUC_units(:,1),AUC_units(:,2),[-1:0.02:1],'.k')
xlabel('AUC pert period')
ylabel('AUC post pert period')
%%
% show shuffling process
figure
h = histogram(AUC_shuffled(:),[0:0.02:1],'FaceColor',[0.5 0.5 0.5],'Normalization','probability','EdgeColor','none')
hold on
h1 = histogram(AUC_units(:,1),[0:0.02:1],'FaceColor',[1 0 0],'Normalization','probability','EdgeColor','none')
%plot([p_sh(1) p_sh(1)], [0 max(h.BinCounts/length(AUC_shuffled(:)))*1.1],'k') 
plot([p_sh(2) p_sh(2)], [0 max(h1.BinCounts/length(AUC_units(:,1)))*1.1],'k--','LineWidth',3)
xlabel('AUC of shuffled data'); ylabel('probability')
ll = legend({['AUC of shuffled data n=' num2str(sum(AM_UOI)) 'k'],['AUC of units n=' num2str(sum(AM_UOI))], '99% shuffled data'})
ll.Color = 'none'; ll.EdgeColor = 'none';
box off
set(gca,'FontSize',13,'TickDir','out');
%%

% save indexes of units responsive to visual perturbation
idx_temp_selection = (AUC_units(:,2)>=prctile(AUC_shuffled(:),95));
% get the units indexes in the main list
idx_temp_total= find(AM_UOI==1);
Units_idx_PertResp = idx_temp_total(idx_temp_selection);
AM_UOI(:,end+1) = zeros(length(AM_UOI),1);
AM_UOI(Units_idx_PertResp,end) = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on
plot([0 1],[0 1],'--k','LineWidth',1)

%% EXAMPLE UNITS
% plot units for which you showing the AUC

%% sort the cells by their response to perturbation (reference for heat map)
[B,I,PertTrial_Ex] = PlotSortingOrder(ProjectData,SelectedResponses,AM_Param,AM_Speed,Units_Sel,Diff);
% I = descending order taking into account:[pert response; 2s vis stim respo; post-pert reponse]
I = I(~isnan(B(:,1)));
B = reshape(B(~isnan(B(:))),length(B(~isnan(B(:))))/3,3);
%%
figure
param_sel = AM_Param(:,SelectedCells,:); UI = 407;
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
     
%% plot modulation indexes
param_sel = AM_Param(:,SelectedCells,:);
Param = squeeze(param_sel(:,:,3));
% Depth of modulation
r= 1; clear DM DepthMod
for i = 1:size(Param,2)
    % select trial responses
    Responses = Param(:,i);
    % scrub responses
    Responses = logical(Responses(~isnan(Responses)));
    UnitResponse = mean(squeeze(UnitResponses_smooth(Responses==1,i,:)))*60;

    Dir = sign(sum(UnitResponse(po_i(1):po_i(1)+trialSide_samples))-sum(UnitResponse(po_i(1)-trialSide_samples:po_i(1))));
    switch Dir
        case 1
            DM(i) = (max(UnitResponse(po_i(1):po_i(1)+trialSide_samples))-max(UnitResponse(po_i(1)-trialSide_samples:po_i(1))))/...
                (max(UnitResponse(po_i(1):po_i(1)+trialSide_samples))+max(UnitResponse(po_i(1)-trialSide_samples:po_i(1))));
        case -1
            DM(i) = (min(UnitResponse(po_i(1):po_i(1)+trialSide_samples))-min(UnitResponse(po_i(1)-trialSide_samples:po_i(1))))/...
                (min(UnitResponse(po_i(1):po_i(1)+trialSide_samples))+min(UnitResponse(po_i(1)-trialSide_samples:po_i(1))));
    end
    
    if AUC_units(i,1)>p_sh(2)
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
figure
histogram(DM_sh,[-1:0.05:1])
save('DM_pert_shuffled.mat','DM_sh')

% figure
% h = polarhistogram(DepthMod(:,1)*pi,[-pi:0.36:pi],'FaceColor',[1 0 0],'Normalization','probability','EdgeColor','none')
% %h.DisplayStyle = 'stairs';
% set(gca,'ThetaZeroLocation','top','ThetaDir','clockwise');
% set(gca,'RTickLabel',[])
% thetaticks(0:45:315); set(gca,'ThetaTickLabel',{'0'; '0.25'; '0.5'; '0.75'; '-1 | 1'; '-0.75'; '-0.5'; '-0.25'})
% hold on
% title('DM_p_e_r_t probability')
% %xlabel('Depth of Modulation'); ylabel('probability')
% title(['Sign. responsive units n=' num2str(size(DepthMod,1))])
% box off
% set(gca,'FontSize',13);
% text(-0.5,0.1,['n=' num2str(sum(DepthMod(:,1)<=0))],'FontSize',13)
% text(0.5,0.1,['n=' num2str(sum(DepthMod(:,1)>0))],'FontSize',13)

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

[h,p]=kstest2(DepthMod(DepthMod(:,1)>0,1),DM_h(DM_h>0))
[h,p]=kstest2(DepthMod(DepthMod(:,1)<=0,1),DM_h(DM_h<=0))


%% Plot Dm VS AUC
figure
plot(DM,AUC_units(:,1),'k.')
hold on
plot([-1 1], [pp_sh(2) pp_sh(2)], 'k:')
hold on
plot([prctile(DM_sh(:),1) prctile(DM_sh(:),1)], [0 1], 'k:')
plot([prctile(DM_sh(:),99) prctile(DM_sh(:),99)], [0 1], 'k:')
xlabel('Depth of modulation'); ylabel('AUC');
box off
grid on
set(gca,'FontSize',13,'FontName','Arial','TickDir','out');
xlim([-1 1]); ylim([0.4 1])


%% entire population (1493 units --> 830 units) reference
% show me relationship to vis.stim. pert. onset and pert. offset
idx_1 = find(AUC_units(:,1)>p_sh(2));
idx_2 = AUC_units(:,2)>pp_sh(2);
RespUnits = false(sum(AM_UOI),2); 
RespUnits(idx_1,1) = true; % responsive to perturbation onset
RespUnits(idx_2,2) = true; % responsive to perturbation offset
% recordings of responsive units
unique(squeeze(param_sel(1,find(AUC_units(:,1)>p_sh(2)),1)))

load('DM_visresp_shuffled.mat')
% The 5 metrics for measuring visual responses are
% 1) mean FR of first 4 seconds of stim
% 2) mean FR of 1s before stim onset
% 3) DM
% 4) MI
% 5) avg baseline FR
DM_vis = squeeze(FR_UOI_VIS(3,:,SelectedCells));
DM_vis_sh = [prctile(DM_vis(2:end,:),1); prctile(DM_vis(2:end,:),99)];
RespUnits(:,3) = DM_vis(1,:)<=DM_vis_sh(1,:) | DM_vis(1,:)>=DM_vis_sh(2,:);

RespUnits_Sel = RespUnits; 
SumUnits(1) = sum(RespUnits_Sel(:,1) & RespUnits_Sel(:,2));
SumUnits(2) = sum(RespUnits_Sel(:,1) & RespUnits_Sel(:,3));
SumUnits(3) = sum(RespUnits_Sel(:,1) & RespUnits_Sel(:,2) & RespUnits_Sel(:,3));

% show MI pert VS MI of vis.stim. 
figure
%plot((DM_vis(1,:)), (DM) ,'.')
set(0, 'DefaultFigureRenderer', 'opengl');
hold on
plot((DM(RespUnits(:,3))),(DM_vis(1,RespUnits(:,3))) ,'b.') % vis. stim. responsive
grid on
set(gca,'box','off','TickDir','out');
xlim([-1 1]); ylim([-1 1]);
xlabel('Perturbation MI'); ylabel('Visual stim. MI');

figure
%plot((DM_vis(1,:)), (DM) ,'.')
set(0, 'DefaultFigureRenderer', 'opengl');
hold on
plot((DepthMod(:,1)), (DM_vis(1,DepthMod(:,3))),'ro') % perturbation responsive
grid on
set(gca,'box','off','TickDir','out')
xlim([-1 1]); ylim([-1 1]);
xlabel('Perturbation MI'); ylabel('Visual stim. MI');

figure
plot((DM_vis(1,:)),'.')

figure
subplot(2,2,1)
plot((DepthMod(:,1)), (DM_vis(1,DepthMod(:,3))),'ro') % perturbation responsive
hold on
plot((DM(RespUnits(:,3))),(DM_vis(1,RespUnits(:,3))) ,'b.') % vis. stim. responsive
grid on
set(gca ,'box','off','TickDir','out','YAxisLocation','left','XAxisLocation','top');
xlabel('Perturbation MI'); ylabel('Visual Stim. MI');
xlim([-1 1]); ylim([-1 1]);

subplot(2,2,3)
histogram((DM(RespUnits(:,3))),-1:0.05:1,'Normalization','probability','FaceColor','b','EdgeColor',[1 0 0])
xlabel('Perturbation MI'); ylabel('Probability');
set(gca,'box','off','TickDir','out'); grid on; xlim([-1 1]);
set(gca,'view',[0 -90]);

subplot(2,2,2)
histogram((DM_vis(1,RespUnits(:,3))),-1:0.05:1,'Normalization','probability','FaceColor','b','EdgeColor',[1 0 0])
xlabel('Visual Stim. MI'); ylabel('Probability');
set(gca,'box','off','TickDir','out','YAxisLocation','right','XAxisLocation','top'); grid on; xlim([-1 1]);
set(gca,'view',[90 -90]);


%%
figure
hold on
grid on
for i = 1:length(DM)
    if AUC_units(i)>p_sh(2)
        plot(i,DM(i),'.r')
        plot(i,DM(i),'or')
    else
        plot(i,DM(i),'ok')
    end
end
% Modulation index
figure
hold on
for i = 1:size(DiffResponsePercent_Pert,2)
    UnitResponse = mean(squeeze(UnitResponses_smooth(param_sel(:,i,3)==1,i,:)))*60;
    MI = (max(UnitResponse(po_i(1):po_i(1)+trialSide_samples)))/...
        (max(UnitResponse(po_i(1):po_i(1)+trialSide_samples))+max(UnitResponse(po_i(1)-trialSide_samples:po_i(1))));
    if AUC_units(i)>p_sh
        plot(i,nanmean(DiffResponsePercent_Pert(Param(:,i)==1,i,6)),'.r')
        plot(i,nanmean(DiffResponsePercent_Pert(Param(:,i)==1,i,6)),'or')
    else
        plot(i,nanmean(DiffResponsePercent_Pert(Param(:,i)==1,i,6)),'ok')
    end
end


