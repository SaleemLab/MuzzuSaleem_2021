%% Tomaso Muzzu - UCL - 13/09/2019

% plot basic responses of units for visual perturbation experiments

function PlotResponses(ProjectData,SelectedResponses,AM_Param,AM_Speed,Units_Sel,Diff)


%% sort the cells by their response to perturbation 
[B,I,PertTrial_Ex] = PlotSortingOrder(ProjectData,SelectedResponses,AM_Param,AM_Speed,Units_Sel,Diff);
% I = descending order taking into account:[pert response; 2s vis stim respo; post-pert reponse]
I = I(~isnan(B(:,1)));
B = reshape(B(~isnan(B(:))),length(B(~isnan(B(:))))/3,3);

% example perturbation trial
PertTrial_Ex(1) = 1;
PertTrial_Ex(2) = min(find(ProjectData.Trial_data{1,1}.PerturbationON==1))
% PertTrial_Ex = [recording nr, trial nr];


%% 
figure
trialSide_samples = 60;
trialSide_seconds = 1;
TimeLine = linspace(min(ProjectData.Trial_data{PertTrial_Ex(1),1}.time_in_trial{PertTrial_Ex(2)}), ...
                    max(ProjectData.Trial_data{PertTrial_Ex(1),1}.time_in_trial{PertTrial_Ex(2)}),...
                    size(SelectedResponses{1},2)) -trialSide_seconds; % -1 seconds
set(gcf, 'PaperUnits', 'centimeters ');
%set(gcf, 'Position', [0 0 20 30]);
set(gcf, 'Renderer', 'painters');
subplot(10,1,[5:10])
colormap(jet(256))
if Diff==0
    imagesc(TimeLine,1:size(SelectedResponses{1}(I,:),1),SelectedResponses{1}(I,:));
    %control
    %imagesc(TimeLine,1:size(SelectedResponses{2}(I,:),1),SelectedResponses{2}(I,:));
else
    imagesc(TimeLine,1:size(SelectedResponses{1}(I,:),1),(abs((SelectedResponses{1}(I,:))-(SelectedResponses{2}(I,:)))));
end
colorbar
set(gca,'fontsize',13,'TickDir','out');
xlabel('seconds','fontsize',13); ylabel('units','fontsize',13);
xlim([-0.3 max(TimeLine)-0.3]);
box off

% plot also contrast and TF
subplot(10,1,1)
Contrast = ProjectData.Trial_data{PertTrial_Ex(1),1}.contrast_in_trial{PertTrial_Ex(2),1}(1:end-1);
TF = ProjectData.Trial_data{PertTrial_Ex(1),1}.TF_in_trial{PertTrial_Ex(2),1}(1:end-1);
%TF = [max(TF)*(ones(1,po_i(1))) zeros(1,po_i(2)-po_i(1)) max(TF)*(ones(1,length(TF)-po_i(2)))];
h1 = plot(TimeLine,Contrast,'k','LineWidth',2)
hold on
h2 = plot(TimeLine,TF,'b','LineWidth',2)
hold on
plot([0 0],[0 max(TF)],'--r','LineWidth',0.75);
hold on
TrialEnd = TimeLine(min(find(Contrast(length(Contrast)/2:end)==0))+round(length(Contrast)/2));
plot([TrialEnd TrialEnd ], [0 max(TF)],'--r','LineWidth',0.75);
hold on
plot([TimeLine(min(find(TF==0))), TimeLine(max(find(TF==0)))], [max(TF) max(TF)],'--b','LineWidth',1);
xlim([-0.3 max(TimeLine)-0.3]);
ylim([0 max(TF)+0.4]);
box off
set(gca,'TickDir','out','XTickLabel',[],'YTick',[0 max(Contrast) max(TF)],'YTickLabel',[0 max(Contrast) max(TF)])    
ll = legend({'Contrast', 'TF'},'location','northwest','fontsize',10)
ll.Color = 'none'; ll.EdgeColor = 'none'; 
set(gca,'FontSize',13);
colorbar
title(['Response during perturbation trials'])


% plot avg responses to normal and perturbation trials
subplot(10,1,[2:4])
colormap(jet(256))
if Diff == 0
    h1 = shadedErrorBar(TimeLine(1:size(SelectedResponses{1},2)),...
        smooth(mean(SelectedResponses{1}(~isnan(sum(SelectedResponses{1},2)),:),1),5),...
        std(SelectedResponses{1}(~isnan(sum(SelectedResponses{1},2)),:),1)/sqrt(length(SelectedResponses{1}(~isnan(sum(SelectedResponses{1},2)),:))), ...
        'lineprops',{'r-','markerfacecolor','r'});
    % mean(NormFR_allcells_sem(~isnan(sum(NormFR_allcells,2)),:),1), ...
    %
    hold on
    h2 = shadedErrorBar(TimeLine(1:size(SelectedResponses{2},2)),...
        smooth(mean(SelectedResponses{2}(~isnan(sum(SelectedResponses{2},2)),:),1),5),...
        std(SelectedResponses{2}(~isnan(sum(SelectedResponses{2},2)),:),1)/sqrt(length(SelectedResponses{2}(~isnan(sum(SelectedResponses{2},2)),:))), ...
        'lineprops',{'k-','markerfacecolor','k'});
    % mean(NormFR_allcells_off_sem(~isnan(sum(NormFR_allcells_off,2)),:) ,1))
else
    Difference = (abs(mean(SelectedResponses{1}(~isnan(sum(SelectedResponses{1},2)),:),1)-mean(SelectedResponses{2}(~isnan(sum(SelectedResponses{2},2)),:),1)));
    DifferenceSTD = std(abs((SelectedResponses{1}(~isnan(sum(SelectedResponses{1},2)),:))-(SelectedResponses{2}(~isnan(sum(SelectedResponses{2},2)),:))));
    h1 = shadedErrorBar(TimeLine(1:size(SelectedResponses{2},2)),...
        smooth(Difference,5),...
        DifferenceSTD/sqrt(size(SelectedResponses{1},1)), ...
        'lineprops',{'k-','markerfacecolor','r'});
end
hold on
plot([0 0],[0 1],'--r','LineWidth',0.75);
hold on
plot([TrialEnd TrialEnd ], [0 1],'--r','LineWidth',0.75);
hold on
plot([TimeLine(min(find(TF==0))), TimeLine(min(find(TF==0)))], [0 1],'--b','LineWidth',1);
hold on
plot([TimeLine(max(find(TF==0))), TimeLine(max(find(TF==0)))], [0 1],'--b','LineWidth',1);
colorbar
set(gca,'FontSize',13);
ylabel('norm. response','fontsize',13);
box off
set(gca,'TickDir','out','XTickLabel',[]);
TotNrPertTrials = 0;
xlim([-0.3 max(TimeLine)-0.3]);
YtraceMax = max(max(mean(SelectedResponses{2}(~isnan(sum(SelectedResponses{2},2)),:),1)), max(mean(SelectedResponses{1}(~isnan(sum(SelectedResponses{1},2)),:),1)));
YtraceMin = min(min(mean(SelectedResponses{2}(~isnan(sum(SelectedResponses{2},2)),:),1)), min(mean(SelectedResponses{1}(~isnan(sum(SelectedResponses{1},2)),:),1)));
%
if Diff == 0
    ll = legend({'perturbation ON', 'perturbation OFF'},'fontsize',10);
    ll.Color = 'none'; ll.EdgeColor = 'none';
  %  ylim([YtraceMax-0.3 YtraceMax+0.1])
else
    ll = legend('abs(pertON-pertOFF)','fontsize',10);
    ll.Color = 'none'; ll.EdgeColor = 'none';
   % ylim([0 max(Difference)+0.1])
end


end