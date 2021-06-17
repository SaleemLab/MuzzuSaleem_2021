%% Tomaso Muzzu - UCL - 13/11/2019

% plot mean responses curves for different grating directions

function PlotMeanCurveAngle(ProjectData,SelectedResponsesSEL,AM_Param,AM_Speed,Units_Sel,AngleCombo)

figure
set(gcf, 'PaperUnits', 'centimeters ');
%set(gcf, 'Position', [0 0 20 30]);
set(gcf, 'Renderer', 'painters');
trialSide_samples = 60;
trialSide_seconds = 1;
PertTrial_Ex(1) = 1;
PertTrial_Ex(2) = min(find(ProjectData.Trial_data{1,1}.PerturbationON==1))
TimeLine = linspace(min(ProjectData.Trial_data{PertTrial_Ex(1),1}.time_in_trial{PertTrial_Ex(2)}), ...
                    max(ProjectData.Trial_data{PertTrial_Ex(1),1}.time_in_trial{PertTrial_Ex(2)}),...
                    size(SelectedResponsesSEL{1},2)) -trialSide_seconds; % -1 seconds
Contrast = ProjectData.Trial_data{PertTrial_Ex(1),1}.contrast_in_trial{PertTrial_Ex(2),1}(1:end-1);
TF = ProjectData.Trial_data{PertTrial_Ex(1),1}.TF_in_trial{PertTrial_Ex(2),1}(1:end-1);
TrialEnd = TimeLine(min(find(Contrast(length(Contrast)/2:end)==0))+round(length(Contrast)/2));

h1 = shadedErrorBar(TimeLine(1:size(SelectedResponsesSEL{1},2)),...
        smooth(mean(SelectedResponsesSEL{1}(~isnan(sum(SelectedResponsesSEL{1},2)),:),1),5),...
        std(SelectedResponsesSEL{1}(~isnan(sum(SelectedResponsesSEL{1},2)),:),1)/sqrt(size(SelectedResponsesSEL{1}(~isnan(sum(SelectedResponsesSEL{1},2)),:),1)), ...
        'lineprops',{'r-','markerfacecolor','r'});
   
hold on
h2 = shadedErrorBar(TimeLine(1:size(SelectedResponsesSEL{2},2)),...
        smooth(mean(SelectedResponsesSEL{2}(~isnan(sum(SelectedResponsesSEL{2},2)),:),1),5),...
        std(SelectedResponsesSEL{2}(~isnan(sum(SelectedResponsesSEL{2},2)),:),1)/sqrt(size(SelectedResponsesSEL{2}(~isnan(sum(SelectedResponsesSEL{2},2)),:),1)), ...
        'lineprops',{'k-','markerfacecolor','k'});
hold on
plot([0 0],[0 1],'--k','LineWidth',0.75);
hold on
plot([TrialEnd TrialEnd ], [0 1],'--k','LineWidth',0.75);
hold on
plot([TimeLine(min(find(TF==0))), TimeLine(min(find(TF==0)))], [0 1],'--r','LineWidth',1);
hold on
plot([TimeLine(max(find(TF==0))), TimeLine(max(find(TF==0)))], [0 1],'--r','LineWidth',1);

set(gca,'FontSize',13);
ylabel('norm. response','fontsize',13);
xlabel('seconds','fontsize',13);
box off
set(gca,'TickDir','out');
TotNrPertTrials = 0;
xlim([-0.8 max(TimeLine)-0.3]);
YtraceMax = max(max(mean(SelectedResponsesSEL{2}(~isnan(sum(SelectedResponsesSEL{2},2)),:),1)), max(mean(SelectedResponsesSEL{1}(~isnan(sum(SelectedResponsesSEL{1},2)),:),1)));
YtraceMin = min(min(mean(SelectedResponsesSEL{2}(~isnan(sum(SelectedResponsesSEL{2},2)),:),1)), min(mean(SelectedResponsesSEL{1}(~isnan(sum(SelectedResponsesSEL{1},2)),:),1)));

ll = legend({['Direction ' num2str(AngleCombo(1,:))], ['Direction ' num2str(AngleCombo(2,:))]},'fontsize',10);
ll.Color = 'none'; ll.EdgeColor = 'none';

end
