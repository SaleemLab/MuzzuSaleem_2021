%% Tomaso Muzzu - UCL - 13/09/2019

% plot basic responses of units for visual perturbation experiments

function PlotResponsesColormap(ProjectData,SelectedResponses,SelectedResponses1,AM_Param,AM_Speed,Units_Sel,Diff)


%% sort the cells by their response to perturbation 
[B,I,PertTrial_Ex] = PlotSortingOrder(ProjectData,SelectedResponses1,AM_Param,AM_Speed,Units_Sel,Diff);
% I = descending order taking into account:[pert response; 2s vis stim respo; post-pert reponse]
I = I(~isnan(B(:,1)));
B = reshape(B(~isnan(B(:))),length(B(~isnan(B(:))))/3,3);

% example perturbation trial
PertTrial_Ex(1) = 1;
PertTrial_Ex(2) = min(find(ProjectData.Trial_data{1,1}.PerturbationON==1))
% PertTrial_Ex = [recording nr, trial nr];


%% 
%figure
trialSide_samples = 60;
trialSide_seconds = 1;
TimeLine = linspace(min(ProjectData.Trial_data{PertTrial_Ex(1),1}.time_in_trial{PertTrial_Ex(2)}), ...
                    max(ProjectData.Trial_data{PertTrial_Ex(1),1}.time_in_trial{PertTrial_Ex(2)}),...
                    size(SelectedResponses{1},2)) -trialSide_seconds; % -1 seconds
set(gcf, 'PaperUnits', 'centimeters ');
%set(gcf, 'Position', [0 0 20 30]);
set(gcf, 'Renderer', 'painters');


if Diff==0
    imagesc(TimeLine,1:size(SelectedResponses{1}(I,:),1),SelectedResponses{1}(I,:));
    %control
    %imagesc(TimeLine,1:size(SelectedResponses{2}(I,:),1),SelectedResponses{2}(I,:));
else
    imagesc(TimeLine,1:size(SelectedResponses{1}(I,:),1),(abs((SelectedResponses{1}(I,:))-(SelectedResponses{2}(I,:)))));
end

set(gca,'fontsize',13,'TickDir','out');
xlabel('seconds','fontsize',13); ylabel('units','fontsize',13);
xlim([-0.3 max(TimeLine)-0.3]);
box off



end