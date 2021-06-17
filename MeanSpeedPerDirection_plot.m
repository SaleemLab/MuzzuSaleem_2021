%% Tomaso Muzzu - UCL - 16/01/2019

% function to plot speed information of a single or multiple sessions


function MeanSpeedPerDirection_plot(Array,ES,StimulusParams,FileNames, RunningThreshold, MaxSpeed)

% Reminder of how Array is structured
% Array{i,1} = i; % order of trials of display
% Array{i,2} = ES.OrienTrial(i); % grating direction
% Array{i,3} = ind_1(i); % order of trials by ranking of choice
% Array{i,4} = SpeedTrial(i,:); % speed of every trial
% Array{i,5} = ES.PertTrial(i); % perturbation
% Array{i,6} = ind_2(cc);
% Array{i,7} = SpeedTrial_PertStartEnds(cc,:);

if ~exist('MaxSpeed','var')
    MaxSpeed = 50;
end

% find number of sessions store in Array
sessions = 0; clear sessions_i
for i=1:size(Array,1)
    if Array{i,1}==1
        sessions = sessions + 1;
        sessions_i(sessions) = i;
    end
end

clear SpeedTrial DirTrial RankTemp Rank Rank_i
SpeedTrial = cell2mat(Array.Speed_in_trial);
DirTrial = Array.Direction;
RankTemp = Array.Avg_Speed4Rank;
[Rank Rank_i] = sort(RankTemp);
clear SpeedTrial_sorted
SpeedTrial_sorted = SpeedTrial(Rank_i,:);
DirTrial_sorted = DirTrial(Rank_i);
Directions = unique(DirTrial);
MeanSpeedTrial = mean(SpeedTrial_sorted(:,60:end-60),2);
MeanSpeedTrial(MeanSpeedTrial>50) = NaN;
clear trials2plot trials2plot_groups
trials2plot = 0; trials2plot_groups = 1;
for i = 1:length(Directions)
    trials2plotTmp = find(MeanSpeedTrial>RunningThreshold & MeanSpeedTrial<MaxSpeed & DirTrial_sorted==Directions(i));
    trials2plot = [trials2plot; trials2plotTmp];
    trials2plot_groupsTmp = i*ones((length(trials2plotTmp)),1);
    trials2plot_groups = [trials2plot_groups; trials2plot_groupsTmp];
    n_info(i) = length(trials2plotTmp);
end
trials2plot = trials2plot(2:end);  trials2plot_groups = trials2plot_groups(2:end);

figure
set(gcf,'Renderer', 'painters', 'Position', [-1000 0 1050 500])
hold off
% produce a plot that shows the means and SE of the speed for the different
notBoxPlot(MeanSpeedTrial(trials2plot),trials2plot_groups,...
            'style','patch','jitter',0.4,'interval','tInterval','markMedian',true)
hold on
for k1 = 1:size(n_info,2)
    text(k1-0.15,2, sprintf('N = %d', n_info(k1)), 'FontSize',10);
end
set(gca,'TickDir','out','xTickLabels',Directions,'box', 'off')
ylim([0 max(MeanSpeedTrial(trials2plot))+5]);
xlabel('Grating direction (^o)'); ylabel('Mean speed (cm/s)')

a = axes;
MouseName_init = strfind(FileNames{1,1},'\M');
t1 = title(['Speed profiles during ' num2str(length(trials2plot)) ' trials of ' num2str(length(FileNames)) ' sessions, '  ...
             num2str(RunningThreshold) '<speed <' num2str(MaxSpeed) 'cm/s, '...
            FileNames{1,1}(MouseName_init+1:MouseName_init+1+6) 'TM' FileNames{1,1}(MouseName_init+1+11:MouseName_init+1+12)]);
a.Visible = 'off'; % set(a,'Visible','off');
t1.Visible = 'on'; % set(t1,'Visible','on');


%EOF
end