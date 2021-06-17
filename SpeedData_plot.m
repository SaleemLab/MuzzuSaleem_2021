%% Tomaso Muzzu - UCL - 16/01/2019

% function to plot speed information of a single or multiple sessions


function SpeedData_plot(Array,ES,StimulusParams,FileNames, RunningThreshold, MaxSpeed)

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

% show contrast and TF dynamics
% show different grating directions
% take multiple sessions and look at the different

figure
set(gcf,'Renderer', 'painters', 'Position', [-1000 0 1050 900])
subplot(5,1,1)
if isfield(ES,'PertStartsEnds')
    Pert_Start = min(find(StimulusParams(:,3)==0));
    Pert_Stop = max(find(StimulusParams(:,3)==0));
    yyaxis left
    plot([StimulusParams(Pert_Start,1) StimulusParams(Pert_Stop,1)], [max(StimulusParams(:,3)) max(StimulusParams(:,3))])
    hold on
end
yyaxis left
plot(StimulusParams(1:min(find(StimulusParams(2:end,1)==0)),1),StimulusParams(1:min(find(StimulusParams(2:end,1)==0)),3))
xlim([min(StimulusParams(1:min(find(StimulusParams(2:end,1)==0)),1)) max(StimulusParams(1:min(find(StimulusParams(2:end,1)==0)),1))]); 
ylim([0 3.2]);
ylabel('TF'); 
yyaxis right
plot(StimulusParams(1:min(find(StimulusParams(2:end,1)==0)),1),StimulusParams(1:min(find(StimulusParams(2:end,1)==0)),2))
xlim([min(StimulusParams(1:min(find(StimulusParams(2:end,1)==0)),1)) max(StimulusParams(1:min(find(StimulusParams(2:end,1)==0)),1))]); 
ylim([0 1]);
ylabel('contrast'); set(gca,'TickDir','out');
box off
%legend({'TF','Contrast'})
subplot(5,1,2:5)
hold off
%$set(gca,'zscale','log')
colormap(parula(256))
MeanSpeedTrial = mean(SpeedTrial_sorted(:,60:end-60),2);
trials2plot = find(MeanSpeedTrial>RunningThreshold & MeanSpeedTrial<MaxSpeed);
imagesc(StimulusParams(1:min(find(StimulusParams(2:end,1)==0)),1),linspace(1,length(trials2plot),length(trials2plot)), ...
    SpeedTrial_sorted(trials2plot,:));
colorbar
xlabel('seconds'); ylabel('trials');

a = axes;
MouseName_init = strfind(FileNames{1,1},'\M');
t1 = title(['Speed profiles during ' num2str(length(trials2plot)) ' trials of ' num2str(length(FileNames)) ' sessions, '  ...
             num2str(RunningThreshold) 'cm/s < speed <' num2str(MaxSpeed) 'cm/s, '...
            FileNames{1,1}(MouseName_init+1:MouseName_init+1+6) 'TM' FileNames{1,1}(MouseName_init+1+11:MouseName_init+1+12)]);
a.Visible = 'off'; % set(a,'Visible','off');
t1.Visible = 'on'; % set(t1,'Visible','on');




%EOF
end