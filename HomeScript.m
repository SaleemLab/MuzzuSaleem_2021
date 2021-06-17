%% Tomaso Muzzu - UCL - 08 Jan 2018

%% homebase script
clearvars -except StatsM
close all

FileNames = SelectSessions;
        
for i = 1:length(FileNames)
    
    ES_ = getES(FileNames(i));
    [Array_, StimulusParams_] = ...
        GetSpeedData(ES_,'Direction',[0:45:315],...
        'Running_Th',3, ...
        'Pert_Status',2, ... % if 1 -> trials with perturbation; if 0 -> trials without perturbation
        'Rank_mode','AVGspeed',... % 'AVGspeed','TrialStart' speed at trial onset or 'PertStart' for speed at pert onset
        'Norm_mode','none', ... % 'MaxSpeedTrial', 'MaxSpeedSession', 'z_score', 'none'
        'log_speed',0, ... % 1 if computing speed in log scale
        'AvgSpeedPeriod', [0 +1], ... % max -1 and +1 from events specified at Rank_mode
        'Smoothing',0.6); % in seconds, size of Gaussian filter
    if i == 1
        Array = Array_;
        StimulusParams = StimulusParams_;
        ES = ES_;
    else
        Array = [Array; Array_];
        StimulusParams = [StimulusParams; StimulusParams_];
        ES = [ES; ES_];
    end
    clear Array_ StimulusParams_ ES_
end

%% BEHAVIOUR
% show speed profiles with different running thresholds values
for rth = [1 5 10 15]
    MaxSpeed = 50;
    SpeedData_plot(Array,ES,StimulusParams,FileNames,rth,MaxSpeed);
end

% show mean speed for different grating directions
for rth = [1 5 10 15]
    MaxSpeed = 50;
    MeanSpeedPerDirection_plot(Array,ES,StimulusParams,FileNames,rth,MaxSpeed);
end

% show mean speed for different grating directions
rth = 3; MaxSpeed = 50;
if exist('StatsM','var')
    StatsM(end+1,:) = GeneralSessionStats_plot(Array,ES,StimulusParams,FileNames,rth,MaxSpeed);
else
    StatsM = GeneralSessionStats_plot(Array,ES,StimulusParams,FileNames,rth,MaxSpeed);
end

MultiMouseStats_Behaviour(StatsM,rth);

%% NEURAL RESPONSES
% show cell responses to perturbation
rth = 3; MaxSpeed = 50;
Direction = unique(Array.Direction);
for i = 1:length(Direction)
    PSTH_AllCells_pert_angle(Array,ES.SpikesTrial,ES,FileNames{1,1},ES.RecordingOI,Direction(i), rth, MaxSpeed,StimulusParams)
end

rth = 3; MaxSpeed = 50;
Direction = unique(Array.Direction);
for i = 1:length(Direction)
    PSTH_AllCells_pert_angle_moving(Array,ES.SpikesTrial,ES,FileNames{1,1},ES.RecordingOI,Direction(i),rth,MaxSpeed,StimulusParams)
end

% cell by cell responses to perturbation
rth = 3; MaxSpeed = 50;
Direction = unique(Array.Direction);
SingleCell_response_Perturbation(Array,ES.SpikesTrial,ES,FileNames{1,1},ES.RecordingOI,Direction(1), rth, MaxSpeed,StimulusParams)

rth = 3; MaxSpeed = 50;
Direction = unique(Array.Direction);
value2show = 0; % 0 for difference between pert and no pert, 1 for actual response magnitude
SingleCell_response_Perturbation_Color(Array,ES.SpikesTrial,ES,FileNames{1,1},ES.RecordingOI,Direction(1), rth, MaxSpeed,StimulusParams,value2show)
value2show = 1; % 0 for difference between pert and no pert, 1 for actual response magnitude
SingleCell_response_Perturbation_Color(Array,ES.SpikesTrial,ES,FileNames{1,1},ES.RecordingOI,Direction(1), rth, MaxSpeed,StimulusParams,value2show)
                                             
% polar plot for tuning properties of every cell to grating direction and
% to perturbation direction
rth = 3; MaxSpeed = 50;
value2show = 0; % 0 for difference between pert and no pert, 1 for actual response magnitude
Direction = unique(Array.Direction);
SingleCell_response_PolarPlots(Array,ES.SpikesTrial,ES,FileNames{1,1},ES.RecordingOI,Direction, rth, MaxSpeed,StimulusParams,value2show)


% tuning properties of the neurons
% speed

% grating direction


% figure
% plot(ES.Time,ES.Contrast/max(ES.Contrast));
% hold on
% plot(ES.Time,ES.TF/max(ES.TF))
% hold on
% plot(ES.Time(ES.trialStartsEnds(1:2:end)+1),ones(length(ES.trialStartsEnds(1:2:end)),1)*0.4,'or')
% hold on
% plot(ES.Time(ES.trialStartsEnds(:,1)),ones(length(ES.trialStartsEnds(:,1)),1)*0.3,'or') % trial onsets
% hold on
% plot(ES.Time(ES.trialStartsEnds(:,2)),ones(length(ES.trialStartsEnds(:,2)),1)*0.3,'ob') % trial offsets
% 
