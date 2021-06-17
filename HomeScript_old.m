%% Tomaso Muzzu - UCL - 08 Jan 2018

%% homebase script

ES = getES;

GetSpeedData

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   PLOTTING FUNCTIONS
% 4) scan through the cells and plot basic PSTH
for i = 1:size(SpikesTrial,1)   
    PSTH_basic(i,SpikesTrial,ES,FileName,RecordingOI);
end

% 5) scan through the cells and plot basic PSTH distinguishing trials with perturbation
for i = 1:size(SpikesTrial,1)   
    PSTH_basic_pert(i,SpikesTrial,ES,FileName,RecordingOI);
end

% 6) plot all the cells and plot basic PSTH distinguishing trials with
% perturbation and select only a specific direction of stimulus
for i = 1:size(SpikesTrial,1)  
    for Direction = unique(ES.OrienTrial)
        PSTH_basic_pert_angle(i,SpikesTrial,ES,FileName,RecordingOI,Direction);
    end
end

% 7) plot all the cells and plot basic PSTH distinguishing trials with
% perturbation and select only a specific direction of stimulus
for Direction = unique(ES.OrienTrial)
    PSTH_AllCells_pert_angle(SpikesTrial,ES,FileName,RecordingOI,Direction);
end

% 8) plot all the cells and plot basic PSTH distinguishing trials with
% perturbation and select only a specific direction of stimulus0
for Direction = unique(ES.OrienTrial)
    PSTH_AllCells_pert_angle_moving(SpikesTrial,ES,FileName,RecordingOI,Direction);
end


% 9) plot the speed signal and how it changes depending on the time of the
% trial.
for Direction = unique(ES.OrienTrial)
    SpeedPlot_pert_angle(SpikesTrial,ES,FileName,RecordingOI,Direction) %Direction, Running threshold, Perturbation?
end


figure
plot(ES.Time,ES.Contrast/max(ES.Contrast));
hold on
plot(ES.Time,ES.TF/max(ES.TF))
hold on
plot(ES.Time(ES.trialStartsEnds(1:2:end)+1),ones(length(ES.trialStartsEnds(1:2:end)),1)*0.4,'or')
hold on
plot(ES.Time(ES.trialStartsEnds(:,1)),ones(length(ES.trialStartsEnds(:,1)),1)*0.3,'or') % trial onsets
hold on
plot(ES.Time(ES.trialStartsEnds(:,2)),ones(length(ES.trialStartsEnds(:,2)),1)*0.3,'ob') % trial offsets

