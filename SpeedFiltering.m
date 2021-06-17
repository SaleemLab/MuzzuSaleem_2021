%% Tomaso Muzzu - UCL - 16/09/2019

% select trials based on running speed measures


function [SelectedTrials_OI_final SelectedTrials_OI_control_final TotTrials] = SpeedFilteringONOFF(AM_Param, AM_SpeedResponse, AM_SpeedResponseControl, Trials_PertON,Trials_PertOFF, Speed_TH, minTrials)

if ~exist('minTrials','var')
    minTrials = 5; % minimum number of trials per condition per recording
end

% AM_SpeedResponse(:,:,1) = mean speed -0.5s to +2s at visual stimulus onset
% AM_SpeedResponse(:,:,2) = mean speed -0.5s to end of perturbation period
% AM_SpeedResponse(:,:,3) = mean speed -0.5s to +1s at perturbation offset
% AM_SpeedResponse(:,:,4) = mean speed across entire trial
clear SelectedTrials_temp  SelectedTrials_OI  SelectedTrials_OI_control
for i = 1:size(AM_SpeedResponse,3)
    % select trials with perturbation and running
    SelectedTrials_temp = (AM_SpeedResponse(:,:,i)>Speed_TH); % pert ON
    %SelectedTrials_temp(SelectedTrials_temp(:)==0) = nan;
    SelectedTrials_OI{i,1} = SelectedTrials_temp;
    % select trials with perturbation and stationary
    SelectedTrials_temp = (AM_SpeedResponse(:,:,i)<Speed_TH); % pert ON
    %SelectedTrials_temp(SelectedTrials_temp(:)==0) = nan;
    SelectedTrials_OI{i,2} = SelectedTrials_temp;
    % select trials without perturbation and running
    SelectedTrials_temp = (AM_SpeedResponseControl(:,:,i)>Speed_TH); % pert OFF
    %SelectedTrials_temp(SelectedTrials_temp(:)==0) = nan;
    SelectedTrials_OI_control{i,1} = SelectedTrials_temp;
    % select trials without perturbation and stationary
    SelectedTrials_temp = (AM_SpeedResponseControl(:,:,i)<Speed_TH); % pert OFF
    %SelectedTrials_temp(SelectedTrials_temp(:)==0) = nan;
    SelectedTrials_OI_control{i,2} = SelectedTrials_temp;
end
    
%% select only the recordings in which there are at least 'minTrials' running and 10 stationary trials
% count the number of trials for each condition
clear TotTrials
for i = 1 : size(AM_SpeedResponse,3)
    TotTrials{i,1} = [nansum(SelectedTrials_OI{i,1},1) ; nansum(SelectedTrials_OI{i,2},1)]; % pert ON
    TotTrials{i,2} = [nansum(SelectedTrials_OI_control{i,1},1) ; nansum(SelectedTrials_OI_control{i,2},1)];  % pert OFF
end
clear SelectedUnits_OI SelectedUnits_OI_control
for i = 1 : size(AM_SpeedResponse,3)
    SelectedUnits_OI(i,:) = TotTrials{i,1}(1,:)>minTrials & TotTrials{i,1}(2,:)>minTrials; % sum of the pert ON run and still
    SelectedUnits_OI_control(i,:) = TotTrials{i,2}(1,:)>minTrials & TotTrials{i,2}(2,:)>minTrials; % sum of the pert OFF run and still
end
clear SelectedTrials_OI_final SelectedTrials_OI_control_final
for i = 1 : size(AM_SpeedResponse,3) % 1=VS onset; 2=pertOnset; 3=pertOffset; 4=alltrialAVG.
    % pertON & run trials of recs with enough trials for both conditions
    SelectedTrials_OI_final{i,1} = SelectedTrials_OI{i,1} & SelectedUnits_OI(i,:); 
    % pertON & stationary trials of recs with enogh trials for both conditions
    SelectedTrials_OI_final{i,2} = SelectedTrials_OI{i,2} & SelectedUnits_OI(i,:); 
    % pertOFF & run trials of recs with enogh trials for both conditions
    SelectedTrials_OI_control_final{i,1} = SelectedTrials_OI_control{i,1} & SelectedUnits_OI_control(i,:); 
    % pertOFF & stationary trials of recs with enogh trials for both conditions
    SelectedTrials_OI_control_final{i,2} = SelectedTrials_OI_control{i,2} & SelectedUnits_OI_control(i,:);      
end


end



