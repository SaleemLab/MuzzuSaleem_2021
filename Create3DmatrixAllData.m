


clear all
DataFolder = ['X:\DATA\PROJECTS'];
FileNames = uigetfile_n_dir(DataFolder,'Select recording file');
load(FileNames{1,1});

trialSide_seconds = 1; % take 60 samples before and after trial
trialSide_samples = 60; % take 60 samples before and after trial
BonvisionFR = 60;
TrialsTotal = 0; UnitsTotal = 0;
for i = 1:size(ProjectData,1)
     sessionTrials(i) = size(ProjectData.Trial_data{i,1},1);
     sessionUnits = size(ProjectData.Units_Info{i,1},1);
     UnitsTotal = UnitsTotal + sessionUnits;
end
TrialTimeLength = length(ProjectData.Trial_data{1,1}.Speed_in_trial{1,1});
clear AM_UnitResponses AM_Speed AM_Param
AM_UnitResponses = repmat(nan,[max(sessionTrials), UnitsTotal, TrialTimeLength-1]);
AM_Speed = repmat(nan,[max(sessionTrials), UnitsTotal, TrialTimeLength]);
AM_Param = repmat(nan,[max(sessionTrials), UnitsTotal, 5]);

cellcount = 1; G_trialcount = 1;
Accel_Stimulus = 1;
for rr = 1:size(ProjectData,1) % scan through the sessions
    for i = 1:size(ProjectData.Units_Info{rr,1},1) % scan through all units
        trialcount = 1; p_i = 1;
        for j = 1:size(ProjectData.Trial_data{rr,1},1)
            SingletrialSpikes = ProjectData.Trial_data{rr,1}.SpikesTrial{j,1}{1,i};
            TrialStart = ProjectData.Session_data{rr,1}.ACInfo{1,1}.trialStartsEnds(j,1);
            TrialEnd = ProjectData.Session_data{rr,1}.ACInfo{1,1}.trialStartsEnds(j,2);
            edges = linspace(TrialStart-trialSide_seconds,TrialEnd+trialSide_seconds,TrialTimeLength); 
            if ~isempty(SingletrialSpikes)
                [SpikeCount_temp, Edges] = histcounts(SingletrialSpikes,edges);
            end
            AM_UnitResponses(j,cellcount,:) = SpikeCount_temp*BonvisionFR;            
            
            if Accel_Stimulus == 0
                % save information regarding trials
                % recording id
                AM_Param(j,cellcount,1) = rr;
                % direction
                AM_Param(j,cellcount,2) = ProjectData.Trial_data{rr,1}.Direction(j,1);
                % perturbation ON OFF
                AM_Param(j,cellcount,3) = ProjectData.Trial_data{rr,1}.PerturbationON(j,1);
                
                if ProjectData.Trial_data{rr,1}.PerturbationON(j,1)
                    % perturbation onset
                    AM_Param(j,cellcount,4) = ProjectData.Session_data{rr,1}.ACInfo{1,1}.PertStartsEnds(p_i,1)-TrialStart;
                    % perturbation offset
                    AM_Param(j,cellcount,5) = ProjectData.Session_data{rr,1}.ACInfo{1,1}.PertStartsEnds(p_i,2)-TrialStart;
                    p_i= p_i+1;
                end
            else
                % save information regarding trials
                % recording id
                AM_Param(j,cellcount,1) = rr;
                % starting TF
                AM_Param(j,cellcount,2) = ProjectData.Trial_data{rr,1}.Starting_TF(j,1);
                % ending TF
                AM_Param(j,cellcount,3) = ProjectData.Trial_data{rr,1}.Ending_TF(j,1);
                % dTF
                AM_Param(j,cellcount,4) = ProjectData.Trial_data{rr,1}.dTF(j,1);
                % index start of accel
                AM_Param(j,cellcount,5) = ProjectData.Trial_data{rr,1}.Accel_start_end(j,1);
                % index end of accel
                AM_Param(j,cellcount,6) = ProjectData.Trial_data{rr,1}.Accel_start_end(j,2);
            end
            
            AM_Speed(j,cellcount,1:length(ProjectData.Trial_data{rr,1}.Speed_in_trial{j,1})) = ProjectData.Trial_data{rr,1}.Speed_in_trial{j,1};
            
            AM_Time(j,cellcount,1:length(ProjectData.Trial_data{rr,1}.time_in_trial{j,1})) = ProjectData.Trial_data{rr,1}.time_in_trial{j,1};
                        
        end 
        cellcount = cellcount + 1;
    end
    rr
end


FileSep_i = strfind(FileNames{1},filesep);     
if Accel_Stimulus ==0
    % original data
    save([FileNames{1}(1:FileSep_i(end)) 'UnitsResponse.mat'], 'AM_UnitResponses', 'AM_Speed', 'AM_Param', 'AM_Time','-v7.3')
    % control data from last batch of three mice
    save([FileNames{1}(1:FileSep_i(end)) 'UnitsResponseCTRL.mat'], 'AM_UnitResponses', 'AM_Speed', 'AM_Param', 'AM_Time','-v7.3')
else
    % original data
    save([FileNames{1}(1:FileSep_i(end)) 'UnitsResponseAccelCTRL.mat'], 'AM_UnitResponses', 'AM_Speed', 'AM_Param', 'AM_Time','-v7.3')
end











