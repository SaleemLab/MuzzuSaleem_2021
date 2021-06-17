%% Tomaso Muzzu - UCL - 23 September 2019

%% Script to save data for visual perturbation altogheter in a single file

% load multiple sessions from one animal
% load one by one with a for loop

function UpdateAllDataFileAcceleration

clearvars -except StatsM
close all

FileNames = SelectSessions;
        
for i = 1:length(FileNames)
    % create temporary single session Struct/Table 
    ES_ = getES(FileNames(i));
    
%     figure
%     plot(ES_.Time,ES_.TF/max(ES_.TF))
%     hold on
%     plot(ES_.trialStartsEnds(:,1),ones(length(ES_.trialStartsEnds(:,1)),1)*0.5,'r.')
%     plot(ES_.trialStartsEnds(:,2),ones(length(ES_.trialStartsEnds(:,2)),1)*0.6,'k.')
%     
    [Array_, StimulusParams_] = ...
        GetSpeedDataAccelStim(ES_,...
        'Norm_mode','none', ... % 'MaxSpeedTrial', 'MaxSpeedSession', 'z_score', 'none'
        'log_speed',0, ... % 1 if computing speed in log scale
        'Smoothing',0.6); % in seconds, size of Gaussian filter for speed
    
    % remove OE channels and Time array from ES_.ACInfo
    ES_.ACInfo.Data = [];
    ES_.ACInfo.Timestamps = [];
    ES_.ACInfo.RotEncSpeed = [];
    ES_.ACInfo.RotEncPos = [];
      
    % save info regarding mouse name, recording nr and date
    pattern1 = 'SUBJECTS\';
    pattern2 = '\Processed';
    LetterStart = strfind(FileNames{i},pattern1);
    LetterEnd = strfind(FileNames{i},pattern2);
    MouseName = FileNames{i}(LetterStart+length(pattern1):LetterEnd-1);
    RecordingDate = ES_.ACInfo.ExpDate; 
    ProjectData = array2table(cell(1,1),'VariableNames',{'Mouse_name'});
    ProjectData.Mouse_name(1) = {MouseName};
    ProjectData.Recording_date(1) = {RecordingDate};
    
    % create table out of ES data
    Fields = fieldnames(ES_);
    clear SessionData
    SessionData = array2table(cell(1,size(fieldnames(ES_),1)-1),'VariableNames',Fields(1:end-1));
    for k = 1:size(fieldnames(ES_),1)-1 % the last one contains the spikes for each trial-
        SessionData{1,k} = {ES_.(Fields{k})};
    end
    % save('SessionData.mat','SessionData','-v7.3')
    % save table with session info in the first row of the table chunk
    ProjectData.Session_data(1) ={SessionData}; 
    
    % save spike responses of all cells for each trial
    for k = 1:size(ES_.SpikesTrial,2)
        Array_.SpikesTrial(k) = {ES_.SpikesTrial(:,k)'}; 
    end
    
    % save cell info into table
    UnitInfo = {'Spiketimes' 'UnitType' 'Shank' 'Channel' 'k_ID'};
    Cell_Data = array2table(cell(size(ES_.SpikeInfo,2)-1,size(ES_.SpikeInfo,1)-1),'VariableNames',UnitInfo);
    ES_.SpikeInfo(:,1);
    for k = 2:size(ES_.SpikeInfo,2)
        for j = 1:size(ES_.SpikeInfo,1)-1
            Cell_Data{k-1,j} = ES_.SpikeInfo(j,k);
        end
    end
    ProjectData.Units_Info(1) = {Cell_Data}; 
    
    ProjectData.Trial_data(1) = {Array_};
      
    % concatenate sessions in a single table
    if i == 1
        Array = ProjectData;
    else
        Array = [Array; ProjectData];
    end
    clear Array_ StimulusParams_ ES_ ProjectData
end

% here check if needs appending to an already saved file or creating a new file
% further check if I add a new session from an animals with sessions
% already saved.

ProjectFolder = uigetfile_n_dir('X:\DATA\PROJECTS','Select project folder');

if exist([ProjectFolder{1} filesep 'AllDataHere.mat'], 'file') == 2
    load([ProjectFolder{1} filesep 'AllDataHere.mat']);
    [Skip Array] = CheckFiles(ProjectData, Array);
    if ~Skip
        ProjectData = [ProjectData; Array];
        ProjectData = OrderArray(ProjectData);
        save([ProjectFolder{1} filesep 'AllDataHere.mat'], 'ProjectData', '-v7.3');
    else
        fprintf('Recordings already saved. File already up to date.\n')
    end
else % if file does not exit, create new file
    ProjectData = Array;
    save([ProjectFolder{1} filesep 'AllDataHere.mat'], 'ProjectData', '-v7.3');
end

end

function ProjectData = OrderArray(ProjectData)

clear RecOrder RecOrder_i
r = 1;
for ls = 1:size(ProjectData.Mouse_name,1)
    if ~isempty (ProjectData.Mouse_name{ls})
        RecOrder(r,1) = ProjectData.Mouse_name(ls);
        RecOrder(r,2) = ProjectData.Recording_date(ls);
        RecOrder_i(r) = ls;
        r = r + 1;
    end
end
[B, index] = sortrows(RecOrder);

clear NewArray TempArray
for i = 1:length(RecOrder_i)
    if i == 1
        NewArray = ProjectData(RecOrder_i(index(i)):RecOrder_i(index(i+1))-1, :);
    elseif index(i) == length(RecOrder_i)
        TempArray = ProjectData(RecOrder_i(index(i)):end, :);
        NewArray = [NewArray; TempArray];        
    else
        TempArray = ProjectData(RecOrder_i(index(i)):RecOrder_i(min(find(RecOrder_i>RecOrder_i(index(i)))))-1, :);
        NewArray = [NewArray; TempArray];
    end
end

ProjectData = NewArray;
end


function [Skip Array] = CheckFiles(ProjectData, Array)

clear RecOrder RecOrder_i RecOrder_A RecOrder_i_A
r = 1;
for ls = 1:size(ProjectData.Mouse_name,1)
    if ~isempty (ProjectData.Mouse_name{ls})
        RecOrder(r,1) = ProjectData.Mouse_name(ls);
        RecOrder(r,2) = ProjectData.Recording_date(ls);
        RecOrder_i(r) = ls;
        r = r + 1;
    end
end

r = 1;
for ls = 1:size(Array.Mouse_name,1)
    if ~isempty (Array.Mouse_name{ls})
        RecOrder_A(r,1) = Array.Mouse_name(ls);
        RecOrder_A(r,2) = Array.Recording_date(ls);
        RecOrder_i_A(r) = ls;
        r = r + 1;
    end
end
append_i = 0; c = 1;
for ls = 1:size(RecOrder_A,1)
    for lr = 1:size(RecOrder,1)
        if strcmp(RecOrder_A{ls,1},RecOrder{lr,1})
            if strcmp(RecOrder_A{ls,2},RecOrder{lr,2})
                append_is(c) = ls; % these are the recordings in common
                c = c+1;
                append_i = 1;
            end
        end
    end 
end
if append_i == 1 % if there recordings with the same name
    RecNew_i = RecOrder_i_A;
    RecNew_i(append_is) = [];
    clear TempArray NewArray
    % check if not all the new recordings have already been saved
    if length(RecNew_i)<length(RecOrder_i_A)
        for k = 1:length(RecNew_i)            
            if k == 1 & k == length(RecNew_i)
                NewArray = Array(RecNew_i(k):end, :);
            elseif k == 1
                NewArray = Array(RecNew_i(k):RecNew_i(min(find(RecNew_i>RecNew_i(k))))-1, :);
            elseif k == length(RecNew_i)
                TempArray = Array(RecNew_i(k):end, :);
                NewArray = [NewArray; TempArray];
            else
                TempArray = Array(RecNew_i(k):RecNew_i(min(find(RecNew_i>RecNew_i(k))))-1, :);
                NewArray = [NewArray; TempArray];
            end
        end
        Skip = 0; % do not skip the update of the file and copy the new array!
        Array = NewArray;
    else
        Skip = 1;
    end
    
else % if all recordings are new
    Skip = 0;
end


end