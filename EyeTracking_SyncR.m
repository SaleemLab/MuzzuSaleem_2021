%% Tomaso Muzzu - UCL - 04 / 05 / 2020

% load eye tracking data and sync it with ephys

% only for naive mice

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
% 2) Get Bonvision Data
if ~exist('ProjectData','var')
    [ProjectData AM_UnitResponses AM_Param AM_Speed AM_UOI SelectedResponses AM_UnitResponses_smooth] = LoadDataALL;
end
cellcount = 1;
for i = 1:size(ProjectData,1)
    % get recording info
    FileNR = ProjectData.Session_data{i,1}.RecordingOI{1};
    ProjectData.Session_data{i,1}.MetaData{1,1}.FoldersList{FileNR}
    Filesep_i = strfind(ProjectData.Session_data{i,1}.MetaData{1,1}.FoldersList{FileNR},filesep);
    DataFolder = ProjectData.Session_data{i,1}.MetaData{1,1}.FoldersList{FileNR}(1:Filesep_i(4)-1);
    % load CSV file
    FileName = uigetfile_n_dir(DataFolder,'Select CSV file');
    FileName
    VisStimLog = csvread(FileName{1,1},1,0);
    fid = fopen(FileName{1,1});
    Columns = textscan(fid,'%s',1);
    fclose(fid);
    VisStimLog_Header = strsplit(Columns{1,1}{1,1},',');
    for j = 1: length(VisStimLog_Header)
        Dot_i = strfind(VisStimLog_Header{1,j},'.');
        if length(Dot_i)==2
            VisStimLog_Header{j} = [VisStimLog_Header{1,j}(Dot_i(1)+1:Dot_i(2)-1) '_' VisStimLog_Header{1,j}(Dot_i(2)+1:end)];
        elseif length(Dot_i)==1
            VisStimLog_Header{j} = VisStimLog_Header{1,j}(Dot_i(1)+1:end);
        elseif isempty(Dot_i)
            VisStimLog_Header{j} = 'TimeOfDay';
        end
    end
    % coielate the two aiays of time of the day
    T2sync_eye = VisStimLog(:,7);
    T2sync_dTF = ProjectData.Session_data{i,1}.TimeOfDay{1,1};
    T2sync_eye_lims = [min(T2sync_eye) max(T2sync_eye)];
    T2sync_dTF_lims = [min(T2sync_dTF) max(T2sync_dTF)];
    TrimInfo(1) = sign(T2sync_eye_lims(1)-T2sync_dTF_lims(1));
    TrimInfo(2) = sign(T2sync_eye_lims(2)-T2sync_dTF_lims(2));
    if TrimInfo(1)<0
        [v Trim(1)] = min(abs(T2sync_eye-T2sync_dTF_lims(1))); % trim eyetracking signal
    else
        [v Trim(1)] = min(abs(T2sync_dTF-T2sync_eye_lims(1))); % trim dTf signal
    end
    if TrimInfo(2)<0
        [v Trim(2)] = min(abs(T2sync_dTF-T2sync_eye_lims(2))); % trim dTF signal
    else
        [v Trim(2)] = min(abs(T2sync_eye-T2sync_dTF_lims(2))); % trim  eyetracking signal
    end 
    if TrimInfo(1)<0    
        T2sync_eye = T2sync_eye(Trim(1):end);
    else
        T2sync_dTF = T2sync_dTF(Trim(1):end);
    end
    if TrimInfo(2)<0
        T2sync_dTF = T2sync_dTF(1:Trim(2));
    else
        T2sync_eye = T2sync_eye(1:Trim(2));
    end
    T2sync_eye_lims = [min(T2sync_eye) max(T2sync_eye)];
    T2sync_dTF_lims = [min(T2sync_dTF) max(T2sync_dTF)];
    TimeDiff = diff(T2sync_eye_lims)-diff(T2sync_dTF_lims)
    % interpolate eyetracking data based on dTF timestamps
    if diff(T2sync_eye_lims)-diff(T2sync_dTF_lims)>50
        printf('Alignment issues');
    else
        % make eye tracking data as long as the DTF signal
        ET_interp = interp1(T2sync_eye, ...
                            VisStimLog(find(VisStimLog(:,end)-T2sync_eye_lims(1)==0):...
                                        find(VisStimLog(:,end)-T2sync_eye_lims(2)==0),:),...
                            T2sync_dTF);
        ET_interp(end,:) = [];
        % create 3D aiay of eye tracking data as you did for speed
        trialSide_seconds = 1; % take 60 samples before and after trial
        trialSide_samples = 60; % take 60 samples before and after trial
        for ll = 1:size(ProjectData.Units_Info{i,1},1) % scan through all units
            trialcount = 1; p_i = 1;
            for j = 1:size(ProjectData.Trial_data{i,1},1)
                TrialStart = ProjectData.Session_data{i,1}.trialStartsEnds{1,1}(j,1);
                TrialEnd = ProjectData.Session_data{i,1}.trialStartsEnds{1,1}(j,2);                 
                TempEyeTracking = interp1(1:length(ET_interp(max(TrialStart-60,1):min(TrialEnd+60,length(ET_interp)),:)), ...
                                                ET_interp(max(TrialStart-60,1):min(TrialEnd+60,length(ET_interp)),:), ...
                                                1:560);
                AM_EyeTracking(j,cellcount,1:length(TempEyeTracking),:) = TempEyeTracking;
            end
            cellcount = cellcount + 1;
        end
        
    end  
    % create table with eye tracking data and save it in ProjectData
    ProjectData.EyeTracking{i} = array2table(ET_interp,'VariableNames',VisStimLog_Header);
end

'X:\DATA\PROJECTS\VisPerturbation\'

save(['X:\DATA\PROJECTS\VisPerturbation\' 'UnitsResponseEyeCTRL.mat'], 'AM_UnitResponses', 'AM_Speed', 'AM_Param','AM_EyeTracking','-v7.3')
ProjectDataCTRL = ProjectData;
save(['X:\DATA\PROJECTS\VisPerturbation\' 'AllDataHereEyeCTRL.mat'],'ProjectDataCTRL','-v7.3')

