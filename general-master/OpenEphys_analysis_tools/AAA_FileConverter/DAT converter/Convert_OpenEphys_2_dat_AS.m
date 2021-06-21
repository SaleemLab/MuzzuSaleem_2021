% function Convert_OpenEphys_2_dat(MouseName, ExpDate, Recordings)
%
% Openephys opener and dat converter 
%
% Workflow of the program:
% 1 - select files to concatenate
% 2 - load the files into memory, channel per channel
% 3 - correct timestamps if needed
% 4 - re-map channels following the specs of the probe used
% 5 - rearrange data and save as .dat file for spike sorting
%
% History
%  Tomaso Muzzu - UCL 27/02/2018: wrote it
%  SGS and SZ - UCL 29/10/2019: added capacity to pass through recordings rather than user select

function Convert_OpenEphys_2_dat_AS(MouseName, ExpDate, Recordings, OutputExpName)
tic
% Some defaults
% % % AS added for new recordings
ePhysRecordingChannelRoot = '103_CH*.continuous';
nChannels = 32;

%%
% 1 - select data to convert and path to save
% Set up path to this subject/exp date if known
SubjectFolder = fullfile('X:','DATA','SUBJECTS'); % nb: may want to make root folder flexible for differnet network mappings
% If we also know the recording list already, can ignore user input
if exist('Recordings','var') && ~isempty(Recordings)
    FoldersList = cell(1,length(Recordings));
    for thisRecording = 1:length(Recordings)
        FoldersList{1,thisRecording} = fullfile(SubjectFolder,MouseName,'ePhys',ExpDate,Recordings{thisRecording});
    end
else
    % If we know mouse/expdate but not session, can set up
    if exist('MouseName','var') && ~isempty(MouseName) && exist('ExpDate','var') && ~isempty(ExpDate)
        SubjectFolder = fullfile(SubjectFolder,MouseName,'ePhys',ExpDate);
    end
    FoldersList = uigetfile_n_dir(SubjectFolder); % nb: will start from root if mouse/expdate are not known
end
% Specify folder where to save files
if ~exist('OutputExpName','var') || isempty(OutputExpName)
    fprintf('Select target folder for saving .dat files please. \n\n\n');
    TargetPath = uigetdir(SubjectFolder);
    prompt = 'Insert name of recording please. File chunks will be ordered by number \n';
    ExpName = input(prompt,'s');
    OutFileName = fullfile(TargetPath,ExpName);
    OutFileName = [OutFileName '.dat'];
else
    if exist('MouseName','var') && ~isempty(MouseName) && exist('ExpDate','var') && ~isempty(ExpDate) && exist('Recordings','var') && ~isempty(Recordings)
        OutFileName = fullfile(SubjectFolder,MouseName,'ePhys',ExpDate,'kilosort',sprintf('%s_%s_%s',MouseName, ExpDate,OutputExpName));
        OutFileName = [OutFileName '.dat'];
    end
end

%%% 
% OK, now do everything
for i = 1:length(FoldersList)
    %% 2 - load the files into memory, channel per channel
    FileList = dir(fullfile(FoldersList{i},ePhysRecordingChannelRoot));
    fileNum = zeros(1,nChannels);
    for iFile = 1:nChannels
        temp = FileList(iFile).name;
        fileNum(iFile) = str2num(temp(7:end-11));
    end
    [~, fileOrder] = sort(fileNum);
    % check how many channels were used and load them onto memory
   
    for j = 1:length(FileList)
        data(:,j) = load_open_ephys_data_faster(fullfile(FileList(j).folder,FileList(fileOrder(j)).name));  % NB should preallocate data for speed
        % printing filename here as well for visual inspection -MM 2019-06
        fprintf('\nLoading channel %01d of %01d, file %01d of %01d\nFilename = %s\n',j,length(FileList),i,length(FoldersList),FileList(fileOrder(j)).name); 
    end
    
    % save the lengths of the current recording file
    lims(i) = length(data);
    
    %% 3 - remap the channels (depends on the probe used)
     if i == 1
        ProbeBrand = input('Did you use a CamNeuroTech E series probe? Please type 1 for YES \n','s')
     end
    Ch_map = remap_32ch(str2num(ProbeBrand));
    
    Ch_map' % printing chan map here for visual inspection -MM 2019-06
    
    data_ch_sorted = data(:,Ch_map);
    clear data
    
    %% 4 - linearise data to 1D array
    Channels = int16(data_ch_sorted);
    clear data_ch_sorted
    
    % Compute length of vector for single channel. To save recording as:
    % [Time1Channel1, Time1Channel2, .. Time1ChannelN, Time2Channel1, ...]
    V_Length = size(Channels,1) * size(Channels,2);
    % Reshape matrix into array
    Channels_1Darray = reshape(Channels',1,V_Length);
    clear Channels
    
    %% 4 - save the ephys data as a dat file at int16 precision 
    % Save the dat file
    save2dat(Channels_1Darray, i, OutFileName) 
    clear Channels_1Darray
    
end    

% Save the meta data to a matlab file - AS 19/03/2018
matFileName = strcat(OutFileName, '_meta.mat');
save(matFileName, 'lims', 'FoldersList', 'Ch_map', 'OutFileName');


timepast = toc;
fprintf(['\n\nFINISHED: It took ' num2str(timepast) ' seconds. \n\nHappy Clustering!!!\n\n\n']);

end
    
    