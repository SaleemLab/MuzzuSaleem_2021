%% Tomaso Muzzu - UCL 27/02/2018

% Openephys opener and dat converter 

% Workflow of the program:
% 1 - select files to concatenate
% 2 - load the files into memory, channel per channel
% 3 - correct timestamps if needed
% 4 - re-map channels following the specs of the probe used
%    ***added option to choose channels here*** --2019-05 MM
%       ex. ChanOpt=[1:16]
% 5 - rearrange data and save as .dat file for spike sorting

function Convert_OpenEphys_2_dat_channel_option(MouseName, ExpDate, ChanOpt)
tic

%% 1 - select data to convert
if nargin == 2
    SubjectFolder = ['X:\DATA\SUBJECTS\' MouseName '\ePhys\' ExpDate];
else
    SubjectFolder = 'X:\DATA\SUBJECTS\';
end
fprintf('Please select all the open ephys folders to concatenate. \n\n\n');
FoldersList = uigetfile_n_dir(SubjectFolder);

for i = 1:length(FoldersList)
    %% 2 - load the files into memory, channel per channel
    FileList = dir([FoldersList{i} filesep '101_CH*.continuous']);
    % check how many channels were used and load them onto memory
    for j = 1:length(FileList)
        [data(:,j), timestamps, info] = load_open_ephys_data_faster([FileList(j).folder filesep FileList(j).name]);
        fprintf(['\nLoading channel ' num2str(j) ' of ' num2str(length(FileList)) ', file ' num2str(i) ' of ' num2str(length(FoldersList)) '.\n']);
    end
    clear FileList timestamps info
    
    % save the lengths of the current recording file
    lims(i) = length(data);
    
    %% 3 - remap the channels (depends on the probe used)
    if i == 1
        ProbeBrand = input('Please type 1 for NeuroNexus or 2 for CamNeuroTech \n','s')
    end
    Ch_map = remap_32ch(str2num(ProbeBrand));
    data_ch_sorted = data(:,Ch_map);
    clear data
    
    %% 4 - linearise data to 1D array
    Channels = int16(data_ch_sorted);
    clear data_ch_sorted
    % Compute length of vector for single channel. To save recording as:
    % [Time1Channel1, Time1Channel2, .. Time1ChannelN, Time2Channel1, ...]
    
    % select only the channels chosen by user (2019-05 MM)
    Channels=Channels(:,ChanOpt);
    
    V_Length = size(Channels,1) * size(Channels,2);
    % Reshape matrix into array
    Channels_1Darray = reshape(Channels',1,V_Length);
    clear Channels
    
    %% 4 - save the ephys data as a dat file at int16 precision 
    % Select folder where to save files
    if i == 1
        fprintf('Select target folder for saving .dat files please. \n\n\n');
        TargetPath = uigetdir(SubjectFolder);
        prompt = 'Insert name of recording please. File chunks will be ordered by number \n';
        ExpName = input(prompt,'s');
        FileName = [TargetPath filesep ExpName  '.dat'];
    end
    save2dat(Channels_1Darray, i, FileName) 
    clear Channels_1Darray
    i
end    

% Save the meta data to a matlab file - AS 19/03/2018
matFileName = strcat(FileName, '_meta.mat');
save(matFileName, 'lims', 'FoldersList', 'Ch_map', 'FileName');

timepast = toc;
fprintf(['\n\nFINISHED: It took ' num2str(timepast) ' seconds. \n\nHappy Clustering!!!\n\n\n']);

end
    
    