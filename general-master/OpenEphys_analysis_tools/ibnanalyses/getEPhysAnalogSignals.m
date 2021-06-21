% Analog channel reader. Loads the photodiode, speed, and sync pulse signals for one recording
% 
% History
% Tomaso Muzzu - UCL - 05/04/2018: Wrote it
% SGS and SZ 29/10/2019: adjusted to allow fileanme pass through

function ACInfo = getEPhysAnalogSignals(MouseName,ExpDate,Recording)
global DIRS
ePhysAnalogChannelRoot = '101_ADC*.continuous';

% Select the folders where the OpenEphys data are stored
% If we also know the recording list already, can ignore user input
if exist('MouseName','var') && ~isempty(MouseName) && exist('ExpDate','var') && ~isempty(ExpDate) && exist('Recording','var') && ~isempty(Recording)
    FoldersList = {fullfile(DIRS.Subjects,MouseName,'ePhys',ExpDate,Recording)};
else
    %FoldersList = uigetfile_n_dir(DIRS.multichanspikes,'Select the folder(s) of interest');
    FoldersList = uigetfile_n_dir(DIRS.ePhys,'Select the single folder of interest'); %SZ 29/10/19 updated with new folder name
end

% Now get the files
FileList = dir(fullfile(FoldersList{1},ePhysAnalogChannelRoot));

% check how many channels were used and load them onto memory
for j = 1:length(FileList)
    % Actually loads the analog channel
    [data(:,j), timestamps(:,j), info] = load_open_ephys_data_faster([FileList(j).folder filesep FileList(j).name]);
    fprintf(['\nLoading channel ' num2str(j) ' of ' num2str(length(FileList)) '.']);
    ACInfo.AnChannelsOE(j) = str2num(info.header.channel(end));  % physical channel used
end
fprintf('\n')
ACInfo.Data = data;
ACInfo.Timestamps = timestamps(:,j);
ACInfo.SamplingRateOE = info.header.sampleRate;    % sampling rate in Hz
% check if there is more than one recording in the file (appended)
SepIndexes = find(diff(ACInfo.Timestamps)>1);
if ~isempty(SepIndexes)
    ACInfo.SessionStarts = SepIndexes;
end

% %%%%%%%% Keep sampling rate and 'channel' as metadata recorded into
% the ACInfo
% data(:,1) = sync pulse signal
% data(:,2) = photodiode signal
% data(:,3) = signal A from rotary encoder
% data(:,4) = signal B from rotary encoder





