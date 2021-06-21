%% Tomaso Muzzu - UCL 27/02/2018
% 2019-06 AS modified to fix channel ordering
% 2019-10 MM modified to take multiple probes 
%             -input channel allocation for each probe as vector (e.g.ProbeChans=[1:32 ; 33:64])
%             -input save file names for each probe as string cell array (e.g.SaveFileSuffix={'CA1','V1'})

% Openephys opener and dat converter 

% Workflow of the program:
% 1 - select files to concatenate
% 2 - load the files into memory, channel per channel
% 3 - correct timestamps if needed
% 4 - re-map channels following the specs of the probe used
% 5 - rearrange data and save as .dat file for spike sorting

function Convert_OpenEphys_2_dat_multiprobe(FolderPath, MouseName, ExpDate, ProbeChans, SaveFileSuffix)
tic

%% 1 - select data to convert

% if nargin > 2
%     SubjectFolder = ['X:\DATA\SUBJECTS\' MouseName '\ePhys\' ExpDate];
% else
%     SubjectFolder = 'X:\DATA\SUBJECTS\';
% end

SubjectFolder=fullfile(FolderPath,MouseName,'ePhys',ExpDate);

FoldersList = uigetfile_n_dir(SubjectFolder);

for i = 1:length(FoldersList)
    %% 2 - Identify number of probes
%     if i==1
%     ProbeNum = input('How many probes did you use?');
%     ProbeChans=input('Which channels belong to which probe? ex.[1:32 ; 33:64]');
%     end
    
    FileList = dir([FoldersList{i} filesep '101_CH*.continuous']);
    for ii=1:length(FileList)
    FileName{ii}=FileList(ii).name;
    end
    FileNameSorted = natsortfiles(FileName);
    
   for iprobe=1:size(ProbeChans,1)
    %% 3 - load the files into memory, channel per channel
    
%     for iFile = 1:length(ProbeChans(iprobe,:))
%         temp = FileList(ProbeChans(iprobe,iFile)).name;
%         fileNum(iFile) = str2num(temp(7:end-11));
%     end
%     [~, fileOrder] = sort(fileNum);
    
    % check how many channels were used and load them onto memory
    for j = 1:length(ProbeChans(iprobe,:))
        ch=ProbeChans(iprobe,j); % current channel
%         [data(:,j), timestamps, info] = load_open_ephys_data_faster([FileList(j).folder filesep FileList(fileOrder(j)).name]);
        [data(:,j), timestamps, info] = load_open_ephys_data_faster([FileList(ch).folder filesep FileNameSorted{ch}]);

        % printing filename here as well for visual inspection -MM 2019-06
%         fprintf(['\nLoading channel ' num2str(j) ' of ' num2str(length(FileList)) ', file ' num2str(i) ' of ' num2str(length(FoldersList)) '. \nFilename = '  FileList(fileOrder(j)).name  '\n' ]); 
        fprintf(['\nLoading channel ' num2str(ch) ' of ' num2str(length(FileList)) ', file ' num2str(i) ' of ' num2str(length(FoldersList)) '. \nFilename = '  FileNameSorted{ch}  '\n' ]); 
    end
    
    clear timestamps info
    
    % save the lengths of the current recording file
    lims(i) = length(data);
    
    %% 4 - remap the channels (depends on the probe used)
%      if i == 1
%         ProbeBrand = input('Did you use a CamNeuroTech E series probe? Please type 1 for YES \n','s');
%      end
    ProbeBrand='1';
    Ch_map = remap_32ch(str2num(ProbeBrand));
    Ch_map' % printing chan map here for visual inspection -MM 2019-06
    %data2use=data(:,ProbeChans(iprobe,:));
    %data_ch_sorted = data2use(:,Ch_map);
    data_ch_sorted = data(:,Ch_map);
    clear data2use
    %% 5 - linearise data to 1D array
    Channels = int16(data_ch_sorted);
    clear data_ch_sorted
    % Compute length of vector for single channel. To save recording as:
    % [Time1Channel1, Time1Channel2, .. Time1ChannelN, Time2Channel1, ...]
    V_Length = size(Channels,1) * size(Channels,2);
    % Reshape matrix into array
    Channels_1Darray = reshape(Channels',1,V_Length);
    clear Channels
    
    %% 6 - save the ephys data as a dat file at int16 precision 
    % Select folder where to save files
%     if i == 1
%         fprintf('Select target folder for saving .dat files please. \n\n\n');
%         TargetPath = uigetdir(SubjectFolder);
%         prompt = ['Insert name of recording for probe ' num2str(iprobe) ' please. File chunks will be ordered by number \n'];
%         ExpName = input(prompt,'s');
%         FileNames{iprobe} = [TargetPath filesep ExpName  '.dat'];
%     end
    ExpName = [MouseName '_' ExpDate '_' SaveFileSuffix{iprobe}];
    FileNames{iprobe} = [SubjectFolder filesep ExpName  '.dat'];
    disp(['saving to ' FileNames{iprobe} ])
    save2dat(Channels_1Darray, i, FileNames{iprobe}) 
    clear Channels_1Darray
   end
   clear data FileList
   disp(['finished processing folder' num2str(i)])
end    

for iprobe=1:size(ProbeChans,1)
% Save the meta data to a matlab file - AS 19/03/2018
matFileName{iprobe} = strcat(FileNames{iprobe}, '_meta.mat');
Filename=FileNames{iprobe};
save(matFileName{iprobe}, 'lims', 'FoldersList', 'Ch_map', 'Filename');
end

timepast = toc;
fprintf(['\n\nFINISHED: It took ' num2str(timepast) ' seconds. \n\nHappy Clustering!!!\n\n\n']);

end
    
    