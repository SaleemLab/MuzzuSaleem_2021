%% Tomaso Muzzu - get RF maps from LFP - 29 May 2020

%% Steps
% 1) Load sparse noise data from processed data
% 2) Load raw traces directly from OpenEphys
% 3) Filter the raw traces
% 4) Synchronise raw traces with Bonvision data
% 5) Compute the RF maps and save them in ProjectData

%% load data
if ~exist('ProjectData','var')
    [ProjectData AM_UnitResponses AM_Param AM_Speed AM_UOI SelectedResponses AM_UnitResponses_smooth] = LoadDataALL;
end
if filesep == '\' 
    addpath(['X:' filesep 'CODE' filesep 'DEV' filesep 'general' filesep 'SparseNoise']);
    addpath(['X:' filesep 'CODE' filesep 'STABLE' filesep 'OpenEphys_analysis_tools']);
else
    addpath(['/mnt/pfs09' filesep 'CODE' filesep 'DEV' filesep 'general' filesep 'SparseNoise']);
    addpath(['/mnt/pfs09' filesep 'CODE' filesep 'STABLE' filesep 'OpenEphys_analysis_tools']);
end

if ~isfield(ProjectData.SparseNoise{end,1},'LFP_1kHz')
    for i = 1:size(ProjectData,1)
        if ~isempty(ProjectData.SparseNoise{i,1})
            % load SN file directly from Processed Data
            [~,hostName] = system('hostname'); hostName = hostName(1:end-1);
            if ~strcmp(hostName, 'saleem12')
                DataFolder = 'X:\DATA\SUBJECTS';
            else
                DataFolder = 'X:\ibn-vision\DATA\SUBJECTS';
            end
            % 1) Load sparse noise data from processed data
            ProjectData.SparseNoise{i,1}.MetaData.FoldersList{ProjectData.SparseNoise{i,1}.ROI}
            FileName = uigetfile_n_dir([DataFolder filesep ProjectData.Mouse_name{i,1} filesep 'Processed'],'Select processed file of sparse noise');
            load(FileName{1});
            
            % 2) Load raw traces directly from OpenEphys
            FileList = dir([ES.MetaData.FoldersList{ES.ROI} filesep '101_CH*.continuous']);
            nChannels = 32;
            fileNum = zeros(1,nChannels);
            for iFile = 1:nChannels
                temp = FileList(iFile).name;
                fileNum(iFile) = str2num(temp(7:end-11));
            end
            [~, fileOrder] = sort(fileNum);
            % check how many channels were used and load them onto memory
            clear data
            data = zeros(length(ES.ACInfo.Timestamps),nChannels);
            for j = 1:length(FileList)
                data(:,j) = load_open_ephys_data_faster(fullfile(FileList(j).folder,FileList(fileOrder(j)).name));  % NB should preallocate data for speed
                % printing filename here as well for visual inspection -MM 2019-06
                fprintf('\nLoading channel %01d of %01d, file %01d of %01d\nFilename = %s\n',j,length(FileList),j,nChannels,FileList(fileOrder(j)).name);
            end
            
            % 3) save the LFP traces
            ProjectData.SparseNoise{i,1}.LFP_1kHz = resample(data,1,30);
            ProjectData.SparseNoise{i,1}.Time_1kHz = resample(ES.ACInfo.Timestamps,1,30);
            clear ES
        end
    end
save(['X:\DATA\PROJECTS\VisPerturbation' filesep 'AllDataHere_SN_LFP.mat'],'ProjectData','-v7.3');   
end
        
% clear SparseNoiseData
% k = 1;
% for i = 1:size(ProjectData,1)
%     if ~isempty(ProjectData.SparseNoise{i,1})
%         SparseNoiseData{k,1} = ProjectData.SparseNoise{i,1};     
%         k = k + 1;
%     end
% end
% save(['E:\OneDrive - University College London\SparseNoise_Haider' filesep 'SparseNoiseData.mat'],'SparseNoiseData','-v7.3');   
        
        
      
        
        
        
        
        
        
        
        
        
        
        
        