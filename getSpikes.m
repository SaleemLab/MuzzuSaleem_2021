function SpikeData = getSpikes(directory)

    ePhysData = uigetfile_n_dir(directory,'select Kwik file file');
        
    FileNameKwik = ePhysData{1,1};
    clear spkTimes spkClus cellIDs ncells
    for i = 0:1
        spkTimes{i+1} = hdf5read(FileNameKwik, ['/channel_groups/' num2str(i) '/spikes/time_samples']);
        spkClus{i+1} = hdf5read(FileNameKwik, ['/channel_groups/' num2str(i) '/spikes/clusters/main']);
        cellIDs(i+1,1:length(unique(spkClus{i+1}))) = unique(spkClus{i+1});
        ncells(i+1)  = length(unique(spkClus{i+1}));
    end
    cellIDs=cellIDs(:,2:end);
    ncells = ncells-1;
    clear Channel_Info        
    for i_shank = 0:1
        temp =  h5info(FileNameKwik,  ['/channel_groups/' num2str(i_shank) '/channels']);
        for i_channel=0:length(temp.Groups)-1
            TempChannel = h5info(FileNameKwik, ['/channel_groups/' num2str(i_shank) '/channels/' num2str(i_channel+(i_shank*(length(temp.Groups))))]);
            Channel_Info{i_channel+1+(i_shank*(length(temp.Groups)-1)),1} = TempChannel.Attributes(1).Value;
            Channel_Info{i_channel+1+(i_shank*(length(temp.Groups)-1)),2} = TempChannel.Attributes(2).Value;
        end
    end
    temp =  h5info(FileNameKwik,  ['/recordings/0']);
    sampleRate = double((temp.Attributes(3).Value));
    for i = 1:length(spkTimes)
        spkTimes{i} = double(spkTimes{i});
    end
    
    noise_list = [];
    clear SpikeInfo
    SpikeInfo{1,1} = 'Spiketimes';
    SpikeInfo{2,1} = 'Unit type';
    SpikeInfo{3,1} = 'Shank';
    SpikeInfo{4,1} = 'Channel';
    SpikeInfo{5,1} = 'ID';
    SpikeInfo{6,1} = 'Frame rate'; SpikeInfo{6,2} = sampleRate; 
    j=2;
    for i = 1:length(spkTimes)
        TimeStamps = cell2mat(spkTimes(1,i));
        spkCluster = cell2mat(spkClus(1,i));
        for icell = 1:ncells(i)
            temp = h5info(FileNameKwik, ['/channel_groups/' num2str(i-1) '/clusters/main/' num2str(cellIDs(i,icell)) '/']);
            for idx = 1:length(temp.Attributes)
                if strcmp(temp.Attributes(idx).Name, 'cluster_group')
                    ilabel = temp.Attributes(idx).Value;
                    break
                end
                idx = idx + 1;
            end
            ilabel = temp.Attributes(idx).Value;
            temp = h5info(FileNameKwik, ['/channel_groups/' num2str(i-1) '/cluster_groups/main']);
            labelType = eval(['temp.Groups(' num2str(ilabel+1) ').Attributes(end).Value']);
            if ~strcmp(labelType,'Noise')    
                SpikeInfo{1,j} = TimeStamps(find(spkCluster==cellIDs(i,icell))); % spike times
                SpikeInfo{2,j} = labelType; % unit type
                SpikeInfo{3,j} = i; % shank
                SpikeInfo{4,j} = 0; % channel
                SpikeInfo{5,j} = cellIDs(i,icell); % ID
                SpikeInfo{6,j} = sampleRate;
                j = j+1;
            end
        end
    end
   
    %% load metadata info
    MetaData = load([FileNameKwik(1:end-4) 'dat_meta.mat']);
    
    %% save SpikeData
    SpikeData.SpikeInfo = SpikeInfo;  
    SpikeData.MetaData = MetaData;
    
    %save([char(directory) filesep 'SpikeData.mat'], 'SpikeData', '-v7.3');
end

