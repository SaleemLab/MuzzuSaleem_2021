% name of dat file
ProjectData.Session_data{1,1}.MetaData{1,1}.FileName

% cluster ID
ProjectData.Units_Info{1,1}.k_ID

INPUT
gwfparams.dataDir = '/path/to/data/';    % KiloSort/Phy output folder
gwfparams.fileName = 'data.dat';         % .dat file containing the raw 
gwfparams.dataType = 'int16';            % Data type of .dat file (this should be BP filtered)
gwfparams.nCh = 32;                      % Number of channels that were streamed to disk in .dat file
gwfparams.wfWin = [-40 41];              % Number of samples before and after spiketime to include in waveform
gwfparams.nWf = 2000;                    % Number of waveforms per unit to pull out
gwfparams.spikeTimes =    [2,3,5,7,8,9]; % Vector of cluster spike times (in samples) same length as .spikeClusters
gwfparams.spikeClusters = [1,2,1,1,1,2]; % Vector of cluster IDs (Phy nomenclature)   same length as .spikeTimes


 for rec = 1:size(ProjectData,1)
    clear SpikeTimes_def SpikeCluster_def
    Temp_i_r = ProjectData.Session_data{rec,1}.RecordingOI{1};
    StartTime_samples = ProjectData.Session_data{rec,1}.MetaData{1,1}.lims(Temp_i_r-1);
    for j = 1:size(ProjectData.Units_Info{rec,1},1)
        SpikeTimes_adj_temp = ProjectData.Units_Info{rec,1}.Spiketimes{j,1}*30000; 
        SpikeCluster_temp = single(ProjectData.Units_Info{rec,1}.k_ID{j,1}) * ones(size(SpikeTimes_adj_temp,1),1);
        if j == 1
            SpikeTimes_def = SpikeTimes_adj_temp;
            SpikeCluster_def = SpikeCluster_temp;
        else
            SpikeTimes_def = [SpikeTimes_def; SpikeTimes_adj_temp];
            SpikeCluster_def = [SpikeCluster_def; SpikeCluster_temp]; 
        end
    end
                   
    pos_sep = strfind(ProjectData.Session_data{rec,1}.MetaData{1,1}.FileName,filesep);
    if strfind(ProjectData.Session_data{rec,1}.MetaData{1,1}.FileName,'kilosort')
        gwfparams.dataDir = ProjectData.Session_data{rec,1}.MetaData{1,1}.FileName(1:pos_sep(end));
    else
        gwfparams.dataDir = [ProjectData.Session_data{rec,1}.MetaData{1,1}.FileName(1:pos_sep(end)) 'kilosort\'];
    end
    gwfparams.fileName = ProjectData.Session_data{rec,1}.MetaData{1,1}.FileName(pos_sep(end)+1:end);
    gwfparams.dataType = 'int16';
    gwfparams.nCh = 32;
    gwfparams.wfWin = [-60 61];              % -60 samples to +60 (30 samples = 1 ms)
    gwfparams.nWf = 500;                    % Number of waveforms per unit to pull out
    gwfparams.spikeTimes = SpikeTimes_def; % Vector of cluster spike times (in samples) same length as .spikeClusters
    gwfparams.spikeClusters = SpikeCluster_def;  % Vector of cluster IDs (Phy nomenclature)   same length as .spikeTimes
    
    clear SpikeTimes_def SpikeCluster_def
    
    wf = getWaveForms(gwfparams);
    clear gwfparams
    rec
    % save file every rec individually
    save(['Waveforms_' num2str(rec) '_f.mat'], 'wf' ,'-v7.3')
    clear wf
end



OUTPUT
wf.unitIDs                               % [nClu,1]            List of cluster IDs; defines order used in all wf.* variables
wf.spikeTimeKeeps                        % [nClu,nWf]          Which spike times were used for the waveforms
wf.waveForms                             % [nClu,nWf,nCh,nSWf] Individual waveforms
wf.waveFormsMean                         % [nClu,nCh,nSWf]     Average of all waveforms (per channel)
                                         % nClu: number of different clusters in .spikeClusters
                                         % nSWf: number of samples per waveform