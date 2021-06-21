% Tomaso Muzzu - University College London - Neural Coding Lab 15/03/2014

%% Function to read files .kwik generated with the Klustasuite
% http://klusta.readthedocs.io/en/latest/kwik/

function ReadKwik

% Select the .kwik file containing spike times and cluster numbers
[FileName,PathName,FilterIndex] = uigetfile('Select .kwik');
FullPathKwik = [PathName FileName];

% Save all time_samples and clusters from the 4 shanks
NrShanks = 4;
for i = 0:NrShanks-1
    TimeSamples(i+1) = {h5read(FullPathKwik, ['/channel_groups/' num2str(i) '/spikes/time_samples'])};
    Clusters(i+1) = {h5read(FullPathKwik, ['/channel_groups/' num2str(i) '/spikes/clusters/main'])};
end

%% OBSOLETE. WAVEFORMS ARE NOT SAVED ANYMORE. CAN DO IT MANUALLY THOUGH!
% % Select the .kwx file containing waveforms and features
% FullPathKwx = FullPathKwik(1:end-1);
% FullPathKwx(end-2:end)= 'kwx';
% % Save all filtered waveforms from the 4 shanks
% for i = 0:NrShanks-1
%     WaveformsFilt(i+1) = {h5read(FullPathKwx, ['/channel_groups/' num2str(i) '/waveforms_filtered'])};
%     WaveformsRaw(i+1) = {h5read(FullPathKwx, ['/channel_groups/' num2str(i) '/waveforms_raw'])};
% end

%% TO BE CHECKED WITH AMAN WHETHER WE WANT TO DO THIS
% % Read text file with info about nice cells (ordered by shank and channel)
% [InfoFileName,PathInfoFile,FilterIndex] = uigetfile('Select .txt');
% filename = [ PathInfoFile '\' InfoFileName];
% SCUInfo = dlmread(filename, '\t', 7, 0);
% ChannelsPerShank = 8;
% for i=1:size(SCUInfo,1)
%     if size(SCUInfo,1)==1
%         SCUInfo(2) = ChannelsPerShank - mod(SCUInfo(1)*ChannelsPerShank,SCUInfo(2));
%         break
%     else
%         SCUInfo(i,2) = ChannelsPerShank - mod(SCUInfo(i,1)*ChannelsPerShank,SCUInfo(i,2));
%     end
% end


%% Save the following information from the selected cells
% Clusters with source cell, 
% TimeSamples with timing of spikes
% WaveformsFilt with filtered waveforms
% WaveformsRaw with raw waveforms

% One matlab cell with times of spikes, raw and filtered waveforms from
% main channel 
for i = 1:size(SCUInfo,1)
    if size(SCUInfo,1)==1
        ShankOfInterest = SCUInfo(1);
        ChannelOfInterest = SCUInfo(2);
        UnitID = SCUInfo(3);
        CellType = SCUInfo(4); % 1=single, 2=MUA
    else
        ShankOfInterest = SCUInfo(i,1);
        ChannelOfInterest = SCUInfo(i,2);
        UnitID = SCUInfo(i,3);
        CellType = SCUInfo(i,4);  % 1=single, 2=MUA
    end
    CurInfo = [ShankOfInterest,ChannelOfInterest+(ShankOfInterest*ChannelsPerShank), ...
                UnitID, CellType];
    % Get time stamps and clusters of all spikes in second shanks
    TimeStampsID = cell2mat(TimeSamples(1,ShankOfInterest));
    ClusterNR = cell2mat(Clusters(1,ShankOfInterest));
    % Get time stamps of spike of cell of interest UnitID
    TimeStampsOI = (TimeStampsID(find(ClusterNR==UnitID)));
    % Get waveforms from shank of interest 
    AllWaveformsFilt = cell2mat(WaveformsFilt(1,ShankOfInterest));
    AllWaveformsRaw = cell2mat(WaveformsRaw(1,ShankOfInterest));
    % Get waveforms only from channel of interest
    WaveformsFiltID = AllWaveformsFilt(ChannelOfInterest,:,ClusterNR==UnitID);
    WaveformsRawID = AllWaveformsRaw(ChannelOfInterest,:,ClusterNR==UnitID);
    % Store the info of all cells in a single multifiedl cell
    Data{i}= {CurInfo,TimeStampsOI,WaveformsFiltID,WaveformsFiltID};
end

% Reorganise Data so that first level of cells are the shanks
for i = 1:NrShanks
    j = 1;
    for k = 1:length(Data)    
        if SCUInfo(k,1)==i
            DataS{i}{j}=Data{k};
            j=j+1;
        end
    end
end

% Find not empty cells in EPhys_Data
for i  = 1:length(DataS)
    R(i) = ~isempty(DataS{1,i});
end
DataS = DataS(find(R));
% The cell structure 'DataS' has as many sub-cells as the number of shanks;
% each subcells has as many neurons as were selected. In the 3rd level, 
% each sub-cell has 4 cells that are:
% Data{1,ShanksNR}{1,NeuronNR}(1,1) = array of 4 elements with [ShankOfInterest,ChannelOfInterest,UnitID, CellType];
% Data{1,ShanksNR}{1,NeuronNR}(1,2) = array with time stamps of every spike of the neuron UnitID
% Data{1,ShanksNR}{1,NeuronNR}(1,3) = array with filtered waveforms ordered as timestamps
% Data{1,ShanksNR}{1,NeuronNR}(1,4) = array with raw waveforms ordered as timestamps of same index

% Save these data into a .mat file easily accessible from Matlab at later
% times
SavingPath = [PathName(1:end-7) 'C_Matlab_Analysis\'];
prompt = '\n\n\n\nInsert number of experiment please\n';
ExpNumber = input(prompt,'s');
save([SavingPath 'EPhys_Data__' ExpNumber '.mat'], 'DataS');

end