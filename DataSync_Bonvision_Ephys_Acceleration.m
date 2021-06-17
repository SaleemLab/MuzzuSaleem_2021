%% Tomaso Muzzu - UCL - 23/09/2019

%% Main script to generate a file with synchronised data of visual stimuli, behaviour and ephys
% 1) Get data from csv files saved from Bonsai
% 2) Get photodiode signal
% 3) Sync it with the timestamps of Bonsai and save these info
% 4) Save spikes, single unit info and compute FR

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% 1) Get data from csv files saved from Bonsai
clear all
close all
CurrentFolder = cd; 
%DataFolder = [CurrentFolder(1:end-(length(CurrentFolder)-3)) 'Archive - saleemlab' filesep 'Data' ];
DataFolder = 'X:\DATA\SUBJECTS';
FileName = uigetfile_n_dir(DataFolder,'Select CSV file for TFdTF stimulus');
FileName
VisStimLog = csvread(FileName{1,1},1,0);
if sum(sign(diff(VisStimLog(:,11))))<0 % adjustment for jan-march '19 when the rotary encoder was moved on the right.
    VisStimLog(:,11) = -VisStimLog(:,11)+max(VisStimLog(:,11));
end
fid = fopen(FileName{1,1});
Columns = textscan(fid,'%s',1);
fclose(fid);
VisStimLog_Header = strsplit(Columns{1,1}{1,1},',');
for i = 1:length(VisStimLog_Header)
    ES.(VisStimLog_Header{i}) = VisStimLog(:,i);
end
ES.Time = ES.Time-min(ES.Time);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
% 2) Get analog signals
FileSep_i = strfind(FileName{1},filesep);
EphysDataFolder = [FileName{1}(1:FileSep_i(end-2)) 'ePhys'];
ACInfo = OE_get_AC_Signals(EphysDataFolder);
ACInfo = OE_Analyse_PD_Signal_Acceleration(ACInfo);
ACInfo = OE_Compute_Speed_Signal(ACInfo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% 3) Synchronise Bonsai data with the timing of OpenEphys
% resample and normalise signal
if size(ACInfo.Data,2)>5
    normPD = ACInfo.Data(:,4)/max(ACInfo.Data(:,4));
else
    normPD = ACInfo.Data(:,2)/max(ACInfo.Data(:,2));
end
plot_TS_temp = resample(ACInfo.Timestamps,1,30);
plot_PD_temp = resample(normPD,1,30);

% plot PD signal under Temporal Frequency time series to verify alignment
figure
plot(plot_TS_temp(plot_TS_temp>ACInfo.trialStartsEnds(1,1) & plot_TS_temp<ACInfo.trialStartsEnds(end,2))-ACInfo.trialStartsEnds(1,1),...
     plot_PD_temp(plot_TS_temp>ACInfo.trialStartsEnds(1,1) & plot_TS_temp<ACInfo.trialStartsEnds(end,2)))
hold on
plot(ACInfo.trialStartsEnds(:,1)-ACInfo.trialStartsEnds(1,1),ones(length(ACInfo.trialStartsEnds(:,1)),1)*0.5,'sr')
hold on
plot(ACInfo.trialStartsEnds(:,2)-ACInfo.trialStartsEnds(1,1),ones(length(ACInfo.trialStartsEnds(:,1)),1)*0.5,'sk')
hold on
plot(ES.Time,ES.TF/max(ES.TF))

TrialExs_i = [ACInfo.trialStartsEnds(:,1)-ACInfo.trialStartsEnds(1,1) ,...
    ACInfo.trialStartsEnds(:,2)-ACInfo.trialStartsEnds(1,1)];

% simply convert to bonvision time units the starts and stops of the trial
ES.trialStartsEnds = TrialExs_i;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
% 4) Save spikes, single unit info 
clear SpikeData
kilosort = 1;
if kilosort == 1
    SpikeData = get_ks_Spikes(EphysDataFolder);
else
    SpikeData = getSpikes(EphysDataFolder);
end

%% 
% 5) Save altogether
ES.SpikeInfo = SpikeData.SpikeInfo;
ES.MetaData = SpikeData.MetaData;
ES.ACInfo = ACInfo;

% Save all these in a file (this will be added with the spikes later)
SavingDir = uigetfile_n_dir(FileName{1}(1:FileSep_i(end-1)),'Choose save directory'); 
if kilosort == 1
    save([SavingDir{1,1} filesep FileName{1}(FileSep_i(end-3)+1:FileSep_i(end-2)-1) '_TFdTF_AC_Spikes_ks_' ES.ACInfo.ExpDate '.mat'], 'ES', '-v7.3');
else
    save([SavingDir{1,1} filesep SavingDir{1,1}(end-15:end-10) '_TFdTF_AC_Spikes_' ES.ACInfo.ExpDate '.mat'], 'ES', '-v7.3');
end

clear all
close all
                                                

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MANUAL CORRECTIONS

% Manual correction for mouse TM_03, second recording.
ES.Time(TrialExs_i(59:61,:))-ES.Time(TrialExs_i(1,1))
TrialExs_i(61:end+1,:)=TrialExs_i(60:end,:);
[v, ind] = min(abs(ES.Time-ES.Time(TrialExs_i(1,1))-(trialON(60)-trialON(1))));
ES.Time(ind)-ES.Time(TrialExs_i(1,1))
TrialExs_i(60,1) = ind;
[v, ind] = min(abs(ES.Time-ES.Time(TrialExs_i(1,1))-(trialOFF(60)-trialON(1))));
ES.Time(ind)-ES.Time(TrialExs_i(1,1))
TrialExs_i(60,2) = ind;



