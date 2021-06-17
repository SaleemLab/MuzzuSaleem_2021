%% Tomaso Muzzu - UCL - 12/12/2018

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
FileName = uigetfile_n_dir(DataFolder,'Select CSV file');
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
% 2) Get analog signals
FileSep_i = strfind(FileName{1},filesep);
EphysDataFolder = [FileName{1}(1:FileSep_i(end-2)) 'ePhys'];
ACInfo = OE_get_AC_Signals(EphysDataFolder);
ACInfo = OE_Analyse_PD_Signal(ACInfo);
ACInfo = OE_Compute_Speed_Signal(ACInfo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% 3) Synchronise Bonsai data with the timing of OpenEphys
clear TrialExs TwoDiffContrast
TwoDiffContrast = diff(diff(ES.Contrast));
TrialExs = find((TwoDiffContrast/max(TwoDiffContrast))>0.1);
display('verify signal from Bonvision log')
% find the wrong trial onset/offsets
j = 1; k = 1; clear TrialExs_i
for i = 1:length(TrialExs)
    % trial onsets
    if TrialExs(i)+400<length(ES.Contrast)
        if sum(ES.Contrast(TrialExs(i)+1:TrialExs(i)+420)>0)>=360
            if i == 1
                TrialExs_i(j,1) = TrialExs(i)+1;
                j = j+1;
            else
                if TrialExs(i)-TrialExs_i(j-1,1)>450
                    TrialExs_i(j,1) = TrialExs(i)+1;
                    j = j+1;
                end
            end
        end
    end
%     % trial offsets
%     if i>=2
%         if TrialExs(i)>420
%             if sum(ES.Contrast(TrialExs(i)-420:TrialExs(i)+1)>0)>=380
%                 if k==1
%                     TrialExs_i(k,2) = TrialExs(i)+1;
%                     k = k+1;
%                 else
%                     if TrialExs(i)-TrialExs_i(k-1,2)>450
%                         TrialExs_i(k,2) = TrialExs(i)+1;
%                         k = k+1;
%                     end
%                 end
%             end
%         end
%     end
end

for i = 2:length(TrialExs_i)
    % trial offsets
    TrialExs_i(i-1,2) = TrialExs(max(find(TrialExs<TrialExs_i(i,1)))-1);
end
TrialExs_i(end,2) = TrialExs(end);

% figure
% plot(ES.Time,ES.Contrast/max(ES.Contrast));
% hold on
% plot(ES.Time,ES.TF/max(ES.TF))
% hold on
% plot(ES.Time(3:end),TwoDiffContrast/max(TwoDiffContrast))
% hold on
% plot(ES.Time(TrialExs(1:2:end)+1),ones(length(TrialExs(1:2:end)),1)*0.4,'or')
% hold on
% plot(ES.Time(TrialExs_i(:,1)),ones(length(TrialExs_i(:,1)),1)*0.3,'ok') % trial onsets
% hold on
% plot(ES.Time(TrialExs_i(:,2)),ones(length(TrialExs_i(:,2)),1)*0.3,'*k') % trial offsets

DiffTF = diff(ES.TF);
PertExs_i = [find(DiffTF<-1), find(DiffTF>1)];
% remove the trials in which pert happens during zero contrast trials
ZeroContr_i = unique(ES.Trial(ES.TF<=0.1 & ES.Contrast<=0.1));
PertExs_i_t = unique(ES.Trial(ES.TF==0));
for i = length(ZeroContr_i):-1:1
    PertExs_i(find(PertExs_i_t==ZeroContr_i(i)),:) = [];
end
% [MinTrialDur MinTrialDur_i] = min(TrialExs_i(:,2)-TrialExs_i(:,1));
% if MinTrialDur<mean((TrialExs_i(:,2)-TrialExs_i(:,1)))*0.5
%     fprintf('\nPlease check the trial on/offsets are correctly detected!\n');
%     % overalapping of one start with an end of trial
%     TrialExs = [TrialExs(1:MinTrialDur_i*2-1); TrialExs(MinTrialDur_i*2+1:end)];
%     TrialExs_i = [TrialExs(1:2:end-1)+1, TrialExs(2:2:end)+1];
% end
% 
% [MinPertDur MinPertDur_i] = min(PertExs_i(:,2)-PertExs_i(:,1));
% if MinPertDur<mean((PertExs_i(:,2)-PertExs_i(:,1)))*0.5
%     fprintf('\nPlease check the perturbation on/offsets are correctly detected!\n');
% end
    
% plot data to verify that on/offsets of trials and perturbation are
% correctly detected:
% downsamle timestamps and photodiode signals
VerifySignals(ACInfo,ES,TrialExs_i,PertExs_i); 

%%
ES.trialStartsEnds = TrialExs_i;
ES.PertStartsEnds = PertExs_i;

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
    save([SavingDir{1,1} filesep FileName{1}(FileSep_i(end-3)+1:FileSep_i(end-2)-1) '_BV_AC_Spikes_ks_' ES.ACInfo.ExpDate '.mat'], 'ES', '-v7.3');
else
    save([SavingDir{1,1} filesep SavingDir{1,1}(end-15:end-10) '_BV_AC_Spikes_' ES.ACInfo.ExpDate '.mat'], 'ES', '-v7.3');
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



