%% Tomaso Muzzu - UCL - 23/09/2019

%% detect PD changes from PD signal acquired with OpenEhpys 
% version for the TF acceleration stimulus

%% photodiode signal: save original signal downsampled at 1kHz, timestamps of stimulus onsets and offsets
function [ACInfo] = OE_Analyse_PD_Signal_Acceleration(ACInfo,directory)

% data(:,1) = sync pulse signal
% data(:,2) = photodiode signal
% data(:,3) = signal A from rotary encoder
% data(:,4) = signal B from rotary encoder
% data(:,5) = Bonvision state signal

% normalise signal
if size(ACInfo.Data,2)>5
    normPD = ACInfo.Data(:,4)/max(ACInfo.Data(:,4));
else
    normPD = ACInfo.Data(:,2)/max(ACInfo.Data(:,2));
end
[temp_pks, temp_locs] = findpeaks(normPD);
GMModel = fitgmdist(normPD,2);
ThresholdPeaks = max(GMModel.mu);
pks = temp_pks(temp_pks>ThresholdPeaks);
locs = temp_locs(temp_pks>ThresholdPeaks);
clear temp_pks temp_locs
% verify position of peaks
%         figure
%         plot(ACInfo.Data(:,2)/max(ACInfo.Data(:,2)));
%         hold on
%         plot(locs,pks,'ro');

% find first peaks signalling the onset of the white square
PksTimeDifference = diff(locs);

temp_VSOnsets = find(PksTimeDifference>0.1*ACInfo.SamplingRateOE);
VSOnsetsIndecesOFF = [locs(1); locs(temp_VSOnsets+1)]; % save indexes of stimulus onsets
% find last peaks signalling the offset of the white square
VSOnsetsIndecesON = [locs(temp_VSOnsets-1); locs(end)]; % save indexes of stimulus offsets
stimOFF = ACInfo.Timestamps(VSOnsetsIndecesON(1:end-1));
stimON = ACInfo.Timestamps(VSOnsetsIndecesOFF(2:end));
trialOFF = sort(stimOFF);
trialON = sort(stimON);

% verify the events are detected correctly:
fprintf('Please verify that all trial and perturbation onsets/offsets are detected correctly\n');
plot_TS_temp = resample(ACInfo.Timestamps,1,30);
plot_PD_temp = resample(normPD,1,30);
figure
plot(plot_TS_temp,plot_PD_temp)
hold on
plot(trialON,ones(1,length(trialON))*0.35,'r.','MarkerSize',18)
hold on
plot(trialOFF,ones(1,length(trialOFF))*0.4,'k.','MarkerSize',18)

ACInfo.GrayScreen = [trialOFF(1) , trialON(1) ; trialOFF(end) , trialON(end)]
ACInfo.trialStartsEnds = [trialON(2:end-1), trialOFF(3:end)];

end




