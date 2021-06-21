%% Tomaso Muzzu - UCL - 10/04/2018

%% photodiode signal: save original signal downsampled at 1kHz, timestamps of stimulus onsets and offsets
function [stimON, stimOFF, ACInfo] = getPhotoDiodeSignal(animal, iseries, iexp, photoChannel,verbose)
global DIRS
if ~exist('photoChannel','var') || isempty(photoChannel)
    photoChannel=2; %SZ 29/10/19 Specify here the analog channel for the photodiode signal NB was 2 but SZ needs 4
end
if ~exist('verbose','var') || isempty(verbose)
    verbose=1; %include plots
end
% Build path to folder
DIRname  = fullfile(DIRS.ePhys,num2str(iseries),num2str(iexp));
% load the mat file containing stim ON and OFF events recorded with OE

if (nargout<3) & exist([DIRname filesep 'AC_Info.mat'],'file')
    load(fullfile(DIRname,'AC_Info.mat'),'stimON','stimOFF');
    % else load the analog channels now, find when the timestamps of the VS
    % onsets and offsets
else
    fprintf('\nAnalogue channels have not been converted yet or some error reading the saved file. Please select folder containing relevant OE recording for recovering ACs info.');
    ACInfo = getEPhysAnalogSignals(animal, iseries, iexp);
    % data(:,1) = sync pulse signal
    % data(:,2) = photodiode signal
    % data(:,3) = signal A from rotary encoder
    % data(:,4) = signal B from rotary encoder
    
    
    % resample and normalise signal
    normdata = ACInfo.Data(:,photoChannel)/max(ACInfo.Data(:,photoChannel));
    
    % Plot the PD signal --added 21/08/18 MM
    if verbose
        figure;
        subplot(211); plot(ACInfo.Data(:,photoChannel))
        subplot(212); plot(normdata)
    end
    
    
    [temp_pks,temp_locs] = findpeaks(normdata);
    pks = temp_pks(temp_pks>0.5);
    locs = temp_locs(temp_pks>0.5);
    clear temp_pks temp_locs
    % find first peaks signalling the onset of the white square
    PksTimeDifference = diff(locs);
    temp_VSOnsets = find(PksTimeDifference>900);
    VSOnsetsIndecesON = [locs(1); locs(temp_VSOnsets+1)]; % save indexes of stimulus onsets
    % find last peaks signalling the offset of the white square
    VSOnsetsIndecesOFF = [locs(temp_VSOnsets-1); locs(end)]; % save indexes of stimulus offsets
    if ~isempty(find(diff(VSOnsetsIndecesON)>3*ACInfo.SamplingRateOE)) || ~isempty(find(diff(VSOnsetsIndecesOFF)>3*ACInfo.SamplingRateOE))
        warning('Error during photodiode signal extraction. Please check it manually [SGS and SZ - not sure what this is trying to pick up]');
    end
    ACInfo.stimON = ACInfo.Timestamps(VSOnsetsIndecesON)-min(ACInfo.Timestamps);
    ACInfo.stimOFF = ACInfo.Timestamps(VSOnsetsIndecesOFF)-min(ACInfo.Timestamps);
    stimON = ACInfo.stimON;
    stimOFF = ACInfo.stimOFF;
    % save Onsets and Offsets into relevant folder
    save(fullfile(DIRname,'AC_Info.mat'),'-struct', 'ACInfo');
end

% verify detection of stimulus presentation
%     figure
%     plot(ACInfo.Timestamps-min(ACInfo.Timestamps),ACInfo.Data(:,2))
%     hold on
%     plot(ACInfo.Timestamps(VSOnsetsIndecesON)-min(ACInfo.Timestamps),ACInfo.Data(VSOnsetsIndecesON,2),'s')
%     hold on
%     plot(ACInfo.Timestamps(VSOnsetsIndecesOFF)-min(ACInfo.Timestamps),ACInfo.Data(VSOnsetsIndecesOFF,2),'s')
%
end

