%% Tomaso Muzzu - UCL - 16/01/2018
%  function to get ES data from Bonvision data

function ES = getES(FileName)

load(FileName{1,1})

% 1) find start of trials
for i = 1:length(ES.trialStartsEnds)
    if isfield(ES,'PertStartsEnds')
        if ~isempty(ES.PertStartsEnds(ES.Time(ES.PertStartsEnds(:,1))>ES.Time(ES.trialStartsEnds(i,1)) & ...
                ES.Time(ES.PertStartsEnds(:,1))<ES.Time(ES.trialStartsEnds(i,2)),1)) % it is a trial with perturbation?
            ES.PertTrial(i) = 1;
        else
            ES.PertTrial(i) = 0;
        end
    end
%     if sum(ES.Orientation==180)>1
%         ES.OrienTrial(i) = round(mean(ES.Orientation(ES.trialStartsEnds(i,1):ES.trialStartsEnds(i,2))));
%     else
%         ES.OrienTrial(i) = round(mean(ES.Orientation(ES.trialStartsEnds(i,1):ES.trialStartsEnds(i,2))/(2*pi)*360));
%     end
end
% 2) compute position and speed in cm from the Bonvision data
Pos_ = ES.Wheel;
jumps = find(abs(diff(Pos_))>60000);
if ~isempty(jumps)
    clear Position
    for i=1:length(jumps)
        if i == 1
            Position(1:jumps(i)) = Pos_(1:jumps(i));
        else
            Position(jumps(i-1)+1:jumps(i)) = Pos_(jumps(i-1)+1:jumps(i)) + Pos_(jumps(i-1))*(i-1);
        end
    end
    Position(jumps(end)+1:length(Pos_)) = Pos_(jumps(end)+1:end) + Pos_(jumps(i))*(i);
else
    Position = Pos_';
end
% Encoder resolution is 1024 ticks per revolution;
Encoder_Res = 1024;
% radius of the wheel is 9 cm;
Wheel_Circum = 8.9*2*pi;
ES.MousePosition = Position'*Wheel_Circum/Encoder_Res;
ES.MouseSpeed = diff(ES.MousePosition)./diff(ES.Time);
ES.MouseSpeed(abs(ES.MouseSpeed)>100) = 0; 
clear Position Pos_ jumps i
if isfield(ES,'MetaData')
    % 3) get the spikes of the recording of interest
    for i = 1:length(ES.MetaData.FoldersList)
        if strcmp(ES.MetaData.FoldersList{i}(end-18:end),FileName{1,1}(end-22:end-4))
            RecordingOI = i;
        end
    end
    clear SessionSpikes
    TimeStart = sum(ES.MetaData.lims(1:RecordingOI-1)/ES.SpikeInfo{end,2});
    TimeEnd = sum(ES.MetaData.lims(1:RecordingOI))/ES.SpikeInfo{end,2};
    for i = 2:length(ES.SpikeInfo)
        if RecordingOI>1
            SessionSpikes{i-1} = (ES.SpikeInfo{1,i}(ES.SpikeInfo{1,i}>TimeStart & ES.SpikeInfo{1,i}<TimeEnd)) - TimeStart;
        else
            SessionSpikes{i-1} = (ES.SpikeInfo{1,i}(ES.SpikeInfo{1,i}>0 & ES.SpikeInfo{1,i}<TimeEnd)) - TimeStart;
        end
    end
    for i = 1:length(SessionSpikes)
        for j = 1:length(ES.ACInfo.trialStartsEnds)
            SpikeTimes = SessionSpikes{i};
            SpikesTrial{i,j} = SpikeTimes(SpikeTimes>ES.ACInfo.trialStartsEnds(j,1)-1 & SpikeTimes<ES.ACInfo.trialStartsEnds(j,2)+1);
        end
    end
    ES.RecordingOI = RecordingOI;
    ES.SpikesTrial = SpikesTrial;
end


end