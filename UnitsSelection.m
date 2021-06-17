%% Tomaso Muzzu - UCL - 12/09/2019

% initial selection on best and/or desired units


%% select best units
function AM_UOI = UnitsSelection(ProjectData,AM_UnitResponses,varargin)

BonvisionFR = 60;
GaussFilter_W = 0.3; % seconds
Width = round(GaussFilter_W*BonvisionFR);
Sigma = Width/3; % standard deviation in number of samples (converted from time in seconds)
x_g = linspace(-Width/2, Width/2, Width);
gaussFilter = exp(-x_g.^2/(2*Sigma^2));
gaussFilter_ = gaussFilter / sum (gaussFilter); % normalize


if nargin == 3 % only FR threshold applied
     
    trialSide_seconds = 1; % take 60 samples before and after trial
    trialSide_samples = 60; % take 60 samples before and after trial
    cellcount = 1; clear MeanUnitFR
    for rr = 1:size(ProjectData,1)
        % define Time series for current experiment
        RecordingOI = ProjectData.Session_data{rr,1}.RecordingOI{:};
        SessionStart = sum(ProjectData.Session_data{rr,1}.MetaData{1,1}.lims(1:RecordingOI-1))/ProjectData.Session_data{rr,1}.SpikeInfo{1,1}{end,2};
        SessionEnd = sum(ProjectData.Session_data{rr,1}.MetaData{1,1}.lims(1:RecordingOI))/ProjectData.Session_data{rr,1}.SpikeInfo{1,1}{end,2};
        
        %SessionStart = ProjectData.Session_data{rr,1}.ACInfo{1,1}.trialStartsEnds(1,1);
        %SessionEnd = ProjectData.Session_data{rr,1}.ACInfo{1,1}.trialStartsEnds(end,2);
        Time_edges = linspace(SessionStart, SessionEnd, (SessionEnd-SessionStart)*BonvisionFR); % 60 is sampling frequency.
        for i = 1:size(ProjectData.Units_Info{rr,1},1)
            UnitSpikesDurSession{cellcount,1} = ProjectData.Units_Info{rr,1}.Spiketimes{i,1}(...
                ProjectData.Units_Info{rr,1}.Spiketimes{i,1}>=SessionStart & ...
                ProjectData.Units_Info{rr,1}.Spiketimes{i,1}<=SessionEnd);
            clear SpikeCount_temp SpikeCount
            [SpikeCount_temp, Edges] = histcounts(UnitSpikesDurSession{cellcount,1},Time_edges);
            
            SpikeCount = SpikeCount_temp*BonvisionFR; % to get Hz
            UnitFR = conv(SpikeCount, gaussFilter_, 'same');
            MeanUnitFR(cellcount,1) = mean(UnitFR); % mean spiking rate
            % Reliability measure: FR is the same between first and
            % last 3rd or recording.
            Time_1_3 = round(length(Time_edges)/3);
            MeanUnitFR(cellcount,2) = mean(UnitFR(1:Time_1_3)); % mean of spiking count in first 3rd of recording
            MeanUnitFR(cellcount,3) = std(UnitFR(1:Time_1_3)); % std of spiking count in first 3rd of recording
            MeanUnitFR(cellcount,4) = mean(UnitFR(Time_1_3*2:end)); % mean of spiking count in last 3rd of recording
            MeanUnitFR(cellcount,5) = std(UnitFR(Time_1_3*2:end)); % std of spiking count in last 3rd of recording
            MeanUnitFR(cellcount,6) = std(UnitFR(Time_1_3+1:Time_1_3*2)); % std of spiking count in middle 3rd of recording
            MeanUnitFR(cellcount,7) = std(UnitFR(Time_1_3+1:Time_1_3*2)); % std of spiking count in middle 3rd of recording
            cellcount = cellcount + 1;  
        end
    end
    clear AM_UOI
    UOI(:,1) = MeanUnitFR(:,1)>=varargin{1}; %FR_Thres
    % the FR is different more than a STD
    UOI(:,2) = MeanUnitFR(:,4)>=max((MeanUnitFR(:,2)-MeanUnitFR(:,3)),0) & MeanUnitFR(:,4)<=(MeanUnitFR(:,2)+MeanUnitFR(:,3));
    % the FR either at the beginning or at the end is <varargin{1}
    UOI(:,3) = MeanUnitFR(:,4)>=varargin{1} & MeanUnitFR(:,2)>=varargin{1} & MeanUnitFR(:,6)>=varargin{1} ;
    % consider the first and last spike times when computing responses 
    
    
    AM_UOI = UOI(:,1) & UOI(:,2) & UOI(:,3);
else
    
   
    
end
    

end
