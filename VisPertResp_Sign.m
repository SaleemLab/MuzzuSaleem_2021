%% Tomaso Muzzu - UCL - 07/10/2019

% statistical test for significance of vis. stim. response 

clear all
DataFolder = ['X:\DATA\PROJECTS'];
FileNames = uigetfile_n_dir(DataFolder,'Select recording file');
load(FileNames{1,1});

% smooth the FR
BonvisionFR = (ProjectData.Session_data{1,1}.Time{1,1}(end)-ProjectData.Session_data{1,1}.Time{1,1}(1))/length(ProjectData.Session_data{1,1}.Time{1,1});
GaussFilter_W = 0.3;
Width = round(GaussFilter_W/BonvisionFR);
Sigma = Width/3; % standard deviation in number of samples (converted from time in seconds)
x_g = linspace(-Width/2, Width/2, Width);
gaussFilter = exp(-x_g.^2/(2*Sigma^2));
gaussFilter_ = gaussFilter / sum (gaussFilter); % normalize

% 1 - select the spikes happening in this recording
% 2 - smooth the 60Hz spike count with a gaussian filter
% 3 - select the responses to the first 4 seconds of vis. stim.
% 4 - select the FR during the Inter-Trial Intervals
% 5 - compute the mean FR, median FR, std of FR, of 4s and ITI
% 6 - compute the DM, MI, peak during 4s, integral of the mean curve
% 7 - repeat step from 2 to 6 1k times. Each time shifts the smoothed
% FR by 20-60 seconds and move the later spikes to the beginning
tic
clear FR_vis FR_UOI_VIS
for i = 1:size(ProjectData,1) % loop through sessions
    trialSide_seconds = 1; 
     
    %% 1 - select the spikes happening in this recording
    RecordingOI = ProjectData.Session_data{i,1}.RecordingOI{:};
    SessionStart = sum(ProjectData.Session_data{i,1}.MetaData{1,1}.lims(1:RecordingOI-1))/ProjectData.Session_data{i,1}.SpikeInfo{1,1}{end,2};
    SessionEnd = sum(ProjectData.Session_data{i,1}.MetaData{1,1}.lims(1:RecordingOI))/ProjectData.Session_data{i,1}.SpikeInfo{1,1}{end,2};
    
    Time_edges = linspace(SessionStart, SessionEnd, (SessionEnd-SessionStart)*60); % 60 is sampling frequency.
    
    %% 2 - smooth the 60Hz spike count with a gaussian filter
    clear SessionSpikes
    for cs = 1:size(ProjectData.Units_Info{i,1},1)
        SessionSpikes_Unit = ProjectData.Units_Info{i,1}.Spiketimes{cs,1}(...
            ProjectData.Units_Info{i,1}.Spiketimes{cs,1}>=SessionStart & ...
            ProjectData.Units_Info{i,1}.Spiketimes{cs,1}<=SessionEnd);
        SessionSpikes_Unit = histcounts(SessionSpikes_Unit,Time_edges);
        SessionSpikes(cs,:) = conv(SessionSpikes_Unit, gaussFilter_, 'same');
    end
    SessionSpikes = SessionSpikes.*60;
    
    %% 3 - select the responses to the first 4 seconds of vis. stim.
    Time_edges = linspace(SessionStart, SessionEnd, (SessionEnd-SessionStart)*60); % 60 is sampling frequency.
    Time_edges = Time_edges-min(Time_edges);
    clear FR_4s
    for j = 1:size(ProjectData.Trial_data{i,1},1)
        [v ind_m(1)] = min(abs(Time_edges-...
                        (ProjectData.Session_data{i,1}.ACInfo{1,1}.trialStartsEnds(j,1)-trialSide_seconds)));
        %[v ind_m(2)] = min(abs(Time_edges-...
        %               (ProjectData.Session_data{i,1}.ACInfo{1,1}.trialStartsEnds(j,2)+trialSide_seconds)));
        FR_4s(j,:,:) = SessionSpikes(:,ind_m(1):ind_m(1)+(trialSide_seconds+4*trialSide_seconds)*60-1);
    end
    
    %% 4 - select the FR during the Inter-Trial Intervals
    clear FR_ITI
    FR_ITI = zeros(size(FR_4s,2),1);
    for j = 1:size(ProjectData.Trial_data{i,1},1)-1
        [v ind_m(1)] = min(abs(Time_edges-...
                        (ProjectData.Session_data{i,1}.ACInfo{1,1}.trialStartsEnds(j,2))));
        [v ind_m(2)] = min(abs(Time_edges-...
                        (ProjectData.Session_data{i,1}.ACInfo{1,1}.trialStartsEnds(j+1,1))));
        FR_ITI = [FR_ITI SessionSpikes(:,ind_m(1):ind_m(2))];
    end
    FR_ITI = FR_ITI(:,2:end);
    
    %% 5 & 6 - compute the mean FR, median FR, std of FR, of 4s and ITI 
    clear VSR
    for k = 1:size(FR_4s,2)
        FR_mean = squeeze(mean(FR_4s(:,k,:)));
        % mean for first 4 seconds
        VRM_stim = mean(FR_mean(trialSide_seconds*60:end));
        VRM_prestim = mean(FR_mean(1:trialSide_seconds*60));
        VSR(k,1) = VRM_stim;
        VSR(k,2) = VRM_prestim;
        VSR(k,3) = (VRM_stim-VRM_prestim)/(VRM_stim+VRM_prestim) ; % DM
        VSR(k,4) = (VRM_stim)/(VRM_stim+VRM_prestim) ; % MI
    end
        
    FR_gray = mean(FR_ITI,2);
            
    %% 7 - repeat step from 2 to 6 1k times. Each time shifts the smoothed
    % FR by 20-60 seconds and move the later spikes to the beginning
    %tic
    for p_i = 1:2
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % shuffle the firing rate of every cell and re-compute FR over the first 4 seconds
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if p_i == 1
            clear FR_vis
            % mean FR of first 4 seconds of stim
            FR_vis(1,p_i,:) = VSR(:,1); 
            % mean FR of 1s before stim onset
            FR_vis(2,p_i,:) = VSR(:,2); 
            % DM
            FR_vis(3,p_i,:) = VSR(:,3); 
            % MI
            FR_vis(4,p_i,:) = VSR(:,4); 
            % avg baseline FR
            FR_vis(5,p_i,:) = FR_gray;
        elseif p_i>=2
            % take all the spiketimes and add an arbitrary time btw 20 and 60
            % seconds. If the spiketimes fall over the end of the recording
            % wrap them around.
            clear SessionSpikes_sh
            Time_shift = round(20/BonvisionFR+rand*(60/BonvisionFR));
            SessionSpikes_sh = [SessionSpikes(:,Time_shift+1:end) SessionSpikes(:,1:Time_shift)];
            
            clear FR_4s
            for j = 1:size(ProjectData.Trial_data{i,1},1)
                [v ind_m(1)] = min(abs(Time_edges-...
                    (ProjectData.Session_data{i,1}.ACInfo{1,1}.trialStartsEnds(j,1)-trialSide_seconds)));
                %[v ind_m(2)] = min(abs(Time_edges-...
                %               (ProjectData.Session_data{i,1}.ACInfo{1,1}.trialStartsEnds(j,2)+trialSide_seconds)));
                FR_4s(j,:,:) = SessionSpikes_sh(:,ind_m(1):ind_m(1)+(trialSide_seconds+4*trialSide_seconds)*60-1);
            end
            
            % 4 - select the FR during the Inter-Trial Intervals
            clear FR_ITI
            FR_ITI = zeros(size(FR_4s,2),1);
            for j = 1:size(ProjectData.Trial_data{i,1},1)-1
                [v ind_m(1)] = min(abs(Time_edges-...
                    (ProjectData.Session_data{i,1}.ACInfo{1,1}.trialStartsEnds(j,2))));
                [v ind_m(2)] = min(abs(Time_edges-...
                    (ProjectData.Session_data{i,1}.ACInfo{1,1}.trialStartsEnds(j+1,1))));
                FR_ITI = [FR_ITI SessionSpikes_sh(:,ind_m(1):ind_m(2))];
            end
            FR_ITI = FR_ITI(:,2:end);
            
            % 5 & 6 - compute the mean FR, median FR, std of FR, of 4s and ITI
            clear VSR
            for k = 1:size(FR_4s,2)
                FR_mean = squeeze(mean(FR_4s(:,k,:)));
                % mean for first 4 seconds
                VRM_stim = mean(FR_mean(trialSide_seconds*60:end));
                VRM_prestim = mean(FR_mean(1:trialSide_seconds*60));
                VSR(k,1) = VRM_stim;
                VSR(k,2) = VRM_prestim;
                VSR(k,3) = (VRM_stim-VRM_prestim)/(VRM_stim+VRM_prestim) ; % DM
                VSR(k,4) = (VRM_stim)/(VRM_stim+VRM_prestim) ; % MI
            end            
            FR_gray = mean(FR_ITI,2);
            
            % mean FR of first 4 seconds of stim
            FR_vis(1,p_i,:) = VSR(:,1); 
            % mean FR of 1s before stim onset
            FR_vis(2,p_i,:) = VSR(:,2); 
            % DM
            FR_vis(3,p_i,:) = VSR(:,3); 
            % MI
            FR_vis(4,p_i,:) = VSR(:,4); 
            % avg baseline FR
            FR_vis(5,p_i,:) = FR_gray;
        end 
    end
    %toc
    if i == 1
        FR_UOI_VIS = FR_vis;
    else
        FR_UOI_VIS = cat(3,FR_UOI_VIS,FR_vis);
    end
    i 
end
toc
% save info re sign. modulated units.
AM_UOI_VIS = FR_UOI_VIS;

save('DM_visresp_shuffled.mat', 'FR_UOI_VIS','-v7.3');

%% FR_VIS is 3 matrix where 
% 1st dim => 5 metrics of visual response
% 2nd dim => 1 actual metric + 1000 computed from shifted data. 
% 3rd dim => 1493 units recorded over 37 recordings.
% The 5 metrics for measuring visual responses are
% 1) mean FR of first 4 seconds of stim
% 2) mean FR of 1s before stim onset
% 3) DM
% 4) MI
% 5) avg baseline FR
FR_VisResp = squeeze(FR_VIS);

figure
plot(FR_VIS(3,:),'.')
hold on
plot(find(RespUnits(:,1)),FR_VIS(3,RespUnits(:,1)),'*g')




