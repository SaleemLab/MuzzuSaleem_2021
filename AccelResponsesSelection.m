%% Tomaso Muzzu - UCL - 24/09/2019

% initial selection on best and/or desired units


function SelectedResponses = AccelResponsesSelection(ProjectData,AM_UnitResponses_smooth,AM_Param,AM_Speed,AM_UOI,Cond)

% conditions = 1 --> nr of recording
% conditions = 2 --> TF at stimulus onset
% conditions = 3 --> TF at simutlus offset
% conditions = 4 --> dTF shown
% conditions = 5 --> index of when Tf starts changing
% conditions = 6 --> index of when Tf finishes changing

% look at perturbation responses regardless of running or grating
% direction
if size(AM_UOI,2)==1
    SelectedCells = AM_UOI;
else
    SelectedCells = AM_UOI(:,1) & AM_UOI(:,2);
end

% find all the TF's at the start of the trials
StarTF_2dmask = squeeze(AM_Param(:,:,2));
start_TF_values = unique(StarTF_2dmask(~isnan(StarTF_2dmask(:))));

% find all the TF's at the end of the trials
EndingTF = squeeze(AM_Param(:,:,3));
End_T = unique(EndingTF(~isnan(EndingTF(:))));

% separate the responses for the different acceleration values
dTF_2Dmask = squeeze(AM_Param(:,:,4));
dTF_values = unique(dTF_2Dmask(~isnan(dTF_2Dmask(:))));

% find trials in which dTF lasts 0.5s and 1s
dTF_dur_2D = (squeeze(AM_Param(:,:,6))-squeeze(AM_Param(:,:,5)))/60; % 60 is the sampling frequency
dTF_long_idx = dTF_dur_2D>0.8;

switch Cond
    case 0
        %% all units together, distinguishing for accel. values
        clear Response_OI_10 Response_OI_05 TrialsResp_max_10 TrialsResp_max_05 SelectedResponses
        for i = 1 : length(dTF_values)
            Accel_Con = dTF_2Dmask==dTF_values(i);
            % look at trials with dTF lasting 1 second
            TrialsSelection_10 = Accel_Con & dTF_long_idx;
            TrialsSelection_05 = Accel_Con & ~dTF_long_idx;
            % generate 3D matrix for including time responses
            t1 = double(repmat(TrialsSelection_10, 1, 1, size(AM_UnitResponses_smooth,3))); % repeating the selection across time
            t2 = double(repmat(TrialsSelection_05, 1, 1, size(AM_UnitResponses_smooth,3))); % repeating the selection across time
            % nan all the zero values
            t1(t1(:)==0) = nan; t2(t2(:)==0) = nan;
            % select only the responses during selected trials
            TrialsResp_OI_10 = t1.*AM_UnitResponses_smooth;
            TrialsResp_OI_05 = t2.*AM_UnitResponses_smooth;
            % compute the mean responses for each unit during trials at specified
            % dTF
            TrialsResp_OI_10_2D = squeeze(nanmean(TrialsResp_OI_10,1));
            TrialsResp_OI_05_2D = squeeze(nanmean(TrialsResp_OI_05,1));
            % normalise the responses of interest
            Response_OI_10(:,:,i) = TrialsResp_OI_10_2D;
            Response_OI_05(:,:,i) = TrialsResp_OI_05_2D;
        end
        
        SelectedResponses{1} = Response_OI_05;
        SelectedResponses{2} = Response_OI_10;
        
    case 1
        %% all units together, distringuishing for accel. values and starting TF
        clear Response_OI SelectedResponses
        for j = 1 : length(start_TF_values)
            for i = 1 : length(dTF_values)
                Accel_Con = dTF_2Dmask==dTF_values(i);
                Start_Con = StarTF_2dmask==start_TF_values(j);
                Stim_Con = Accel_Con & Start_Con;
                % look at trials with dTF lasting 1 second
                TrialsSelection_10 = Stim_Con & dTF_long_idx;
                TrialsSelection_05 = Stim_Con & ~dTF_long_idx;
                % generate 3D matrix for including time responses
                t1 = double(repmat(TrialsSelection_10, 1, 1, size(AM_UnitResponses_smooth,3))); % repeating the selection across time
                t2 = double(repmat(TrialsSelection_05, 1, 1, size(AM_UnitResponses_smooth,3))); % repeating the selection across time
                % nan all the zero values
                t1(t1(:)==0) = nan; t2(t2(:)==0) = nan;
                % select only the responses during selected trials
                TrialsResp_OI_10 = t1.*AM_UnitResponses_smooth;
                TrialsResp_OI_05 = t2.*AM_UnitResponses_smooth;
                % compute the mean responses for each unit during trials at specified
                % dTF
                TrialsResp_OI_10_2D = squeeze(nanmean(TrialsResp_OI_10,1));
                TrialsResp_OI_05_2D = squeeze(nanmean(TrialsResp_OI_05,1));
                % normalise the responses of interest
                Response_OI_10(:,:,i) = TrialsResp_OI_10_2D;
                Response_OI_05(:,:,i) = TrialsResp_OI_05_2D;
            end
            SelectedResponses{j,1} = Response_OI_05;
            SelectedResponses{j,2} = Response_OI_10;
        end
        
    case 2 
        %% all units separately, distringuishing for accel. values only
        clear Response_OI_10 Response_OI_05 SelectedResponses
        for i = 1 : length(dTF_values)
            Accel_Con = dTF_2Dmask==dTF_values(i);
            %Start_Con = StarTF_2dmask==start_TF_values(j);
            Stim_Con = Accel_Con; % & Start_Con;
            % look at trials with dTF lasting 1 second
            TrialsSelection_10 = Stim_Con & dTF_long_idx;
            TrialsSelection_05 = Stim_Con & ~dTF_long_idx;
            % generate 3D matrix for including time responses
            t1 = double(repmat(TrialsSelection_10, 1, 1, size(AM_UnitResponses_smooth,3))); % repeating the selection across time
            t2 = double(repmat(TrialsSelection_05, 1, 1, size(AM_UnitResponses_smooth,3))); % repeating the selection across time
            % nan all the zero values
            t1(t1(:)==0) = nan; t2(t2(:)==0) = nan;
            % select only the responses during selected trials
            TrialsResp_OI_10 = t1.*AM_UnitResponses_smooth;
            TrialsResp_OI_05 = t2.*AM_UnitResponses_smooth;
            SelectedResponses{i,1} = TrialsResp_OI_05;
            SelectedResponses{i,2} = TrialsResp_OI_10;
        end  
        
    case 3
        %% plot single units, distringuishing for accel. values and starting TF
        clear Response_OI
        for i = 1 : length(dTF_values)
            for j = 1 : length(start_TF_values)
                Accel_Con = dTF_2Dmask==dTF_values(i);
                Start_Con = StarTF_2dmask==start_TF_values(j);
                Stim_Con = Accel_Con & Start_Con;
                % generate 3D matrix for including time responses
                t1 = double(repmat(Accel_Con, 1, 1, size(AM_UnitResponses_smooth,3))); % repeating the selection across time
                % nan all the zero values
                t1(t1(:)==0) = nan;
                % select only the responses during selected trials
                TrialsResp_OI = t1.*AM_UnitResponses_smooth;
                % compute the mean responses for each unit during trials at specified
                % dTF
                TrialsResp_OI_2D = squeeze(nanmean(TrialsResp_OI,1));
                % normalise the responses of interest
                TrialsResp_max = nanmax(TrialsResp_OI_2D,[],2);
                Response_OI{j,i} = TrialsResp_OI_2D./TrialsResp_max;
            end
        end
end








end

