%% Tomaso Muzzu - UCL - 12/09/2019

% initial selection on best and/or desired units


function SelectedResponses = DataSelector1(ProjectData,AM_UnitResponses_smooth,AM_Param,AM_Speed,AM_UOI,varargin)

% look at perturbation responses regardless of running or grating
% direction
if size(AM_UOI,2)==1
    SelectedCells = AM_UOI;
else
    SelectedCells = AM_UOI(:,1) & AM_UOI(:,2);
end
% select trials of interest and control trials as well
Trials_PertON = AM_Param(:,:,3)==1; % find indexes where perturbation is on
Trials_PertOFF= AM_Param(:,:,3)==0; % find indexes where perturbation is off
% Trials_PertOFF_lim = logical(zeros(size(Trials_PertOFF,1),size(Trials_PertOFF,2)));
% % select an equal number of trials for the off conditions.
% for i = 1:size(Trials_PertOFF,2)
%     PertOFF_idx = find(Trials_PertOFF(:,i));
%     PertOFF_idx_r = randperm(length(PertOFF_idx));
%     Trials_PertOFF_lim(PertOFF_idx_r(1:sum(Trials_PertON(:,1))),i) = 1;
% end
% Trials_PertOFF = Trials_PertOFF_lim;

if isempty(varargin)
    
    % generate 3D matrix for including time responses
    t1 = double(repmat(Trials_PertON, 1, 1, size(AM_UnitResponses_smooth,3))); % repeating the selection across time
    t2 =  double(repmat(Trials_PertOFF, 1, 1, size(AM_UnitResponses_smooth,3))); % repeating the selection across time
    % nan all the zero values
    t1(t1(:)==0) = nan; t2(t2(:)==0) = nan;
    
    % select only the responses during selected trials (pert ON and OFF separately)
    TrialsResp_OI = t1.*AM_UnitResponses_smooth;
    TrialsRespControl_OI = t2.*AM_UnitResponses_smooth;
    % compute the mean responses for each unit during trials with pert ON
    % and OFF
    TrialsResp_OI_2D = squeeze(nanmean(TrialsResp_OI,1));
    TrialsRespControl_OI_2D = squeeze(nanmean(TrialsRespControl_OI,1));
    % normalise the responses of interest
    TrialsResp_max = max(nanmax(TrialsResp_OI_2D,[],2),nanmax(TrialsRespControl_OI_2D,[],2));
%     Response_OI = TrialsResp_OI_2D./TrialsResp_max;
%     ResponseControl_OI = TrialsRespControl_OI_2D./TrialsResp_max;   
    TrialsResp_min = min(nanmin(TrialsResp_OI_2D,[],2),nanmin(TrialsRespControl_OI_2D,[],2));
    Response_OI = (TrialsResp_OI_2D - TrialsResp_min) ./ (TrialsResp_max - TrialsResp_min);
    ResponseControl_OI = (TrialsRespControl_OI_2D - TrialsResp_min) ./(TrialsResp_max - TrialsResp_min);   

else
    if isnumeric(varargin{1}) % look at perturbation responses regardless of running for specific grating directions
        Angles_OI = varargin{1};
        [AM_AngleResponse AM_AngleResponseControl] = ComputeAngleResponse(AM_Param, AM_Speed,SelectedCells,Trials_PertON,Trials_PertOFF,ProjectData, Angles_OI);
        for i = 1:size(AM_AngleResponse,3)
            % generate 3D matrix for including time responses
            t1 = double(repmat(AM_AngleResponse(:,:,i), 1, 1, size(AM_UnitResponses_smooth,3))); % repeating the selection across time
            t2 =  double(repmat(AM_AngleResponseControl(:,:,i), 1, 1, size(AM_UnitResponses_smooth,3))); % repeating the selection across time
            % nan all the zero values
            t1(t1(:)==0) = nan; t2(t2(:)==0) = nan;
            % select only the responses during selected trials (pert ON and OFF separately)
            TrialsResp_OI = t1.*AM_UnitResponses_smooth;
            TrialsRespControl_OI = t2.*AM_UnitResponses_smooth;
            % compute the mean responses for each unit during trials with pert ON
            % and OFF
            TrialsResp_OI_2D = squeeze(nanmean(TrialsResp_OI,1));
            TrialsRespControl_OI_2D = squeeze(nanmean(TrialsRespControl_OI,1));
            TrialsResp_max = max(nanmax(TrialsResp_OI_2D,[],2),nanmax(TrialsRespControl_OI_2D,[],2));
%             Response = TrialsResp_OI_2D./TrialsResp_max;
%             ResponseControl = TrialsRespControl_OI_2D./TrialsResp_max;
            TrialsResp_min = min(nanmin(TrialsResp_OI_2D,[],2),nanmin(TrialsRespControl_OI_2D,[],2));
            Response = (TrialsResp_OI_2D - TrialsResp_min) ./ (TrialsResp_max - TrialsResp_min);
            ResponseControl = (TrialsRespControl_OI_2D - TrialsRsep_min) ./(TrialsResp_max - TrialsResp_min);   
            Response_OI{i} = Response(SelectedCells,:);
            ResponseControl_OI{i} = ResponseControl(SelectedCells,:);
        end
        
    else % look at perturbation responses during running regardless of grating directions
        Run_TH = varargin{2};
        [AM_SpeedResponse AM_SpeedResponseControl] = ComputeSpeedResponse(AM_Param, AM_Speed,SelectedCells,Trials_PertON,Trials_PertOFF,ProjectData);
        % AM_SpeedResponse(:,:,1) = mean speed -0.5s to +2s at visual stimulus onset
        % AM_SpeedResponse(:,:,2) = mean speed -0.5s to end of perturbation period
        % AM_SpeedResponse(:,:,3) = mean speed -0.5s to +1s at perturbation offset
        % AM_SpeedResponse(:,:,4) = mean speed across entire trial
        minTrials = 4;
        [SelectedTrials_OI  SelectedTrials_OI_control TotTrials] = SpeedFiltering(AM_Speed, AM_SpeedResponse, AM_SpeedResponseControl, Trials_PertON,Trials_PertOFF, Run_TH,minTrials);
        % SelectedTrials_OI{4,2} = 4 rows as above, column 1 = pert ON & run; column 2 = pert ON & static
        % SelectedTrials_OI_control{4,2} = 4 rows as above, column 1 = pert OFF & run; column 2 = pert OFF & static
        
        % generate 3D matrix for including time responses
        t1 = double(repmat(SelectedTrials_OI{2,1}, 1, 1, size(AM_UnitResponses_smooth,3))); % repeating the selection across time
        t2 =  double(repmat(SelectedTrials_OI{2,2}, 1, 1, size(AM_UnitResponses_smooth,3))); % repeating the selection across time
        % nan all the zero values
        t1(t1(:)==0) = nan; t2(t2(:)==0) = nan;
        
        % select only the responses during selected trials (pert ON and OFF separately)
        TrialsResp_OI = t1.*AM_UnitResponses_smooth;
        TrialsRespControl_OI = t2.*AM_UnitResponses_smooth;
        % compute the mean responses for each unit during trials with pert ON
        % and OFF
        TrialsResp_OI_2D = squeeze(nanmean(TrialsResp_OI,1));
        TrialsRespControl_OI_2D = squeeze(nanmean(TrialsRespControl_OI,1));
        % normalise the responses of interest
        TrialsResp_max = max(nanmax(TrialsResp_OI_2D,[],2),nanmax(TrialsRespControl_OI_2D,[],2));
%         Response_OI = TrialsResp_OI_2D./TrialsResp_max;
%         ResponseControl_OI = TrialsRespControl_OI_2D./TrialsResp_max;
        TrialsResp_min = min(nanmin(TrialsResp_OI_2D,[],2),nanmin(TrialsRespControl_OI_2D,[],2));
        Response_OI = (TrialsResp_OI_2D - TrialsResp_min) ./ (TrialsResp_max - TrialsResp_min);
        ResponseControl_OI = (TrialsRespControl_OI_2D - TrialsResp_min) ./(TrialsResp_max - TrialsResp_min);
        
        % desired function
        % function trialSpeedFilter(AM_Speed, start_time, stop_time, thres_speed, 'type')
        % % type can be meanSpeed / maxSpeed / median / adaptive
        % end
    end
    
end

% select only responses of units of interest
clear SelectedResponses
if iscell(Response_OI)
    SelectedResponses{1} = Response_OI;
    SelectedResponses{2} = ResponseControl_OI;
else
    SelectedResponses{1} = Response_OI(SelectedCells,:);
    SelectedResponses{2} = ResponseControl_OI(SelectedCells,:);
end

end


%% Template with for loop (obselete from 13/09/2019)

% find perturbation and non-perturbation trials
%     SelectedTrials{1} = find(AM_Param(:,:,3)==1); % find indexes where perturbation is on
%     SelectedTrials{2} = find(AM_Param(:,:,3)==0); % find indexes where perturbation is off
%     
%     ReferSize = [size(AM_Param,1),size(AM_Param,2)];
%     [I_el(:,1),I_el(:,2)] = ind2sub(ReferSize,SelectedTrials{1}); % from linear to 2D indexes
%     [I_el_ctrl(:,1),I_el_ctrl(:,2)] = ind2sub(ReferSize,SelectedTrials{2}); % from linear to 2D indexes


% % find total nr of recordings to scan units in batches
% nr_recordings = unique(AM_Param(1,:,1));
% % initialise arrays for the responses of interest
% Response_OI = ones(1,size(AM_UnitResponses,3));
% ResponseControl_OI = ones(1,size(AM_UnitResponses,3));
% 
% clear Response_OI_temp ResponseControl_OI_temp
% for i = 1:length(nr_recordings)
%     Response_OI_temp = squeeze(mean(...
%         AM_UnitResponses_smooth(I_el(I_el(:,2)==max(find(AM_Param(1,:,1)==i)),1), ...
%         find(AM_Param(1,:,1)==i),...
%         :)...
%         ,1));
%     ResponseControl_OI_temp = squeeze(mean(...
%         AM_UnitResponses_smooth(I_el_ctrl(I_el_ctrl(:,2)==max(find(AM_Param(1,:,1)==i)),1), ...
%         find(AM_Param(1,:,1)==i),...
%         :)...
%         ,1));
%     % find max responses for each unit
%     maxresponse = max(max(Response_OI_temp,[],2), max(ResponseControl_OI_temp,[],2));
%     % normalise the unit responses and control responses by the same value
%     Response_OI_temp = Response_OI_temp./maxresponse ;
%     ResponseControl_OI_temp = ResponseControl_OI_temp./maxresponse;
%     % concatenate the reponses in a single matrix
%     Response_OI = [Response_OI; Response_OI_temp];
%     ResponseControl_OI = [ResponseControl_OI; ResponseControl_OI_temp];
% end
% Response_OI = Response_OI(2:end,:);
% ResponseControl_OI = ResponseControl_OI(2:end,:);
% 
% % select only responses of units of interest
% SelectedResponses{1} = Response_OI(SelectedCells,:);
% SelectedResponses{2} = ResponseControl_OI(SelectedCells,:);


