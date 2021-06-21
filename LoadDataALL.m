%% Tomaso Muzzu - UCL - 11/09/2019

function [ProjectData AM_UnitResponses AM_Param AM_Speed AM_UOI SelectedResponses AM_UnitResponses_smooth AM_EyeTracking] = LoadDataALL

%% Load raw and formatted data
if ~exist('ProjectData','var')
    [ProjectData AM_UnitResponses AM_Param AM_Speed AM_Time AM_UOI AM_EyeTracking] = LoadData;
end

%% Initial selection of best units
% to be run once in theory, repeated in case params need changing
FR_thres = 0.1;
AM_UOI = UnitsSelection(ProjectData,AM_UnitResponses,FR_thres);
% save([ DataFolder filesep 'UnitsResponse.mat'], 'AM_UnitResponses', 'AM_Speed', 'AM_Param', 'AM_UOI','-v7.3')
% ATM, units with mean FR>0.1Hz and comparable FR in first, middle, and last third
% of recording are kept

%% options for control experiments in naive mice
% rec_1day = [1 4 7];
% rec_2day = [2 5 8];
% rec_3day = [3 6 9];
% for i = 1:length(rec_1day)
%     units_1day(:,i) = (AM_Param(1,:,1) == rec_1day(i));
%     units_2day(:,i) = (AM_Param(1,:,1) == rec_2day(i));
%     units_3day(:,i) = (AM_Param(1,:,1) == rec_3day(i));
% end
% AM_UOI(:,2) = sum(units_1day,2);
% AM_UOI(:,3) = sum(units_2day,2);
% AM_UOI(:,4) = sum(units_3day,2);
% 
% AM_UOI(:,2) = ((AM_Param(1,:,1) >= 12));


%% smooth data with Gaussian filter
% 3rd argument is width of Gaussian filter
% 4th OPTIONAL argument is width of std - sigma=1/3 width by defeault;
Gauss_width = 0.3; % in seconds
AM_UnitResponses_smooth = Smooth_AM_FR(ProjectData,AM_UnitResponses,Gauss_width);

%% apply conditions to select responses
% conditions are set to define subset of units (AM_UOI), trials, periods within
% trials based on info store in matrix AM_Param(possible trials, units,
% conditions) where last dimension is 
% conditions = 1 --> nr of recording
% conditions = 2 --> grating direction [0:45:315]
% conditions = 3 --> 1 pert ON, 0 pert OFF
% conditions = 4 --> pert onset
% conditions = 5 --> pert offset

%Cond_Sel = {} ; % {optional('run'), optional(grating dir=[0:45:315])}
Units_Sel = AM_UOI(:,1); % this could include further selection criteria steming from other analysis
% Units_Sel = logical(ones(size(AM_UOI(:,1),1),1));
Run_TH = 2; % running speed threshold is 1 cm/s to be used like this: ..., 'run', Run_TH)
clear SelectedResponses

SelectedResponses = DataSelector(ProjectData,AM_UnitResponses_smooth,AM_Param,AM_Speed(:,:,1:end-1),Units_Sel);%,'run',Run_TH);
% SelectedResponses1 = DataSelector1(ProjectData,AM_UnitResponses_smooth,AM_Param,AM_Speed(:,:,1:end-1),Units_Sel);%,'run',Run_TH);
% SelectedResponses = DataSelector(ProjectData,AM_UnitResponses_smooth,AM_Param,AM_Speed(:,:,1:end-1),Units_Sel,'run',Run_TH);
% SelectedResponses = DataSelector(ProjectData,AM_UnitResponses_smooth,AM_Param,AM_Speed(:,:,1:end-1),Units_Sel,'run',Run_TH);

% % output contains a cell with matrices containing the trial mean responses 
% output contains a cell with matrices containing the trial mean responses 
% of units with their relevant "control" responses: e.g. if we want to look 
% at the responses during perturbation while the mouse runs for any grating 
% direction, SelectedResponses will be a cell with:
% SelectedResponses{1}=[units, trialtime] for the responses specified
% in above conditions and SelectedResponses{2}=[units, trialtime] during
% stationary periods.

%% plot the responses
% if iscell(SelectedResponses{1})
%     for i = 1:length(SelectedResponses{1})
%         clear SelResponse
%         SelResponse{1} = SelectedResponses{1}{i}; % response of interest
%         SelResponse{2} = SelectedResponses{2}{i}; % ctrl response of interest
%         PlotResponses(ProjectData,SelResponse,AM_Param,AM_Speed,Units_Sel,0);
%     end
% else

  

  % 25/05/2021 double input of SelectedResponses for verifying visual effect noted by
  % reviewer -> plot first input and use second input to order the units
  % try different colours
%    cmp = crameri('vik');
%   
%   PlotResponses(ProjectData,SelectedResponses,SelectedResponses,AM_Param,AM_Speed,Units_Sel,0)
%   colormap(cmp)
%   RedWhiteBlue
%   
%   figure
%   plot(cmp(:,1),'r')
%   hold on
%   plot(cmp(:,2),'g')
%   hold on
%   plot(cmp(:,3),'b')
%   tint_factor = 0.15
%   cmp(:,1) = cmp(:,1) + (1 - cmp(:,1)) * tint_factor
%   cmp(:,2) = cmp(:,2) + (1 - cmp(:,2)) * tint_factor
%   cmp(:,3) = cmp(:,3) + (1 - cmp(:,3)) * tint_factor
%   plot(cmp(:,1),'r--')
%   hold on
%   plot(cmp(:,2),'g--')
%   hold on
%   plot(cmp(:,3),'b--')
  
%   
%   ColormapNames = {'parula', 'hsv', 'hot', 'cool', 'spring', ...
%                    'summer', 'autumn', 'winter', 'gray', 'bone', 'copper', 'pink'};
%   CrameriCM =  {'broc', 'cork', 'vik', 'lisbon', 'tofino',  'berlin', 'roma', 'bam', 'vanimo'};
%   CrameriCM =  {'acton','bam','bamO','bamako','batlow','batlowK','batlowW','berlin','bilbao','broc','brocO','buda','bukavu','cork',...
%                 'corkO','davos','devon','fes','grayC','hawaii','imola','lajolla','lapaz','lisbon',...
%                 'nuuk','oleron','oslo','roma','romaO','tofino','tokyo','turku','vanimo','vik','vikO'}          
%   figure
%   for i = 1:size(CrameriCM,2)        
%       ax(i) = subplot(2,6,i)
%       PlotResponsesColormap(ProjectData,SelectedResponses1,SelectedResponses,AM_Param,AM_Speed,Units_Sel,0);
%       
%       cmp = crameri(['-' CrameriCM{i}]);
% %       cmp = ColormapNames{i};
%       
%         colormap(ax(i),cmp)
%       
%       title(['-' CrameriCM{i}]);
% %       title(ColormapNames{i})
%   end
  
  
%  PlotResponses(ProjectData,SelectedResponses,AM_Param,AM_Speed,Units_Sel,0);

  %PlotResponsesRunStill(ProjectData,{SelectedResponses{1,1}(PertRespUnits_pos,:) SelectedResponses{1,2}(PertRespUnits_pos,:)} ,AM_Param,AM_Speed,Units_Sel,0,1)
% end
        
end


