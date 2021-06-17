%% Tomaso Muzzu - UCL - 31/01/2020

% experience related effects on perturbatio responses

% look at session by session and compare responses between 1of3 2of3 3of3
% periods
% look at all units and only perturbation responsive units
% see if there are any effects on the vis. stim. responses

%%
trialSide_samples = 60;
trialSide_seconds = 1;

PertOnsets = AM_Param(:,:,4); % 2D matrix
PertOffsets = AM_Param(:,:,5); % 2D matrix
[v min_i(1)] = max(PertOnsets(:)); % find the latest moment when pert onset happens
[v min_i(2)] = max(PertOffsets(:)); % find the latest moment when pert offset happens

ReferSize = [size(AM_Param,1),size(AM_Param,2)];
[trial_el(1), trial_el(2)] = ind2sub(ReferSize,min_i(1)); % 2D index of minimum pert onset time

Rec = AM_Param(trial_el(1), trial_el(2),1);
% define timeline of example recording
TrialStart = ProjectData.Session_data{Rec,1}.ACInfo{1,1}.trialStartsEnds(trial_el(1),1);
TrialEnd = ProjectData.Session_data{Rec,1}.ACInfo{1,1}.trialStartsEnds(trial_el(1),2)-TrialStart;
TrialStart = 0;
TimeLine = linspace(TrialStart-trialSide_seconds, ...
                    TrialEnd+trialSide_seconds,...
                    size(AM_UnitResponses_smooth,3)); % -1 seconds
[v po_i(1)] = min(abs(TimeLine-min(PertOnsets(:))));
[v po_i(2)] = min(abs(TimeLine-min(PertOffsets(:))));
SelectedCells = AM_UOI;
UnitResponses_smooth = AM_UnitResponses_smooth(:,SelectedCells,:);

%% DEFINE INDEXES FOR FINDING INDIVUAL SESSIONS
for i = 1:size(ProjectData,1)
    UnitCount(i) = size(ProjectData.Units_Info{i},1);
end
Animal_1st_idx = [1 4 7];

%% save the perturbation responses in early, middle and late stages
% load results from significance 
load('AUC_shuffled_CTRL.mat');
AUC_real = AUC_shuffled(:,1);
AUC_thres = prctile(reshape(AUC_shuffled(:,2:end),size(AUC_shuffled(:,2:end),1)*size(AUC_shuffled(:,2:end),2),1) ,99);
PertRespUnits = AUC_real>AUC_thres;
% select only the positively modulated units
SelectedCells = AM_UOI;
param_sel = AM_Param(:,SelectedCells,:);
Param = squeeze(param_sel(:,:,3));
% Depth of modulation
r= 1; clear DM DepthMod
for i = 1:size(Param,2)
    % select trial responses
    Responses = Param(:,i);
    % scrub responses
    Responses = logical(Responses(~isnan(Responses)));
    UnitResponse = mean(squeeze(UnitResponses_smooth(Responses==1,i,:)))*60;
    Dir = sign(sum(UnitResponse(po_i(1):po_i(1)+trialSide_samples))-sum(UnitResponse(po_i(1)-trialSide_samples:po_i(1))));
    switch Dir
        case 1
            DM(i) = (max(UnitResponse(po_i(1):po_i(1)+trialSide_samples))-max(UnitResponse(po_i(1)-trialSide_samples:po_i(1))))/...
                (max(UnitResponse(po_i(1):po_i(1)+trialSide_samples))+max(UnitResponse(po_i(1)-trialSide_samples:po_i(1))));
        case -1
            DM(i) = (min(UnitResponse(po_i(1):po_i(1)+trialSide_samples))-min(UnitResponse(po_i(1)-trialSide_samples:po_i(1))))/...
                (min(UnitResponse(po_i(1):po_i(1)+trialSide_samples))+min(UnitResponse(po_i(1)-trialSide_samples:po_i(1))));
    end
end
PertRespUnits_pos = PertRespUnits==1 & DM'>0; % sum(PertRespUnits_pos)
PertRespUnits_neg = PertRespUnits==1 & DM'<=0; % sum(PertRespUnits_neg)

%% look at mean responses 
% FIRST RECORDING SESSION PER ANIMAL
% 1) look at perturbation responses of pert. resp. units
% 2) look at perturbation responses of all units
% 3) look at vis. stim. responses of pert. resp. units
% 4) look at vis. stim. responses of all units
% 5) differentiate by the angle?
% REPEAT THIS ANALYSIS FOR THE OTHER SESSIONS

% 1)
PertUnits_idx = find(PertRespUnits_pos);
PertUnits_es_idx = param_sel(1,PertUnits_idx,1);
Trials_PertON = param_sel(:,PertRespUnits_pos,3)==1; % find indexes where perturbation is on
Trials_PertOFF= param_sel(:,PertRespUnits_pos,3)==0; % find indexes where perturbation is off
% generate 3D matrix for including time responses
t1 = double(repmat(Trials_PertON, 1, 1, size(AM_UnitResponses_smooth,3))); % repeating the selection across time
t2 =  double(repmat(Trials_PertOFF, 1, 1, size(AM_UnitResponses_smooth,3))); % repeating the selection across time
% nan all the zero values
t1(t1(:)==0) = nan; t2(t2(:)==0) = nan;
% select only the responses during selected trials (pert ON and OFF separately)
TrialsResp_OI = t1.*UnitResponses_smooth(:,PertUnits_idx,:);
% look at responses now
Sessions_Selection = Animal_1st_idx;
Sessions_Selection = unique(PertUnits_es_idx);
param_sel = param_sel(:,PertUnits_idx,:);
k =1;
for i = 1:length(Sessions_Selection)
    Units_OI = find(PertUnits_es_idx==Sessions_Selection(i));
    if ~isempty(Units_OI)
        for j = 1:length(Units_OI)
            PertResp_oneUnit = squeeze(TrialsResp_OI(squeeze(~isnan(TrialsResp_OI(:,Units_OI(j),1))),Units_OI(j),:))*60;
            % first 10, last 10
            % first third, last third
            % total average
            PertOnsets = param_sel(squeeze(~isnan(TrialsResp_OI(:,Units_OI(j),1))),Units_OI(j),4); % 2D matrix
            PertOffsets = param_sel(squeeze(~isnan(TrialsResp_OI(:,Units_OI(j),1))),Units_OI(j),5); % 2D matrix
            trials_nr = 10;
            figure
            h1 = shadedErrorBar(TimeLine,mean(PertResp_oneUnit(1:trials_nr,:)), std(PertResp_oneUnit(1:trials_nr,:))/sqrt(trials_nr),...
                'lineprops',{'Color',[255,127,80]/255,'markerfacecolor',[255,127,80]/255});
            hold on
            h2 = shadedErrorBar(TimeLine,mean(PertResp_oneUnit(end-trials_nr:end,:)), std(PertResp_oneUnit(end-trials_nr:end,:))/sqrt(trials_nr),...
                'lineprops',{'Color',[1 0 0],'markerfacecolor',[1 0 0]});
            hold on
            MaxValue = max(max(mean(PertResp_oneUnit(end-trials_nr:end,:))),max(mean(PertResp_oneUnit(1:trials_nr,:))));
            plot([0 0],[0 MaxValue],'-k','LineWidth',0.75);
            hold on
            plot([TrialEnd TrialEnd], [0 MaxValue],'-k','LineWidth',0.75);
            hold on
            plot([nanmean(PertOnsets(1:trials_nr)), mean(PertOnsets(end-trials_nr:end))], [0 MaxValue],'-b','LineWidth',1);
            hold on
            plot([nanmean(PertOffsets(1:trials_nr)), mean(PertOffsets(end-trials_nr:end))], [0 MaxValue],'-b','LineWidth',1);
            
            MeanCurve = mean(PertResp_oneUnit(1:trials_nr,:));
            MIs(k,1) = (max(MeanCurve(po_i(1):po_i(1)+trialSide_samples))-max(MeanCurve(po_i(1)-trialSide_samples:po_i(1))))/...
                (max(MeanCurve(po_i(1):po_i(1)+trialSide_samples))+max(MeanCurve(po_i(1)-trialSide_samples:po_i(1))));
            
            MeanCurve = mean(PertResp_oneUnit(end-trials_nr:end,:));
            MIs(k,2) = (max(MeanCurve(po_i(1):po_i(1)+trialSide_samples))-max(MeanCurve(po_i(1)-trialSide_samples:po_i(1))))/...
                (max(MeanCurve(po_i(1):po_i(1)+trialSide_samples))+max(MeanCurve(po_i(1)-trialSide_samples:po_i(1))));
            xlabel('Seconds'); ylabel('Spikes/s');
            set(gca,'TickDir','out'); xlim([-0.5 max(TimeLine)])
            ll = legend({'First 10 trials', 'Last 10 trials'},'fontsize',10);
            ll.Color = 'none'; ll.EdgeColor = 'none';
            title([ProjectData.Mouse_name{Sessions_Selection(i)} ...
                   ' unit='  num2str(PertUnits_idx(Units_OI(j))) ...
                   ' AUC='   num2str(AUC_real(PertUnits_idx(Units_OI(j)))) ...
                   ' MIs='   num2str(MIs(k,1)) ',' num2str(MIs(k,2)) ...
                   ]);
            saveas(gcf,['D:\Dropbox\UCL\Tomaso-Aman\VisPerturbation\Fig_3\panels\Naive\examples' filesep 'PertUnits_' num2str(k) '.pdf' ])
            close all
            k = k+1;
        end
    end
end

%% plot MIs for before and after...
figure
UnivarScatter(MIs(:,2)-MIs(:,1),'Label',{'Perturbation units'},'MarkerFaceColor',[ 0.5 0.5 0.5])
ylabel('MI_e_n_d - MI_s_t_a_r_t')


% how to compare responses: 
% modulation index difference for single units
% population responses are not necessary






