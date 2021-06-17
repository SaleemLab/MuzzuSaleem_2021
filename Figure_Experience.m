%% Tomaso Muzzu - UCL - 23/03/2020 - Experience effects

%% Here :
% 1) responses in first trials and last trials - general
% 2) responses in first trials and last trials - agnel 

%% load data
if ~exist('ProjectData','var')
    [ProjectData AM_UnitResponses AM_Param AM_Speed AM_UOI SelectedResponses AM_UnitResponses_smooth] = LoadDataALL;
end

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

%% first 7 animals
if  size(ProjectData,1)>10
    CTRL_exp = 0;
    Animal_1st_idx = [1 5 7 12 15 24 31];
    
    if ~exist('PertResp_units','var')
        % select only perturbation responsive units
        load('AUC_shuffled.mat')
        Sh_responses = AUC_shuffled(:,2:end);
        p_pert_th = prctile(Sh_responses(:),99);
        PertResp_units = (AUC_shuffled(:,1)>p_pert_th);
        % select only pos. modulated perturbation responsive units
        load('DM_pert_shuffled.mat')
        DM = DM_sh(:,1);
        DM_sign_i(:,1) = DM>0;
        DM_sign_i(:,2) = DM<=0;
        % select only pos. modulated perturbation responsive units
        PertRespUnits_pos = PertResp_units & DM_sign_i(:,1);
        PertRespUnits_neg = PertResp_units & DM_sign_i(:,2);
    end
    
else
    CTRL_exp = 1;
    % naive animals
    Animal_1st_idx = [1 4 7];
end

%% direction
if size(AM_UOI,2)==1
    SelectedCells = AM_UOI;
else
    SelectedCells = AM_UOI(:,1) & AM_UOI(:,2);
end
% select trials of interest and control trials as well
Trials_PertON = AM_Param(:,:,3)==1; % find indexes where perturbation is on
Trials_PertOFF= AM_Param(:,:,3)==0; % find indexes where perturbation is off

BonvisionFR = 60; %Hz
trialSide_samples = 60;
trialSide_seconds = 1;

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
param_sel = AM_Param(:,SelectedCells,:);
UnitResponses_smooth = AM_UnitResponses_smooth(:,SelectedCells,:);

%% DEFINE INDEXES FOR FINDING INDIVUAL SESSIONS
for i = 1:size(ProjectData,1)
    UnitCount(i) = size(ProjectData.Units_Info{i},1);
end

%% save the perturbation responses in early, middle and late stages
% load results from significance 
if CTRL_exp == 1
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
end

%% look at mean responses 
% FIRST RECORDING SESSION PER ANIMAL
% 1) look at perturbation responses of pert. resp. units
% 2) look at perturbation responses of all units
% 3) look at vis. stim. responses of pert. resp. units
% 4) look at vis. stim. responses of all units
% 5) differentiate by the angle?
% REPEAT THIS ANALYSIS FOR THE OTHER SESSIONS

%% 1)
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
% Sessions_Selection = unique(PertUnits_es_idx);
param_sel = param_sel(:,PertUnits_idx,:);
k = 1; clear MIs
plot_SU = 0;
for i = 1:length(Sessions_Selection)
    Units_OI = find(PertUnits_es_idx==Sessions_Selection(i));
    if ~isempty(Units_OI)
        for j = 1:length(Units_OI)
            PertResp_oneUnit = squeeze(TrialsResp_OI(squeeze(~isnan(TrialsResp_OI(:,Units_OI(j),1))),Units_OI(j),:))*60;
            Angle_Info = param_sel(squeeze(~isnan(TrialsResp_OI(:,Units_OI(j),1))),Units_OI(j),2);
            % first 10, last 10
            % first third, last third
            % total average
            PertOnsets = param_sel(squeeze(~isnan(TrialsResp_OI(:,Units_OI(j),1))),Units_OI(j),4); % 2D matrix
            PertOffsets = param_sel(squeeze(~isnan(TrialsResp_OI(:,Units_OI(j),1))),Units_OI(j),5); % 2D matrix
            trials_nr = 10;
            Directions = unique(Angle_Info);
            % Normalise each angle response by its max
            clear MaxDirectResp PertResp_oneUnit_norm
            for m = 1:length(Directions)
                MaxDirectResp(m) = max(mean(PertResp_oneUnit(Angle_Info==Directions(m),:))); 
                PertResp_oneUnit_norm(Angle_Info==Directions(m),:) = PertResp_oneUnit(Angle_Info==Directions(m),:)/MaxDirectResp(m);
            end
            %             PertResp_oneUnit_norm = PertResp_oneUnit;
            MeanCurve = mean(PertResp_oneUnit_norm(1:trials_nr,:));
            MIs(k,1) = (max(MeanCurve(po_i(1):po_i(1)+trialSide_samples))-max(MeanCurve(po_i(1)-trialSide_samples:po_i(1))))/...
                (max(MeanCurve(po_i(1):po_i(1)+trialSide_samples))+max(MeanCurve(po_i(1)-trialSide_samples:po_i(1))));
            
            MeanCurve = mean(PertResp_oneUnit_norm(end-trials_nr:end,:));
            MIs(k,2) = (max(MeanCurve(po_i(1):po_i(1)+trialSide_samples))-max(MeanCurve(po_i(1)-trialSide_samples:po_i(1))))/...
                (max(MeanCurve(po_i(1):po_i(1)+trialSide_samples))+max(MeanCurve(po_i(1)-trialSide_samples:po_i(1))));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plot_SU
                figure
                h1 = shadedErrorBar(TimeLine,mean(PertResp_oneUnit_norm(1:trials_nr,:))/max(max(PertResp_oneUnit_norm(1:trials_nr,:))),...
                    std(PertResp_oneUnit_norm(1:trials_nr,:)/max(max(PertResp_oneUnit_norm(1:trials_nr,:))))/sqrt(trials_nr),...
                    'lineprops',{'Color',[255,127,80]/255,'markerfacecolor',[255,127,80]/255});
                hold on
                h2 = shadedErrorBar(TimeLine,mean(PertResp_oneUnit_norm(end-trials_nr:end,:))/max(max(PertResp_oneUnit_norm(end-trials_nr:end,:))), ...
                    std(PertResp_oneUnit_norm(end-trials_nr:end,:)/max(max(PertResp_oneUnit_norm(end-trials_nr:end,:))))/sqrt(trials_nr),...
                    'lineprops',{'Color',[1 0 0],'markerfacecolor',[1 0 0]});
                %             figure
                %             h1 = shadedErrorBar(TimeLine,mean(PertResp_oneUnit(1:trials_nr,:))/max(max(PertResp_oneUnit(1:trials_nr,:))),...
                %                                 std(PertResp_oneUnit(1:trials_nr,:)/max(max(PertResp_oneUnit(1:trials_nr,:))))/sqrt(trials_nr),...
                %                                 'lineprops',{'Color',[255,127,80]/255,'markerfacecolor',[255,127,80]/255});
                %             hold on
                %             h2 = shadedErrorBar(TimeLine,mean(PertResp_oneUnit(end-trials_nr:end,:))/max(max(PertResp_oneUnit(end-trials_nr:end,:))), ...
                %                                 std(PertResp_oneUnit(end-trials_nr:end,:)/max(max(PertResp_oneUnit(end-trials_nr:end,:))))/sqrt(trials_nr),...
                %                                 'lineprops',{'Color',[1 0 0],'markerfacecolor',[1 0 0]});
                hold on
                MaxValue = max(max(mean(PertResp_oneUnit_norm(1:trials_nr,:))/max(max(PertResp_oneUnit_norm(1:trials_nr,:)))),...
                    max(mean(PertResp_oneUnit_norm(end-trials_nr:end,:))/max(max(PertResp_oneUnit_norm(end-trials_nr:end,:)))));
                plot([0 0],[0 MaxValue],'-k','LineWidth',0.75);
                hold on
                plot([TrialEnd TrialEnd], [0 MaxValue],'-k','LineWidth',0.75);
                hold on
                plot([nanmean(PertOnsets(1:trials_nr)), mean(PertOnsets(end-trials_nr:end))], [0 MaxValue],'-b','LineWidth',1);
                hold on
                plot([nanmean(PertOffsets(1:trials_nr)), mean(PertOffsets(end-trials_nr:end))], [0 MaxValue],'-b','LineWidth',1);
                xlabel('Seconds'); ylabel('Spikes/s');
                set(gca,'TickDir','out'); xlim([-0.5 max(TimeLine)])
                ll = legend({'First 10 trials', 'Last 10 trials'},'fontsize',10);
                ll.Color = 'none'; ll.EdgeColor = 'none';
                title([ProjectData.Mouse_name{Sessions_Selection(i)} ...
                    ' unit='  num2str(PertUnits_idx(Units_OI(j))) ...
                    ' AUC='   num2str(AUC_real(PertUnits_idx(Units_OI(j)))) ...
                    ' MIs='   num2str(MIs(k,1)) ',' num2str(MIs(k,2)) ...
                    ]);
                saveas(gcf,['E:\Dropbox\UCL\Tomaso-Aman\VisPerturbation\Fig_3\panels\Naive\examples' filesep 'PertUnits_' num2str(k) '.pdf' ])
                close all
            end
            k = k+1;
        end
    end
end

%% plot change of MIs between first and last 10 trials
% units to exclude 12 and 19
figure
n=0;
for i = 1:length(MIs)
    if CTRL_exp == 1
        if i~=12 & i~=19
            plot([1 2],[MIs(i,1) MIs(i,2)],'k')
            hold on
            n=n+1;
        end
    else
        plot([1 2],[MIs(i,1) MIs(i,2)],'k')
        hold on
        n=n+1;
    end
end
hold on
boxplot(MIs,'notch','on')
set(gca,'XTicklabel',{'first 10 trials', 'last 10 trials'},'box','off','TickDir','out')
ylabel('MI of pert. resp.')
title(['n = ' num2str(n)])

figure
plot(MIs(:,1),MIs(:,2),'or')
xlabel('MI - first 10 trials')
ylabel('MI - last 10 trials')
set(gca,'box','off','TickDir','out','FontSize',13)
xlim([0 1]); ylim([0 1])
hold on
plot([0 1],[0 1],':')

[p h] = kstest(MIs(:,1)-MIs(:,2))

p = kruskalwallis(MIs)





%% plot session by session

if CTRL_exp == 1
    % naive animals 
    Rec_nr = [1 2 3 1 2 3 1 2 3 4];
else
    % standad
    Rec_nr = [1 2 3 4 ...
              1 2 ...
              1 2 3 4 5 ...
              1 2 3 ...
              1 2 3 4 5 6 7 8 9 ...
              1 2 3 4 5 6 7 ...
              1 2 3 4 5 6 7];
end

figure
for s=1:length(Sessions_Selection)
    % plot unit by unit
    clear MIs_
    MIs_ = MIs(find(PertUnits_es_idx==Sessions_Selection(s)),:);
    if CTRL_exp == 1
        subplot(1,length(Sessions_Selection),s)
    else
        subplot(2,ceil(length(Sessions_Selection)/2),s)
    end
    n=0;
    for i = 1:size(MIs_,1)
        if i~=12 & i~=19
            plot([1 2],[MIs_(i,1) MIs_(i,2)],'k')
            hold on
            n=n+1;
        end
    end
%     hold on
%     boxplot(MIs_,'notch','on')
    set(gca,'XTick', [1 2], 'XTicklabel',{'first 10 trials', 'last 10 trials'},'box','off','TickDir','out')
    ylabel('MI of pert. resp.')
    title(['n = ' num2str(n)])
    title([ProjectData.Mouse_name{Sessions_Selection(s)} ...
                    ' rec.nr = ' num2str(Rec_nr(Sessions_Selection(s))) ...
                    ', n = ' num2str(n) ...   
                    ]);
    xlim([0 3]); ylim([0 0.7])
    
end



%%

RecunitCounting = unique(param_sel(1,:,1))
for i = 1:length(RecunitCounting)
    UnitCounting(i) = sum(param_sel(1,:,1)==RecunitCounting(i))
end


