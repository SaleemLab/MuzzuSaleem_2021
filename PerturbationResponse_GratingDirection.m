%% Tomaso Muzzu - UCL - 23/10/2019

% compute pert responses in different gratinf direciton conditions

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

BonvisionFR = 60; %Hz
trialSide_samples = 60;
trialSide_seconds = 1;

%% with components congruent with running direction
AngleCombo = [0 ; 180]; 
AngleCombo = [45 ; 225];
AngleCombo = [90 ; 270];
AngleCombo = [135 ; 315];

AngleCombo = [0 180; 90 270];
AngleCombo = [0 45 315; 135 180 225];

clear SelectedResponses
for a = 1:size(AngleCombo,1)
    Angles_OI = AngleCombo(a,:);
    clear AngleFilter
    for i = 1:length(Angles_OI)
        % find all PertON trials with the grating direction of interest
        AngleFilter_t(:,:,i) =  AM_Param(:,:,2)==Angles_OI(i);
        if i == length(Angles_OI)
            AngleFilter = sum(AngleFilter_t,3);
        end
    end
    clear AM_AngleResponse AM_AngleResponseControl
    AM_AngleResponse = AngleFilter & Trials_PertON;
    % find all PertOFF trials with the grating direction of interest
    AM_AngleResponseControl = AngleFilter & Trials_PertOFF;  
    
    
    % generate 3D matrix for including time responses
    t1 = double(repmat(AM_AngleResponse, 1, 1, size(AM_UnitResponses_smooth,3))); % repeating the selection across time
    t2 =  double(repmat(AM_AngleResponseControl, 1, 1, size(AM_UnitResponses_smooth,3))); % repeating the selection across time
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
    Response = TrialsResp_OI_2D./TrialsResp_max;
    ResponseControl = TrialsRespControl_OI_2D./TrialsResp_max;
    Response_OI = Response(SelectedCells,:);
    ResponseControl_OI = ResponseControl(SelectedCells,:);
    
    SelectedResponses{a,1} = Response_OI;
    SelectedResponses{a,2} = ResponseControl_OI;
    
end

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
    UnitsPertPos = PertResp_units & DM_sign_i(:,1);
    UnitsPertNeg = PertResp_units & DM_sign_i(:,2);
end
%PlotResponses(ProjectData,{SelectedResponses{1,1}(UnitsPertPos,:) SelectedResponses{2,1}(UnitsPertPos,:)},AM_Param,AM_Speed,UnitsPertPos,0); 
%PlotResponses(ProjectData,{SelectedResponses{1,1}(UnitsPertNeg,:) SelectedResponses{2,1}(UnitsPertNeg,:)},AM_Param,AM_Speed,UnitsPertNeg,0); 

PlotMeanCurveAngle(ProjectData,{SelectedResponses{1,1}(UnitsPertPos,:) SelectedResponses{2,1}(UnitsPertPos,:)},AM_Param,AM_Speed,UnitsPertPos,AngleCombo)
PlotMeanCurveAngle(ProjectData,{SelectedResponses{1,1}(UnitsPertNeg,:) SelectedResponses{2,1}(UnitsPertNeg,:)},AM_Param,AM_Speed,UnitsPertNeg,AngleCombo)

%% PLOT CURVES ON SESSION BY SESSION BASIS. ONLY IF NECESSARY
% show different angle pop. responses for stimulus onset
% show different angle pop. responses for perturbation onset
% all units and only pert. responsive units
param_sel = AM_Param(:,AM_UOI,:);
RecordingsNR = unique(param_sel(1,:,1));
for i = 1:length(RecordingsNR)
    RecordingsNR_i(i) = min(find(param_sel(1,:,1)==RecordingsNR(i)));
end
GratingDir = 0:45:315;
AM_UnitResponses_sm = AM_UnitResponses_smooth(:,SelectedCells,:)*BonvisionFR;
clear AngResponse AngResponsePert AngResponsePertOFF
for Rec_OI = 1:length(RecordingsNR_i)
    % responses of all units 
    clear SingleSession
    if Rec_OI == length(RecordingsNR_i)
        SingleSession = AM_UnitResponses_sm(:,RecordingsNR_i(Rec_OI):end,:);
    else
        SingleSession = AM_UnitResponses_sm(:,RecordingsNR_i(Rec_OI):RecordingsNR_i(Rec_OI+1)-1,:);
    end
    for Angle_OI = 1:length(GratingDir) 
        % save the mean responses for each angle of all units in a struct
        AngResponse{Rec_OI,Angle_OI} = SingleSession(param_sel(:,RecordingsNR_i(Rec_OI),2)==GratingDir(Angle_OI),:,:);
        AngResponsePert{Rec_OI,Angle_OI} = SingleSession(param_sel(:,RecordingsNR_i(Rec_OI),2)==GratingDir(Angle_OI) & ...
                                                     param_sel(:,RecordingsNR_i(Rec_OI),3)==1,:,:);
        AngResponsePertOFF{Rec_OI,Angle_OI} = SingleSession(param_sel(:,RecordingsNR_i(Rec_OI),2)==GratingDir(Angle_OI) & ...
                                                     param_sel(:,RecordingsNR_i(Rec_OI),3)==0,:,:);
    end
end

% look at stimulus onset
% mean response - baseline for all angles. All units. Only positive units
% look at perturbation onset
% mean response - baseline for all angles. All units. Only positive units /
% focus on perturbation responsive units but also look rest of units
% look also at peak response to evaluate magnitude.

%% STIMULUS ONSET
% plot all units
figure
set(gcf, 'Position', [-1920 500 3800 400]);
PertTrial_Ex(1) = 1;
PertTrial_Ex(2) = min(find(ProjectData.Trial_data{1,1}.PerturbationON==1))
TimeLine = linspace(min(ProjectData.Trial_data{PertTrial_Ex(1),1}.time_in_trial{PertTrial_Ex(2)}), ...
                    max(ProjectData.Trial_data{PertTrial_Ex(1),1}.time_in_trial{PertTrial_Ex(2)}),...
                    size(AM_UnitResponses_sm,3)) -trialSide_seconds; % -1 seconds
Contrast = ProjectData.Trial_data{PertTrial_Ex(1),1}.contrast_in_trial{PertTrial_Ex(2),1}(1:end-1);
TF = ProjectData.Trial_data{PertTrial_Ex(1),1}.TF_in_trial{PertTrial_Ex(2),1}(1:end-1);
TrialEnd = TimeLine(min(find(Contrast(length(Contrast)/2:end)==0))+round(length(Contrast)/2));
clear PopResponses
for i = 1:length(GratingDir) 
    clear AngleFilter 
    AngleFilter = param_sel(:,:,2) == GratingDir(i);
    AM_AngleResponse = AngleFilter ;
    % generate 3D matrix for including time responses
    t1 = double(repmat(AM_AngleResponse, 1, 1, size(AM_UnitResponses_sm,3))); % repeating the selection across time
    
    % nan all the zero values
    t1(t1(:)==0) = nan;
    % select only the responses during selected trials (pert ON and OFF separately)
    TrialsResp_OI = t1.*AM_UnitResponses_sm;
    TrialsResp_OI_2D = squeeze(nanmean(TrialsResp_OI,1));
    TrialsResp_max = nanmax(TrialsResp_OI_2D,[],2);
    Response = TrialsResp_OI_2D./TrialsResp_max;
    
    subplot(1,12,i)
    h1 = shadedErrorBar(TimeLine(1:size(Response,2)),...
        smooth(mean(Response(~isnan(sum(Response,2)),:),1),5),...
        std(Response(~isnan(sum(Response,2)),:),1)/sqrt(size(Response(~isnan(sum(Response,2)),:),1)), ...
        'lineprops',{'k-','markerfacecolor','k'});
    hold on
    plot([0 0],[0 1],'--r','LineWidth',0.75);
    hold on
    plot([TrialEnd TrialEnd ], [0 1],'--r','LineWidth',0.75);
    hold on
    plot([TimeLine(min(find(TF==0))), TimeLine(min(find(TF==0)))], [0 1],'--b','LineWidth',1);
    hold on
    plot([TimeLine(max(find(TF==0))), TimeLine(max(find(TF==0)))], [0 1],'--b','LineWidth',1);
    %set(gca,'FontSize',13);
    if i == 1
        ylabel('Norm. response (all units)','fontsize',13);
        xlabel('Seconds');
        %set(gca,'XTickLabel',[])
    else
        set(gca,'YTickLabel',[],'XTickLabel',[])
    end 
    box off
    set(gca,'TickDir','out');
    TotNrPertTrials = 0;
    xlim([-0.8 max(TimeLine)-0.3]);
    ylim([0.45 0.6])
    title([ num2str(GratingDir(i)) '^o']);
    
    % compute overall responses
    PopResponses(i,:) = (mean(Response(~isnan(sum(Response,2)),:),1));
end
% plot 
subplot(1,12,[9 10])
YValues = [mean(PopResponses(:,trialSide_samples:trialSide_samples+4*BonvisionFR),2)- ...
        mean(PopResponses(:,1:trialSide_samples),2)]';
plot([-90 -45 GratingDir 360], [YValues(end-1:end) YValues YValues(1)],'.-')
%ylim([0.45 0.6])
xlabel('Grating direction')
title('Stimulus onset responses');
box off
set(gca,'XTick',[-90:90:360],'XTickLabel',[270 0:90:270 0])
set(gca, 'TickDir', 'out')

subplot(1,12,[11 12])
YValues = [mean(PopResponses(:,min(find(TF==0)):max(find(TF==0))),2)- ...
           mean(PopResponses(:,min(find(TF==0))-60:min(find(TF==0))),2)]';
plot([-90 -45 GratingDir 360], [YValues(end-1:end) YValues YValues(1)],'.-')
title('Perturbation responses');
xlabel('Grating direction')
box off
set(gca,'XTick', [-90:90:360],'XTickLabel',[270 0:90:270 0])
set(gca, 'TickDir', 'out')

figure
polarplot(deg2rad(0:45:360), [YValues YValues(1)])


%% plot only pert units
Trials_PertON = param_sel(:,UnitsPertPos,3)==1; % find indexes where perturbation is on
Trials_PertOFF= param_sel(:,UnitsPertPos,3)==0; % find indexes where perturbation is off
figure
set(gcf, 'Position', [-1920 500 3800 400]);
PertTrial_Ex(1) = 1;
PertTrial_Ex(2) = min(find(ProjectData.Trial_data{1,1}.PerturbationON==1))
TimeLine = linspace(min(ProjectData.Trial_data{PertTrial_Ex(1),1}.time_in_trial{PertTrial_Ex(2)}), ...
                    max(ProjectData.Trial_data{PertTrial_Ex(1),1}.time_in_trial{PertTrial_Ex(2)}),...
                    size(AM_UnitResponses_sm,3)) -trialSide_seconds; % -1 seconds
Contrast = ProjectData.Trial_data{PertTrial_Ex(1),1}.contrast_in_trial{PertTrial_Ex(2),1}(1:end-1);
TF = ProjectData.Trial_data{PertTrial_Ex(1),1}.TF_in_trial{PertTrial_Ex(2),1}(1:end-1);
TrialEnd = TimeLine(min(find(Contrast(length(Contrast)/2:end)==0))+round(length(Contrast)/2));
for i = 1:length(GratingDir) 
    clear AngleFilter TrialsResp_OI TrialsResp_OI_control
    AngleFilter = param_sel(:,UnitsPertPos,2) == GratingDir(i);
    %AM_AngleResponse = AngleFilter ;
    AM_AngleResponse = AngleFilter & Trials_PertON;
    % find all PertOFF trials with the grating direction of interest
    AM_AngleResponseControl = AngleFilter & Trials_PertOFF;  
    % generate 3D matrix for including time responses
    t1 = double(repmat(AM_AngleResponse, 1, 1, size(AM_UnitResponses_sm(:,UnitsPertPos,:),3))); % repeating the selection across time
    t2 = double(repmat(AM_AngleResponseControl, 1, 1, size(AM_UnitResponses_sm(:,UnitsPertPos,:),3))); % repeating the selection across time
    % nan all the zero values
    t1(t1(:)==0) = nan; t2(t2(:)==0) = nan;
    % select only the responses during selected trials (pert ON and OFF separately)
    TrialsResp_OI = t1.*AM_UnitResponses_sm(:,UnitsPertPos,:);
    TrialsResp_OI_control = t2.*AM_UnitResponses_sm(:,UnitsPertPos,:);
    TrialsResp_OI_2D = squeeze(nanmean(TrialsResp_OI,1));
    TrialsResp_OI_2D_control = squeeze(nanmean(TrialsResp_OI_control,1));
    TrialsResp_max = max(nanmax(TrialsResp_OI_2D,[],2),nanmax(TrialsResp_OI_2D_control,[],2));
    Response = TrialsResp_OI_2D./TrialsResp_max;
    Response_control = TrialsResp_OI_2D_control./TrialsResp_max;
    
    subplot(1,12,i)
    h1 = shadedErrorBar(TimeLine(1:size(Response,2)),...
        smooth(mean(Response(~isnan(sum(Response,2)),:),1),5),...
        std(Response(~isnan(sum(Response,2)),:),1)/sqrt(size(Response(~isnan(sum(Response,2)),:),1)), ...
        'lineprops',{'r-','markerfacecolor','r'});
    h2 = shadedErrorBar(TimeLine(1:size(Response_control,2)),...
        smooth(mean(Response_control(~isnan(sum(Response_control,2)),:),1),5),...
        std(Response_control(~isnan(sum(Response_control,2)),:),1)/sqrt(size(Response_control(~isnan(sum(Response_control,2)),:),1)), ...
        'lineprops',{'k-','markerfacecolor','k'});
    hold on
    plot([0 0],[0 1],'--r','LineWidth',0.75);
    hold on
    plot([TrialEnd TrialEnd ], [0 1],'--r','LineWidth',0.75);
    hold on
    plot([TimeLine(min(find(TF==0))), TimeLine(min(find(TF==0)))], [0 1],'--b','LineWidth',1);
    hold on
    plot([TimeLine(max(find(TF==0))), TimeLine(max(find(TF==0)))], [0 1],'--b','LineWidth',1);
    %set(gca,'FontSize',13);
    if i == 1
        ylabel('Norm. response (all units)','fontsize',13);
        xlabel('Seconds');
        %set(gca,'XTickLabel',[])
    else
        set(gca,'YTickLabel',[],'XTickLabel',[])
    end 
    box off
    set(gca,'TickDir','out');
    TotNrPertTrials = 0;
    xlim([-0.8 max(TimeLine)-0.3]);
    
    title([ num2str(GratingDir(i)) '^o']);
    
    % compute overall responses
    PopResponses(i,:) = smooth(mean(Response(~isnan(sum(Response,2)),:),1),5);
end
% plot 
subplot(1,12,[9 10])
YValues = [mean(PopResponses(:,trialSide_samples:trialSide_samples+4*BonvisionFR),2)- ...
        mean(PopResponses(:,1:trialSide_samples),2)]';
plot([-90 -45 GratingDir 360], [YValues(end-1:end) YValues YValues(1)],'.-')
%ylim([0.45 0.6])
xlabel('Grating direction')
title('Stimulus onset responses');
box off
set(gca,'XTick',[-90:90:360],'XTickLabel',[270 0:90:270 0])
set(gca, 'TickDir', 'out')

subplot(1,12,[11 12])
YValues = [mean(PopResponses(:,min(find(TF==0)):max(find(TF==0))),2)- ...
           mean(PopResponses(:,min(find(TF==0))-60:min(find(TF==0))),2)]';
plot([-90 -45 GratingDir 360], [YValues(end-1:end) YValues YValues(1)],'.-')
title('Perturbation onset responses');
xlabel('Grating direction')
box off
set(gca,'XTick',[-90:90:360],'XTickLabel',[270 0:90:270 0])
set(gca, 'TickDir', 'out')    
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SESSION BY SESSION 
% All sessions
for i = 1:size(AngResponse,1) 
    figure('Name', ['Recording ' num2str(i)],'NumberTitle','off')
    set(gcf, 'Position', [-1920 500 3800 400]);
    clear PopResponses
    for j = 1:size(AngResponse,2) 
        Response = AngResponsePert{i,j}; % or AngResponse for everything
        Response = squeeze(nanmean(Response,1));
        subplot(1,12,j)
        if size(Response,2)>1
            h1 = shadedErrorBar(TimeLine(1:size(Response,2)),...
                                smooth(mean(Response(~isnan(sum(Response,2)),:),1),5),...
                                std(Response(~isnan(sum(Response,2)),:),1)/sqrt(size(Response(~isnan(sum(Response,2)),:),1)), ...
                                'lineprops',{'k-','markerfacecolor','k'});
        else
            plot(TimeLine,smooth(Response,5),'k-');
        end
        hold on
        plot([0 0],[0 1],'--r','LineWidth',0.75);
        hold on
        plot([TrialEnd TrialEnd ], [0 1],'--r','LineWidth',0.75);
        hold on
        plot([TimeLine(min(find(TF==0))), TimeLine(min(find(TF==0)))], [0 1],'--b','LineWidth',1);
        hold on
        plot([TimeLine(max(find(TF==0))), TimeLine(max(find(TF==0)))], [0 1],'--b','LineWidth',1);
        %set(gca,'FontSize',13);
        if j == 1
            ylabel('Norm. response (all units)','fontsize',13);
            xlabel('Seconds');
            %set(gca,'XTickLabel',[])
        else
            set(gca,'YTickLabel',[],'XTickLabel',[])
        end
        box off
        set(gca,'TickDir','out');
        TotNrPertTrials = 0;
        xlim([-0.8 max(TimeLine)-0.3]);
        ylim([0 0.2])
        title([ num2str(GratingDir(j)) '^o']);
        % compute overall responses
        if size(Response,2)>1
            PopResponses(j,:) = smooth(mean(Response(~isnan(sum(Response,2)),:),1),5);
        else
            PopResponses(j,:) = smooth(Response,5);
        end
    end
    subplot(1,12,[9 10])
    YValues = [mean(PopResponses(:,trialSide_samples:trialSide_samples+4*BonvisionFR),2)- ...
        mean(PopResponses(:,1:trialSide_samples),2)]';
    plot([-90 -45 GratingDir 360], [YValues(end-1:end) YValues YValues(1)],'.-')
    %ylim([0.45 0.6])
    xlabel('Grating direction')
    title('Stimulus onset responses');
    box off
    set(gca,'XTick',[-90:90:360],'XTickLabel',[270 0:90:270 0])
    set(gca, 'TickDir', 'out')
    
    subplot(1,12,[11 12])
    YValues = [mean(PopResponses(:,min(find(TF==0)):max(find(TF==0))),2)- ...
        mean(PopResponses(:,min(find(TF==0))-60:min(find(TF==0))),2)]';
    plot([-90 -45 GratingDir 360], [YValues(end-1:end) YValues YValues(1)],'.-')
    title('Perturbation responses');
    xlabel('Grating direction')
    box off
    set(gca,'XTick',[-90:90:360],'XTickLabel',[270 0:90:270 0])
    set(gca, 'TickDir', 'out')
    saveas(gcf,['D:\Dropbox\UCL\Tomaso-Aman\VisPerturbation\Fig_2\panels\Sessions' filesep 'RecordingPertOnly_' num2str(i) '.png' ])
    close all
end


%% Only sessions with perturbation resp. units
Sessions_OI = unique(param_sel(1,UnitsPertPos,1));
clear Sessions_Units
Sessions_Units(:,2) = param_sel(1,:,1);
Sessions_Units(:,1) = UnitsPertPos;
%find(Sessions_Units(find(Sessions_Units(:,2) == Sessions_OI(i)),1)==1)
for i = 1:length(Sessions_OI) 
    figure('Name', ['Recording ' num2str(i)],'NumberTitle','off')
    set(gcf, 'Position', [-1920 500 3800 400]);
    clear PopResponses
    Units2plot = find(Sessions_Units(find(Sessions_Units(:,2) == Sessions_OI(i)),1)==1);
    for j = 1:size(AngResponse,2) 
        Response = AngResponsePert{Sessions_OI(i),j};
        Response = squeeze(nanmean(Response(:,Units2plot,:),1));
        subplot(1,12,j)
        if size(Response,2)>1
            h1 = shadedErrorBar(TimeLine(1:size(Response,2)),...
                                smooth(mean(Response(~isnan(sum(Response,2)),:),1),5),...
                                std(Response(~isnan(sum(Response,2)),:),1)/sqrt(size(Response(~isnan(sum(Response,2)),:),1)), ...
                                'lineprops',{'r-','markerfacecolor','r'});
        else
            plot(TimeLine,smooth(Response,5),'r-');
        end
        hold on
        plot([0 0],[0 1],'--r','LineWidth',0.75);
        hold on
        plot([TrialEnd TrialEnd], [0 1],'--r','LineWidth',0.75);
        hold on
        plot([TimeLine(min(find(TF==0))), TimeLine(min(find(TF==0)))], [0 1],'--b','LineWidth',1);
        hold on
        plot([TimeLine(max(find(TF==0))), TimeLine(max(find(TF==0)))], [0 1],'--b','LineWidth',1);
        %set(gca,'FontSize',13);
        if j == 1
            ylabel('Norm. response (all units)','fontsize',13);
            xlabel('Seconds');
            %set(gca,'XTickLabel',[])
        else
            set(gca,'YTickLabel',[],'XTickLabel',[])
        end
        box off
        set(gca,'TickDir','out');
        TotNrPertTrials = 0;
        xlim([-0.8 max(TimeLine)-0.3]);
        %ylim([0 0.2])
        title([ num2str(GratingDir(j)) '^o']);
        % compute overall responses
        if size(Response,2)>1
            PopResponses(j,:) = smooth(mean(Response(~isnan(sum(Response,2)),:),1),5);
        else
            PopResponses(j,:) = smooth(Response,5);
        end
    end
    subplot(1,12,[9 10])
    YValues = [mean(PopResponses(:,trialSide_samples:trialSide_samples+4*BonvisionFR),2)- ...
        mean(PopResponses(:,1:trialSide_samples),2)]';
    plot([-90 -45 GratingDir 360], [YValues(end-1:end) YValues YValues(1)],'.-')
    %ylim([0.45 0.6])
    xlabel('Grating direction')
    title('Stimulus onset responses');
    box off
    set(gca,'XTick',[-90:90:360],'XTickLabel',[270 0:90:270 0])
    set(gca, 'TickDir', 'out')
    
    subplot(1,12,[11 12])
    YValues = [mean(PopResponses(:,min(find(TF==0)):max(find(TF==0))),2)- ...
        mean(PopResponses(:,min(find(TF==0))-60:min(find(TF==0))),2)]';
    plot([-90 -45 GratingDir 360], [YValues(end-1:end) YValues YValues(1)],'.-')
    title('Perturbation responses');
    xlabel('Grating direction')
    box off
    set(gca,'XTick',[-90:90:360],'XTickLabel',[270 0:90:270 0])
    set(gca, 'TickDir', 'out')
    saveas(gcf,['D:\Dropbox\UCL\Tomaso-Aman\VisPerturbation\Fig_2\panels\Sessions' filesep 'RecordingPertUnits_' num2str(Sessions_OI(i)) '.png' ])
    close all
end

%% Directional responses of ensembles of pert. responsive units
Sessions_OI = unique(param_sel(1,UnitsPertPos,1));
clear Sessions_Units
Sessions_Units(:,2) = param_sel(1,:,1);
Sessions_Units(:,1) = UnitsPertPos;
clear YValues1 YValues2 
for i = 1:length(Sessions_OI) 
    clear PopResponses
    Units2plot = find(Sessions_Units(find(Sessions_Units(:,2) == Sessions_OI(i)),1)==1);
    for j = 1:size(AngResponse,2) 
        Response = AngResponse{Sessions_OI(i),j};
        Response = squeeze(nanmean(Response(:,Units2plot,:),1));
        if size(Response,2)>1
            PopResponses(j,:) = mean(Response(~isnan(sum(Response,2)),:),1);
        else
            PopResponses(j,:) = Response;
        end
    end
        % stimulus onset
        YValues1(i,:) = [mean(PopResponses(:,trialSide_samples:trialSide_samples+4*BonvisionFR),2)- ...
                    mean(PopResponses(:,1:trialSide_samples),2)]';
        % perturbation onset
        YValues2(i,:) = [mean(PopResponses(:,min(find(TF==0)):max(find(TF==0))),2)- ...
                    mean(PopResponses(:,min(find(TF==0))-60:min(find(TF==0))),2)]';
        EnsembleSize(i) = size(Response,2);
end

% figure
% boxplot(YValues2(find(EnsembleSize~=1),1)-YValues2(find(EnsembleSize~=1),:),'Notch','on','Labels',{'0','45','90','135','180','225','270','315'},'symbol','')
figure
Data2Scatter = YValues2(find(EnsembleSize~=1),1)-YValues2(find(EnsembleSize~=1),:);
UnivarScatter(Data2Scatter(:,2:end),'Label',{'45','90','135','180','225','270','315'},'MarkerFaceColor',[ 0.5 0.5 0.5])

xlabel('Perturbation direction');
ylabel('Differential response');
set(gca,'TickDir','out','box','off')
title('Direction')
ylim([-max(max(abs(YValues2(find(EnsembleSize~=1),1)-YValues2(find(EnsembleSize~=1),:)))) ...
       max(max(abs(YValues2(find(EnsembleSize~=1),1)-YValues2(find(EnsembleSize~=1),:))))]);

figure
YValuesOri =   [[YValues2(find(EnsembleSize~=1),1); YValues2(find(EnsembleSize~=1),5)], ...
                [YValues2(find(EnsembleSize~=1),2); YValues2(find(EnsembleSize~=1),6)], ...
                [YValues2(find(EnsembleSize~=1),3); YValues2(find(EnsembleSize~=1),7)], ...
                [YValues2(find(EnsembleSize~=1),4); YValues2(find(EnsembleSize~=1),8)]];
% boxplot(YValuesOri(:,1) - YValuesOri,'Notch','on','Labels',{'0','45','90','135'},'symbol','')
Data2Scatter = YValuesOri(:,1) - YValuesOri;
UnivarScatter(Data2Scatter(:,2:end),'Label',{'45','90','135'},'MarkerFaceColor',[ 0.5 0.5 0.5])

% hold on
% for i = 1:size(YValuesOri,2)
%     plot(i-0.1+0.2*rand(size(YValuesOri,1),1), YValuesOri(:,1)-YValuesOri(:,i), '.', 'Color',[0.5 0.5 0.5]);
% end
title('Orientation')
xlabel('Perturbation orientation');
ylabel('Differential response');
set(gca,'TickDir','out','box','off')
hold on
ylim([-max(max(abs(YValuesOri(:,1) - YValuesOri))) max(max(abs(YValuesOri(:,1) - YValuesOri)))])


%% UNIT BY UNIT
% plot all perturbation responsive units
% plot responses in rows. Stimulus onset and then perturbation
Sessions_OI = unique(param_sel(1,UnitsPertPos,1));
clear Sessions_Units
Sessions_Units(:,2) = param_sel(1,:,1);
Sessions_Units(:,1) = UnitsPertPos;
for i = 1:length(Sessions_OI) 
    clear PopResponses
    Units2plot = find(Sessions_Units(find(Sessions_Units(:,2) == Sessions_OI(i)),1)==1);
    for k = 1:length(Units2plot)      
        figure
        set(gcf, 'Position', [-1920 500 3800 200]);
        for j = 1:size(AngResponse,2)
            clear StimResponses PertResponses PertResponsesOff
            % save responses for stimulus onset
            Response = AngResponse{Sessions_OI(i),j};
            StimResponses = squeeze(Response(:,Units2plot(k),:));
            % save responses for perturbation onset
            Response = AngResponsePert{Sessions_OI(i),j};
            PertResponses = squeeze(Response(:,Units2plot(k),:));
            Response = AngResponsePertOFF{Sessions_OI(i),j};
            PertResponsesOff = squeeze(Response(:,Units2plot(k),:));
            
            subplot(1,10,j)
            % plot stim onset responses
            shadedErrorBar(TimeLine(1:trialSide_samples+3*BonvisionFR), ...
                nanmean(StimResponses(:,1:trialSide_samples+3*BonvisionFR)),...
                nanstd(StimResponses(:,1:trialSide_samples+3*BonvisionFR))/sqrt(size(StimResponses(:,1:trialSide_samples+3*BonvisionFR),1)),...
                'lineprop',{'k-','markerfacecolor',[0.5 0.5 0.5]});
            hold on
            % plot perturbation responses
            shadedErrorBar(TimeLine(min(find(TF==0))-BonvisionFR:end), ...
                nanmean(PertResponses(:,min(find(TF==0))-BonvisionFR:end)),...
                nanstd(PertResponses(:,min(find(TF==0))-BonvisionFR:end))/sqrt(size(PertResponses(:,min(find(TF==0))-BonvisionFR:end),1)),...
                'lineprop',{'r-','markerfacecolor','r'});
            hold on
            shadedErrorBar(TimeLine(min(find(TF==0))-BonvisionFR:end), ...
                nanmean(PertResponsesOff(:,min(find(TF==0))-BonvisionFR:end)),...
                nanstd(PertResponsesOff(:,min(find(TF==0))-BonvisionFR:end))/sqrt(size(PertResponsesOff(:,min(find(TF==0))-BonvisionFR:end),1)),...
                'lineprop',{'k-','markerfacecolor',[0.5 0.5 0.5]});
            hold on
            plot([0 0],[0 1],'--r','LineWidth',0.75);
            hold on
            plot([TrialEnd TrialEnd ], [0 1],'--r','LineWidth',0.75);
            hold on
            plot([TimeLine(min(find(TF==0))), TimeLine(min(find(TF==0)))], [0 1],'--b','LineWidth',1);
            hold on
            plot([TimeLine(max(find(TF==0))), TimeLine(max(find(TF==0)))], [0 1],'--b','LineWidth',1);
            
        end    
    end
       
end





%%
% select trials of interest and control trials as well
Trials_PertON = AM_Param(:,:,3)==1; % find indexes where perturbation is on
Trials_PertOFF= AM_Param(:,:,3)==0; % find indexes where perturbation is off
% OSI, DSI, etc. for perturbation responses
% find angle direction of each trial
Param = AM_Param(:,AM_UOI,2);
PertON = Trials_PertON(:,AM_UOI);
UnitResponses = AM_UnitResponses_smooth(:,AM_UOI,:)*60;
% for t = 1:size(UnitResponses,2)
%     OneUnitResponse = squeeze(UnitResponses(:,t,:));
%     UnitResponses_n(:,t,:) = OneUnitResponse./nanmax(OneUnitResponse,[],2);
% end
% UnitResponses = UnitResponses_n;
GratingDirections = unique(Param(~isnan(Param(:,1)),1));
clear Gr_Dir_2D t1 TrialsResp_OI
for i = 1:length(GratingDirections)
   Gr_Dir_2D(:,:,i) = Param==GratingDirections(i); 
   Gr_Dir_2D(:,:,i) = Gr_Dir_2D(:,:,i).*PertON;
   t1 = double(repmat(Gr_Dir_2D(:,:,i), 1, 1, size(UnitResponses,3))); % repeating the selection across time
   t1(t1(:)==0) = nan; 
   TrialsResp_OI{i} = t1.*UnitResponses;
end

% go through responses of single untis
ResponseTimes = AM_Param(:,AM_UOI,4:5);
BonvisionFR = 60; %Hz
trialSide_samples = 60;
trialSide_seconds = 1;
clear GratingDirResp GratingDirResp_BL
for k = 1:size(UnitResponses,2)
    clear UnitResp SingleUnitResponse_oneDir
    for i = 1:length(GratingDirections)
        ResponseTime = round(([nanmin(ResponseTimes(:,k,1)) nanmax(ResponseTimes(:,k,2))])*BonvisionFR)+trialSide_samples;
        SingleUnitResponse_oneDir = squeeze(TrialsResp_OI{i}(~isnan(TrialsResp_OI{i}(:,k,1)),k,ResponseTime(1):ResponseTime(2)));
        %SingleUnitResponse_oneDir = SingleUnitResponse_oneDir./max(SingleUnitResponse_oneDir);
        if size(SingleUnitResponse_oneDir,2)==1
            BL = nanmean(squeeze(TrialsResp_OI{i}(~isnan(TrialsResp_OI{i}(:,k,1)),k,ResponseTime(1)-trialSide_samples:ResponseTime(1))));
        else
            BL = nanmean(squeeze(TrialsResp_OI{i}(~isnan(TrialsResp_OI{i}(:,k,1)),k,ResponseTime(1)-trialSide_samples:ResponseTime(1))),2);
        end
        UnitResp(i) = nanmean(nanmean(SingleUnitResponse_oneDir,2))-nanmean(BL);
        UnitResp_BL(i) = nanmean(BL);
    end
    UnitResp(isnan(UnitResp)) = 0;
    %UnitResp = UnitResp+abs(min(UnitResp));
    GratingDirResp(k,:) = UnitResp;
    GratingDirResp_BL(k,:) = UnitResp_BL;
end


%% check whether the unit is negatively modulated by perturbation
% flip the tuning curve if modulation index DM is negative
load('DM_pert_shuffled.mat')
clear DM
DM = DM_sh(:,1);
for d = 1:size(GratingDirResp,1)
    if DM(d)<=0
       GratingDirResp(d,:) = -(GratingDirResp(d,:)-max(GratingDirResp(d,:)));
    end
end

%% Compute OSI and DSI for each unit            
clear Pref_Dir OSI DSI L_ori P_ori L_dir P_dir
for k = 1:size(UnitResponses,2)
    [R_pref ind] = max(GratingDirResp(k,:));
    R_null = GratingDirResp(k,max(mod(ind+4,length(GratingDirections)),1));
    R_ortho1 = GratingDirResp(k,max(mod(ind+2,length(GratingDirections)),1));
    R_ortho2 = GratingDirResp(k,max(mod(ind+6,length(GratingDirections)),1));
    Pref_Dir(k,1) = GratingDirections(ind);% perferred direction; 
    % ORI = (R_pref+R_null)-(R_ortho1+R_ortho2)/(R_pref+R_null)+(R_ortho1+R_ortho2)
    OSI(k,1) = ((R_pref+R_null)-(R_ortho1+R_ortho2))/((R_pref+R_null)+(R_ortho1+R_ortho2)); 
    % DSI = (R_pref-R_null)/(R_pref+R_null)
    DSI(k,1) = (R_pref-R_null)/(R_pref+R_null);
    % circular variance or magnitude of the vector L = 1-CirVar, ORI space
    L_ori(k,1) = norm( GratingDirResp(k,:)*exp(2i*degtorad(GratingDirections))/sum(GratingDirResp(k,:)) ) ;
    P_ori(k,1) = rad2deg(angle(GratingDirResp(k,:)*exp(2i*degtorad(GratingDirections))/sum(GratingDirResp(k,:))))/2;
    if sign(P_ori(k,1))==-1
        P_ori(k,1) = P_ori(k,1)+180;
    end
    % circular variance or magnitude of the vector L = 1-CirVar, DIR space
    L_dir(k,1) = norm( GratingDirResp(k,:)*exp(1i*degtorad(GratingDirections))/sum(GratingDirResp(k,:)) ) ;
    P_dir(k,1) = rad2deg(angle(GratingDirResp(k,:)*exp(1i*degtorad(GratingDirections))/sum(GratingDirResp(k,:))));
    if sign(P_dir(k,1))==-1
        P_dir(k,1) = 360+P_dir(k,1);
    end
end
TuningProps = table(Pref_Dir,OSI,DSI,L_ori,P_ori,L_dir,P_dir);

%% plot direction and orientation distribution
% orientation
figure
histogram(TuningProps.P_ori,0:10:180,'Normalization','probability')
[n_bins edges] = histcounts(TuningProps.P_ori,0:10:180);
figure
polarplot(n_bins)

% direction
figure
histogram(TuningProps.P_dir,0:10:360,'Normalization','probability')
[n_bins edges] = histcounts(TuningProps.P_dir,0:10:360);
figure
polarplot(n_bins)

%% compute significance of response for Orientation tuning : Hotelling's T-squared test
clear GratingDirResp_trial
for k = 1:size(UnitResponses,2)
    clear UnitResp_all SingleUnitResponse_oneDir
    for i = 1:length(GratingDirections)
        ResponseTime = round(([nanmin(ResponseTimes(:,k,1)) nanmax(ResponseTimes(:,k,2))])*BonvisionFR);
        % response to a single direction
        SingleUnitResponse_oneDir = squeeze(TrialsResp_OI{i}(~isnan(TrialsResp_OI{i}(:,k,1)),k,ResponseTime(1):ResponseTime(2)));
        %SingleUnitResponse_oneDir = SingleUnitResponse_oneDir./max(SingleUnitResponse_oneDir);
        % response pre-perturbation or baseline (BL)
        if size(SingleUnitResponse_oneDir,2)==1
            BL = nanmean(squeeze(TrialsResp_OI{i}(~isnan(TrialsResp_OI{i}(:,k,1)),k,ResponseTime(1)-trialSide_samples:ResponseTime(1))));
            % save responses for each direction - BL activity
            UnitResp_all{i} = nanmean(SingleUnitResponse_oneDir-BL);
        else
            BL = nanmean(squeeze(TrialsResp_OI{i}(~isnan(TrialsResp_OI{i}(:,k,1)),k,ResponseTime(1)-trialSide_samples:ResponseTime(1))),2);
            % save responses for each direction - BL activity
            UnitResp_all{i} = nanmean(SingleUnitResponse_oneDir-BL,2);
        end
        % save nr of trials for each angle direction
        trials_nr(i) = size(UnitResp_all{i},1);
    end
    clear TrialResp
    for ll = 1:length(GratingDirections)
        TrialResp(ll,:) = UnitResp_all{ll}(1:min(trials_nr));
        TrialResp(ll,isnan(TrialResp(ll,:))) = 0;
        TrialResp(ll,:) = TrialResp(ll,:)+abs(min(TrialResp(ll,:)));
    end
    GratingDirResp_trial{k} = TrialResp;
end

for i = 1:length(GratingDirResp_trial)
    clear Ori_vec 
    for j = 1:size(GratingDirResp_trial{i},2)
        Ori_vec(j,:) =  [real(GratingDirResp_trial{i}(:,j)'*exp(2i*degtorad(GratingDirections))) ...
            imag(GratingDirResp_trial{i}(:,j)'*exp(2i*degtorad(GratingDirections)))];
    end
    if size(Ori_vec,1)>size(Ori_vec,2)
        p_HT2(i) = HotellingT2Test(Ori_vec,0.01);
    else
        p_HT2(i) = 1;
    end
    TuningProps.p_HT2(i) = p_HT2(i);
end

%% compute significance of response for Direction tuning: direction dot product test
for k = 1:size(GratingDirResp_trial,2)
    clear UnitResp_all SingleUnitResponse_oneDir
    % step 1 : calculate the average orientation vector onto the orientation
    % axis
    OriVec = GratingDirResp(k,:)*exp(2i*degtorad(GratingDirections))/sum(GratingDirResp(k,:));
    OriAxis = acos(real(OriVec)/imag(OriVec));
    % step 2 : calculate the magnitude of the projection of each direction
    % vector onto the orientation axis
    clear Dir_vec ProjMagn
    for j = 1:size(GratingDirResp_trial{k},2)
        Dir_vec(j,:) =  [real(GratingDirResp_trial{1,k}(:,j)'*exp(1i*degtorad(GratingDirections))) ...
            imag(GratingDirResp_trial{1,k}(:,j)'*exp(2i*degtorad(GratingDirections)))];
        ProjMagn(j) = dot(Dir_vec(j,:),[real(OriVec) imag(OriVec)]);
    end
    % step 3 : compute Student's T-test on the distribution of direction dot
    % products
    [h p_tt(k)] = ttest(ProjMagn); 
    TuningProps.p_tt(k) = p_tt(k);
end

sum(p_HT2<0.01)
sum(p_tt<0.01)

%% compare magnitude of responses for angle groups [45 0 315] VS [135 180 225]
% GratingDirections = 0 45 90 135 180 225 270 315
AngleCombo = [1 ; 5]; % [0 ; 180]; 
AngleCombo = [3 ; 7];% [90 ; 270];
AngleCombo = [1 5 ; 3 7];% [0 180; 90 270];
AngleCombo = [1 2 8 ; 4 5 6];% [0 45 315; 135 180 225];
clear GratingDirResp_trial DirGroups
for k = 1:size(UnitResponses,2)
    clear UnitResp_all SingleUnitResponse_oneDir
    for i = 1:length(GratingDirections)
        ResponseTime = round(([nanmin(ResponseTimes(:,k,1)) nanmax(ResponseTimes(:,k,2))])*BonvisionFR);
        % response to a single direction
        SingleUnitResponse_oneDir = squeeze(TrialsResp_OI{i}(~isnan(TrialsResp_OI{i}(:,k,1)),k,ResponseTime(1):ResponseTime(2)));
        %SingleUnitResponse_oneDir = SingleUnitResponse_oneDir./max(SingleUnitResponse_oneDir);
        % response pre-perturbation or baseline (BL)
        if size(SingleUnitResponse_oneDir,2)==1
            BL = nanmean(squeeze(TrialsResp_OI{i}(~isnan(TrialsResp_OI{i}(:,k,1)),k,ResponseTime(1)-trialSide_samples:ResponseTime(1))));
            % save responses for each direction - BL activity
            UnitResp_all{i} = nanmean(SingleUnitResponse_oneDir-BL);
        else
            BL = nanmean(squeeze(TrialsResp_OI{i}(~isnan(TrialsResp_OI{i}(:,k,1)),k,ResponseTime(1)-trialSide_samples:ResponseTime(1))),2);
            % save responses for each direction - BL activity
            UnitResp_all{i} = nanmean(SingleUnitResponse_oneDir-BL,2);
        end
        % save nr of trials for each angle direction
        trials_nr(i) = size(UnitResp_all{i},1);
    end
    DirGroups(k,:) = [mean(cell2mat(UnitResp_all(1,AngleCombo(1,:))')) mean(cell2mat(UnitResp_all(1,AngleCombo(2,:))'))];
    
%   GratingDirResp_trial{k} = TrialResp;
end

figure
boxplot(DirGroups)
xlabel('run dir. and run anti-dir')
figure
histogram(DirGroups(:,1)-DirGroups(:,2))
xlabel('Mean magnitude difference')

% look at tuning curves 
plot_rows = 8;
plot_cols = 12;
for u = 1:plot_rows*plot_cols%size(GratingDirResp,1)
    if mod(u,plot_rows*plot_cols)==1
        figure
        suptitle('perturbation response (norm.)')
    end
    nr_plot = mod(u,plot_rows*plot_cols);
    if nr_plot==0
        subplot(plot_rows,plot_cols,plot_rows*plot_cols)
    else
        subplot(plot_rows,plot_cols,nr_plot)
    end
    polarplot([deg2rad(GratingDirections); 0],[GratingDirResp(u,:) GratingDirResp(u,1)])
    thetaticklabels([])
end

figure
plot(GratingDirResp(14,:))
hold on
plot(-(GratingDirResp(14,:)-max(GratingDirResp(14,:))))




