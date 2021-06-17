%% Tomaso Muzzu - UCL - 19/03/2020 - Figure 1 suppl 1A

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

%%
p_value = 95;

%% first 7 animals
thres = 95;
if  size(ProjectData,1)==37
    CTRL_exp = 0;
    Animal_1st_idx = [1 5 7 12 15 24 31];   
    if ~exist('PertResp_units','var')
        % select only perturbation responsive units
        load('AUC_shuffled.mat')
        Sh_responses = AUC_shuffled(:,2:end);
        p_pert_th = prctile(Sh_responses(:),thres);
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
elseif size(ProjectData,1)==10
    CTRL_exp = 1;
    % naive animals
    Animal_1st_idx = [1 4 7];
    % select only perturbation responsive units
    load('AUC_shuffled_CTRL_1.mat')
    Sh_responses = AUC_shuffled(:,2:end);
    p_pert_th = prctile(Sh_responses(:),thres);
    PertResp_units = (AUC_shuffled(:,1)>p_pert_th);
    % select only pos. modulated perturbation responsive units
    load('DM_pert_shuffled_CTRL.mat')
    DM = DM_sh(:,1);
    DM_sign_i(:,1) = DM>0;
    DM_sign_i(:,2) = DM<=0;
    % select only pos. modulated perturbation responsive units
    PertRespUnits_pos = PertResp_units & DM_sign_i(:,1);
    PertRespUnits_neg = PertResp_units & DM_sign_i(:,2);
else
    % select only perturbation responsive units
    load('AUC_shuffled_CTRL_1.mat')
    AUC_shuffled_CTRL = AUC_shuffled;
    AUC_shuffled_pp_CTRL = AUC_shuffled_pp;
    load('AUC_shuffled.mat')
    AUC_shuffled = cat(1,AUC_shuffled,AUC_shuffled_CTRL);
    AUC_shuffled_pp = cat(1,AUC_shuffled_pp,AUC_shuffled_pp_CTRL);
    Sh_responses = AUC_shuffled(:,2:end);
    p_pert_th = prctile(Sh_responses(:),thres);
    PertResp_units = (AUC_shuffled(:,1)>p_pert_th);
    % select only pos. modulated perturbation responsive units
    load('DM_ALL.mat')
    DM_sign_i(:,1) = DM>0;
    DM_sign_i(:,2) = DM<=0;
    % select only pos. modulated perturbation responsive units
    PertRespUnits_pos = PertResp_units & DM_sign_i(:,1);
    PertRespUnits_neg = PertResp_units & DM_sign_i(:,2);
end


%% analyse speed responses to visual stimuli events

% get the speed signals from the different recordings, both normal batch
% and control batch

% compute the mean curve responses

% conditions = 1 --> nr of recording
% conditions = 2 --> grating direction [0:45:315]
% conditions = 3 --> 1 pert ON, 0 pert OFF
% conditions = 4 --> pert onset
% conditions = 5 --> pert offset

AM_ParamSel = AM_Param(:,AM_UOI,:);
AM_SpeedSel = AM_Speed(:,AM_UOI,1:559);
Rec_OI = unique(AM_ParamSel(1,PertRespUnits_pos,1)); % recordings of interest

% scan every recording and look at the mean speed responses during
% perturbation
for i = 1:length(Rec_OI)
   RecIndex(i) = min(find((AM_ParamSel(1,:,1)==Rec_OI(i))))
end
clear MeanRecSpeed SemRecSpeed
for i = 1:length(RecIndex)
   % get perturbation trials
   PertTrials_1 = AM_ParamSel(:,RecIndex(i),3)==1; 
   PertTrials_0 = AM_ParamSel(:,RecIndex(i),3)==0;
    
   MeanRecSpeed(1,:,i) = nanmean(AM_SpeedSel(PertTrials_1,RecIndex(i),:));
   MeanRecSpeed(2,:,i) = nanmean(AM_SpeedSel(PertTrials_0,RecIndex(i),:));
   SemRecSpeed(1,:,i) = nanstd(AM_SpeedSel(PertTrials_0,RecIndex(i),:))/sqrt(sum(PertTrials_1));
   SemRecSpeed(2,:,i) = nanstd(AM_SpeedSel(PertTrials_0,RecIndex(i),:))/sqrt(sum(PertTrials_0));   
end

% find beginning and end of perturbation trials
BonvisionFR = 60; %Hz  
trialSide_samples = 60; trialSide_seconds = 1; 
PertStartEnds = [round(nanmin(AM_ParamSel(:,RecIndex,4),[],1)*BonvisionFR+trialSide_samples) ;...
                round(nanmax(AM_ParamSel(:,RecIndex,5),[],1)*BonvisionFR+trialSide_samples)];
% convert these values

%% plot heatmap of speed responses for each recording
% create 2D matrix with trials ordered by condition and response strength
for i=1:length(RecIndex)
     % get perturbation trials
   PertTrials_1 = AM_ParamSel(:,RecIndex(i),3)==1; 
   PertTrials_0 = AM_ParamSel(:,RecIndex(i),3)==0;
    
   Trials_1 = squeeze(AM_SpeedSel(PertTrials_1,RecIndex(i),:));
   Trials_0 = squeeze(AM_SpeedSel(PertTrials_0,RecIndex(i),:));
   clear SortingColumns
   SortingColumns(:,1) = mean( Trials_1(:,PertStartEnds(1,i):PertStartEnds(2,i)) ,2);
   SortingColumns(:,2) = mean( Trials_1(:,trialSide_samples:trialSide_samples+4*trialSide_samples) ,2);
   [B,I_p] = sortrows(SortingColumns,[-1 -2]);
   clear SortingColumns
   SortingColumns(:,1) = mean( Trials_0(:,PertStartEnds(1,i):PertStartEnds(2,i)) ,2);
   SortingColumns(:,2) = mean( Trials_0(:,trialSide_samples:trialSide_samples+4*trialSide_samples) ,2);
   [B,I_n] = sortrows(SortingColumns,[-1 -2]);
   
   Trials{i} = cat(1,Trials_1(I_p,:),Trials_0(I_n,:)); 
end 

%%
TimeLine = linspace(-trialSide_seconds,size(MeanRecSpeed,2)/BonvisionFR-trialSide_seconds,size(MeanRecSpeed,2));
figure
suptitle('Speed responses in recordings contataining pert. responsive units')
for i=1:length(RecIndex)
    if length(RecIndex)<10
        subplot(1,length(RecIndex),i)
    else
        subplot(6,8,i)
    end
    imagesc(TimeLine,1:size(Trials{i},1),Trials{i})
    hold on
    plot([min(TimeLine) min(TimeLine)], [1 sum(AM_ParamSel(:,RecIndex(i),3)==1)], 'r', 'LineWidth',4)
    hold on
    plot([min(TimeLine) min(TimeLine)], [sum(AM_ParamSel(:,RecIndex(i),3)==1) ...
        sum(AM_ParamSel(:,RecIndex(i),3)==1)+sum(AM_ParamSel(:,RecIndex(i),3)==0)], 'k', 'LineWidth',4)
    hold on
    plot([TimeLine(PertStartEnds(1,i)) TimeLine(PertStartEnds(2,i)) TimeLine(PertStartEnds(2,i)) TimeLine(PertStartEnds(1,i)) TimeLine(PertStartEnds(1,i))], ...
            [0 0 sum(AM_ParamSel(:,RecIndex(i),3)==1) sum(AM_ParamSel(:,RecIndex(i),3)==1) 0],'r:')
    plot([0 0],[0 size(Trials{i},1)],'k:')
    plot([max(TimeLine)-1 max(TimeLine)-1],[0 size(Trials{i},1)],'k:')
    axis tight
    xlabel('seconds'); ylabel('trials')
    box off
    set(gca,'TickDir','out')
    colorbar
end

%% plot the mean curves for each recordings

clear MeanRecSpeed_grating SemRecSpeed_grating
for i = 1:length(RecIndex)
   MeanRecSpeed_grating(:,i) = nanmean(squeeze(AM_SpeedSel(~isnan(AM_SpeedSel(:,RecIndex(i),1)),RecIndex(i),:)));
   SemRecSpeed_grating(:,i) = nanstd(squeeze(AM_SpeedSel(~isnan(AM_SpeedSel(:,RecIndex(i),1)),RecIndex(i),:)))/sqrt(size(AM_SpeedSel(~isnan(AM_SpeedSel(:,RecIndex(i),1)),RecIndex(i),:),1));
end

figure
TimeLine = linspace(-trialSide_seconds,size(MeanRecSpeed,2)/BonvisionFR-trialSide_seconds,size(MeanRecSpeed,2));
suptitle('Avg. speed responses in recordings contataining pert. responsive units')
for i=1:length(RecIndex)
    if length(RecIndex)<10
        subplot(1,length(RecIndex),i)
    else
        subplot(4,6,i)
    end
    % grating responses
     shadedErrorBar( TimeLine(1:trialSide_samples+60*3),...
        MeanRecSpeed_grating(1:trialSide_samples+60*3,i),...
        SemRecSpeed_grating(1:trialSide_samples+60*3,i),...
        'lineprop',{'k-','markerfacecolor',[0.5 0.5 0.5]});
    % plot perturbation responses separately from non-pert trials
    shadedErrorBar( TimeLine(PertStartEnds(1,i)-60:PertStartEnds(2,i)+60),...
        MeanRecSpeed(1,PertStartEnds(1,i)-60:PertStartEnds(2,i)+60,i),...
        SemRecSpeed(1,PertStartEnds(1,i)-60:PertStartEnds(2,i)+60,i),...
        'lineprop',{'r-','markerfacecolor',[0.5 0.5 0.5]});
%     hold on
%     shadedErrorBar( TimeLine(PertStartEnds(1,i)-60:PertStartEnds(2,i)),...
%         MeanRecSpeed(2,PertStartEnds(1,i)-60:PertStartEnds(2,i),i),...
%         SemRecSpeed(2,PertStartEnds(1,i)-60:PertStartEnds(2,i),i),...
%         'lineprop',{'k-','markerfacecolor',[0.5 0.5 0.5]});
    ylim_mm(1) = min([min(MeanRecSpeed_grating(1:trialSide_samples+60*3,i)-SemRecSpeed_grating(1:trialSide_samples+60*3,i)) ...
                      min(MeanRecSpeed(1,PertStartEnds(1,i)-60:PertStartEnds(2,i)+60,i)-SemRecSpeed(1,PertStartEnds(1,i)-60:PertStartEnds(2,i)+60,i))]);
    ylim_mm(2) = max([max(MeanRecSpeed_grating(1:trialSide_samples+60*3,i)+SemRecSpeed_grating(1:trialSide_samples+60*3,i)) ...
                      max(MeanRecSpeed(1,PertStartEnds(1,i)-60:PertStartEnds(2,i)+60,i)+SemRecSpeed(1,PertStartEnds(1,i)-60:PertStartEnds(2,i)+60,i))]);
    hold on
    plot([TimeLine(trialSide_samples) TimeLine(trialSide_samples)], [0 ylim_mm(2)+5],'k:')
    plot([TimeLine(end-trialSide_samples) TimeLine(end-trialSide_samples)], [0 ylim_mm(2)+5],'k:')
    plot([TimeLine(PertStartEnds(1,i)) TimeLine(PertStartEnds(1,i))], [0 ylim_mm(2)+5],'r:')
    plot([TimeLine(PertStartEnds(2,i)) TimeLine(PertStartEnds(2,i))], [0 ylim_mm(2)+5],'r:')
    xlabel('seconds'); ylabel('cm/s');
    box off
    axis tight
    if i == 1
        ll = legend({'All trials', 'Only trials with Pert.'},'fontsize',10);
        ll.Color = 'none'; ll.EdgeColor = 'none';
    end
    set(gca,'TickDir','out')
    ylim([max(ylim_mm(1)-5,0) ylim_mm(2)+5]);
end

%% plot responses for each session wrt to direction

SP_coords = [95 96 109 110;...
             52 53 66 67; ...
             35 36 49 50; ...
             46 47 60 61; ...
             87 88 101 102; ...
             130 131 144 145; ...
             147 148 161 162; ...
             136 137 150 151];
SP_Polar = [76:79 90:93 104:107 118:121];  

% clear MeanRecSpeed_grating SemRecSpeed_grating
% for i = 1:length(RecIndex)
%    MeanRecSpeed_grating(:,i) = nanmean(squeeze(AM_SpeedSel(~isnan(AM_SpeedSel(:,RecIndex(i),1)),RecIndex(i),:)));
%    SemRecSpeed_grating(:,i) = nanstd(squeeze(AM_SpeedSel(~isnan(AM_SpeedSel(:,RecIndex(i),1)),RecIndex(i),:)))/sqrt(size(AM_SpeedSel(~isnan(AM_SpeedSel(:,RecIndex(i),1)),RecIndex(i),:),1));
% end
% clear MeanRecSpeed SemRecSpeed
% for i = 1:length(RecIndex)
%    % get perturbation trials
%    PertTrials_1 = AM_ParamSel(:,RecIndex(i),3)==1; 
%    PertTrials_0 = AM_ParamSel(:,RecIndex(i),3)==0;
%     
%    MeanRecSpeed(1,:,i) = nanmean(AM_SpeedSel(PertTrials_1,RecIndex(i),:));
%    MeanRecSpeed(2,:,i) = nanmean(AM_SpeedSel(PertTrials_0,RecIndex(i),:));
%    SemRecSpeed(1,:,i) = nanstd(AM_SpeedSel(PertTrials_0,RecIndex(i),:))/sqrt(sum(PertTrials_1));
%    SemRecSpeed(2,:,i) = nanstd(AM_SpeedSel(PertTrials_0,RecIndex(i),:))/sqrt(sum(PertTrials_0));   
% end

for i=1:length(RecIndex)
    
    SpeedResponse = squeeze(AM_SpeedSel(~isnan(AM_SpeedSel(:,RecIndex(i),1)),RecIndex(i),:));
    % get angles of all trials
    Trial_Angle = squeeze(AM_ParamSel(1:size(SpeedResponse,1),RecIndex(i),2));
    % get perturbation trials
    PertTrials_1 = AM_ParamSel(1:size(SpeedResponse,1),RecIndex(i),3)==1; 
    PertTrials_0 = AM_ParamSel(1:size(SpeedResponse,1),RecIndex(i),3)==0;
    % new figure
    figure
    set(gcf,'Position',[100 50 900 675])
    GratingDirectionsOrdered = unique(Trial_Angle);
    for j = 1 : length(GratingDirectionsOrdered)
        subplot(14,14,SP_coords(j,:))
        % visual
        % plot all trials
        Trial2plot = find(Trial_Angle==GratingDirectionsOrdered(j));
        shadedErrorBar(TimeLine(1:trialSide_samples+60*3),...
            mean(SpeedResponse(Trial2plot,1:trialSide_samples+60*3)),...
            std(SpeedResponse(Trial2plot,1:trialSide_samples+60*3))/sqrt(size(SpeedResponse(Trial2plot,1:trialSide_samples+60*3),1)),...
            'lineprop',{'k-','markerfacecolor',[0.5 0.5 0.5]});
        % plot mean
        hold on
        xlim([-0.8 max(TimeLine)-0.5]);
        set(gca,'box','off','TickDir','out')
        % perturbation
        Trial2plot_p = find(Trial_Angle==GratingDirectionsOrdered(j) & PertTrials_1);
        % plot mean
        if length(Trial2plot_p)~=1
            shadedErrorBar(TimeLine(PertStartEnds(1,i)-60:end),...
                mean(SpeedResponse(Trial2plot_p,PertStartEnds(1,i)-60:end),1),...
                std(SpeedResponse(Trial2plot_p,PertStartEnds(1,i)-60:end),1)/sqrt(size(SpeedResponse(Trial2plot_p,PertStartEnds(1,i)-60:end),1)),...
                'lineprop',{'r-','markerfacecolor','r'});
        else
            plot(TimeLine(PertStartEnds(1,i)-60:end),...
                (SpeedResponse(Trial2plot_p,PertStartEnds(1,i)-60:end)),...
                'r-');
        end
        Limit_Y_axis(j) = max([ max(mean(SpeedResponse(Trial2plot,1:trialSide_samples+60*3)))*1.1, ...
                                max(mean(SpeedResponse(Trial2plot_p,PertStartEnds(1,i)-60:end)))*1.1]);
        Trial_nr_p(j) = length(Trial2plot_p);
    end
    clear GratingDirRespV GratingDirResp
    for j = 1 : length(GratingDirectionsOrdered)
        subplot(14,14,SP_coords(j,:))
        ylim([0 max(Limit_Y_axis)]);
        plot([0 0], [0 max(Limit_Y_axis)],'k:'); plot([max(TimeLine)-1 max(TimeLine)-1],[0 max(Limit_Y_axis)],'k:'); % trial start and stop
        plot([TimeLine(PertStartEnds(1,i)) TimeLine(PertStartEnds(1,i))], [0,max(Limit_Y_axis)],'r:'); 
        plot([TimeLine(PertStartEnds(2,i)) TimeLine(PertStartEnds(2,i))], [0,max(Limit_Y_axis)],'r:'); % pert start and stop
%         set(gca,'XTickLabel',[]);
        if j ~=5
%             set(gca,'YTickLabel',[]);
%             set(gca,'xcolor','none','ycolor','none')
        end
        % here compute the max response so you have a tuning curve for
        % speed responses
        Trial2plot = find(Trial_Angle==GratingDirectionsOrdered(j));
        GratingDirRespV(i,j) = mean(mean(mean(SpeedResponse(Trial2plot,1:trialSide_samples+60*3))-mean(SpeedResponse(Trial2plot,1:trialSide_samples),2)));
        GratingDirResp(i,j) = mean(mean(mean(SpeedResponse(Trial2plot,PertStartEnds(1,i):PertStartEnds(2,i)))-mean(SpeedResponse(Trial2plot,PertStartEnds(1,i)-60:PertStartEnds(1,i)),2)));
    end
    %subplot(14,14,SP_Polar);
    hax = axes('Position', [.35, .35, .33, .33]);
    %         bar(hax,y,'EdgeColor','none')
    %         set(hax,'XTick',[])
    polarplot([deg2rad(GratingDirectionsOrdered); 0],max([GratingDirRespV(i,:) GratingDirRespV(i,1)],0),'k')
    % hold on
    % polarplot([deg2rad(GratingDirections); 0],[GratingDirResp(UI_idx,:) GratingDirResp(UI_idx,1)],'r')
    hold on
    GR2Plot = GratingDirResp(i,:);
    GR2Plot(GR2Plot(:)<0) = 0;
    polarplot([deg2rad(GratingDirectionsOrdered); 0],max([GratingDirResp(i,:) GratingDirResp(i,1)],0),'r')
    
    thetaticks(GratingDirectionsOrdered)
    thetaticklabels({''})
    RadiusTicks = get(gca,'RTick');
    clear RadiusTickLabels;
    for str = 1:length(RadiusTicks)-1
        RadiusTickLabels{str} = '';
    end
    RadiusTickLabels{end+1} = num2str(RadiusTicks(end));
    set(gca,'RTickLabel',RadiusTickLabels)
    % % % % % % %% % % % % % %% % % % % % %% % % % % % %
    subplot(14,14, [155 156 157 158 169 170 171 172 183 184 185 186])
     shadedErrorBar( TimeLine(1:trialSide_samples+60*3),...
        MeanRecSpeed_grating(1:trialSide_samples+60*3,i),...
        SemRecSpeed_grating(1:trialSide_samples+60*3,i),...
        'lineprop',{'k-','markerfacecolor',[0.5 0.5 0.5]});
    % plot perturbation responses separately from non-pert trials
    shadedErrorBar( TimeLine(PertStartEnds(1,i)-60:PertStartEnds(2,i)+60),...
        MeanRecSpeed(1,PertStartEnds(1,i)-60:PertStartEnds(2,i)+60,i),...
        SemRecSpeed(1,PertStartEnds(1,i)-60:PertStartEnds(2,i)+60,i),...
        'lineprop',{'r-','markerfacecolor',[0.5 0.5 0.5]});
    ylim_mm(1) = min([min(MeanRecSpeed_grating(1:trialSide_samples+60*3,i)-SemRecSpeed_grating(1:trialSide_samples+60*3,i)) ...
                      min(MeanRecSpeed(1,PertStartEnds(1,i)-60:PertStartEnds(2,i)+60,i)-SemRecSpeed(1,PertStartEnds(1,i)-60:PertStartEnds(2,i)+60,i))]);
    ylim_mm(2) = max([max(MeanRecSpeed_grating(1:trialSide_samples+60*3,i)+SemRecSpeed_grating(1:trialSide_samples+60*3,i)) ...
                      max(MeanRecSpeed(1,PertStartEnds(1,i)-60:PertStartEnds(2,i)+60,i)+SemRecSpeed(1,PertStartEnds(1,i)-60:PertStartEnds(2,i)+60,i))]);
    hold on
    plot([TimeLine(trialSide_samples) TimeLine(trialSide_samples)], [0 ylim_mm(2)+5],'k:')
    plot([TimeLine(end-trialSide_samples) TimeLine(end-trialSide_samples)], [0 ylim_mm(2)+5],'k:')
    plot([TimeLine(PertStartEnds(1,i)) TimeLine(PertStartEnds(1,i))], [0 ylim_mm(2)+5],'r:')
    plot([TimeLine(PertStartEnds(2,i)) TimeLine(PertStartEnds(2,i))], [0 ylim_mm(2)+5],'r:')
    xlabel('seconds'); ylabel('cm/s');
    box off
    axis tight
    ll = legend({'All trials', 'Only trials with Pert.'},'fontsize',10);
    ll.Color = 'none'; ll.EdgeColor = 'none';
    set(gca,'TickDir','out')
    ylim([max(ylim_mm(1)-5,0) ylim_mm(2)+5]);
    
    if length(RecIndex)>10
        saveas(gcf,['E:\Dropbox\UCL\Tomaso-Aman\VisPerturbation\Fig_1\panels\RunStillAnalysis' filesep 'SpeedResponse_Recording_' num2str(i) '.pdf']);
    else
        saveas(gcf,['E:\Dropbox\UCL\Tomaso-Aman\VisPerturbation\Fig_3\panels\Naive\RunStillAnalysis' filesep 'SpeedResponse_Recording_' num2str(i) '.pdf']);
    end
    
    close all
end




