%% Tomaso Muzzu - UCL - 7/7/20 - acceleration stimulus responses

%% load main stimulus data
if ~exist('ProjectData','var')
    [ProjectData AM_UnitResponses AM_Param AM_Speed AM_UOI SelectedResponses AM_UnitResponses_smooth] = LoadDataALL;
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
thres = 95;
CTRL_exp = 1;
% naive animals
Animal_1st_idx = [1 4 7];
% select only perturbation responsive units
load('AUC_shuffled_CTRL_1.mat')
Sh_responses = AUC_shuffled(:,2:end);
p_pert_th = prctile(Sh_responses(:),thres);
PertResp_units = (AUC_shuffled(:,1)>p_pert_th);
% select only pos. modulated perturbation responsive units
load('DM_CTRL.mat')
% DM = DM_sh(:,1);
DM_sign_i(:,1) = DM>0;
DM_sign_i(:,2) = DM<=0;
% select only pos. modulated perturbation responsive units
PertRespUnits_pos = PertResp_units & DM_sign_i(:,1);
PertRespUnits_neg = PertResp_units & DM_sign_i(:,2);
    
%% Load raw and formatted data of acceleration stimulus
[ProjectData_acc AM_UnitResponses_acc AM_Param_acc AM_Speed_acc AM_UOI_acc] = LoadData;

%% Initial selection of best units
% to be run once in theory, repeated in case params need changing
FR_thres = 0.1;
AM_UOI_acc = UnitsSelection(ProjectData_acc,AM_UnitResponses_acc,FR_thres);
% save([ DataFolder filesep 'UnitsResponse.mat'], 'AM_UnitResponses', 'AM_Speed', 'AM_Param', 'AM_UOI','-v7.3')
% ATM, units with mean FR>0.1Hz and comparable FR in first and last third
% of recording are kept

AM_UOI_bi = AM_UOI & AM_UOI_acc;

%% smooth data with Gaussian filter
% 3rd argument is width of Gaussian filter
% 4th OPTIONAL argument is width of std - sigma=1/3 width by defeault;
Gauss_width = 0.15; % in seconds
AM_UnitResponses_smooth_acc = Smooth_AM_FR(ProjectData_acc,AM_UnitResponses_acc,Gauss_width);

%% apply conditions to select responses
% conditions are set to define subset of units (AM_UOI), trials, periods within
% trials based on info store in matrix AM_Param(possible trials, units,
% conditions) where last dimension is 
% conditions = 1 --> nr of recording
% conditions = 2 --> TF at stimulus onset
% conditions = 3 --> TF at simutlus offset
% conditions = 4 --> dTF shown
% conditions = 5 --> index of when Tf starts changing
% conditions = 6 --> index of when Tf finishes changing

%% Selection of best units of main stimulus applied to acceleration protocol
%%%Cond_Sel = {} ; % {optional('run'), optional(grating dir=[0:45:315])}
Units_Sel = AM_UOI(:,1) ; % use best units selected from main stimulus 
% Units_Sel = logical(ones(size(AM_UOI(:,1),1),1));
Run_TH = 3; % running speed threshold is 1 cm/s to be used like this: ..., 'run', Run_TH)

% Cond = 0; % %% all units together, distinguishing for accel. values
% Cond = 1; % %% all units together, distringuishing for accel. values and starting TF
Cond = 2; % %% plot single units, distringuishing for accel. values
% Cond = 3; % %% plot single units, distringuishing for accel. values and starting TF
SelectedResponses = AccelResponsesSelection(ProjectData_acc,AM_UnitResponses_smooth_acc,AM_Param_acc,AM_Speed_acc(:,:,1:end-1),Units_Sel,Cond); %,'run',Run_TH);
% Units_Sel= true(449,1);
% Units_Sel(UnitPerRec(7)+1:UnitPerRec(8)) = 0;
% Units_Sel(UnitPerRec(9)+1:UnitPerRec(10)) = 0;

%%
% AM_param_acc
% conditions = 1 --> nr of recording
% conditions = 2 --> TF at stimulus onset
% conditions = 3 --> TF at simutlus offset
% conditions = 4 --> dTF shown
% conditions = 5 --> index of when Tf starts changing
% conditions = 6 --> index of when Tf finishes changing
AM_UOI = Units_Sel;

% find all the TF's at the start of the trials
StarTF_2dmask = squeeze(AM_Param_acc(:,:,2));
start_TF_values = unique(StarTF_2dmask(~isnan(StarTF_2dmask(:))));
% separate the responses for the different acceleration values
dTF_2Dmask = squeeze(AM_Param_acc(:,:,4));
dTF_values = unique(dTF_2Dmask(~isnan(dTF_2Dmask(:))));
% find trials in which dTF lasts 0.5s and 1s
dTF_dur_2D = (squeeze(AM_Param_acc(:,:,6))-squeeze(AM_Param_acc(:,:,5)))/60; % 60 is the sampling frequency
dTF_long_idx = dTF_dur_2D>0.8;
% find three trials for each condition for which there is dTF>0; dTF=0;
% dTF<0
clear TF_10_i TF_05_i I_10 I_05
for i = 1:length(dTF_values)
    TF_10_i(i) = min(find((AM_Param_acc(:,:,4).*dTF_long_idx)==dTF_values(i) & (AM_Param_acc(:,:,2).*dTF_long_idx)==6));
    TF_05_i(i) = min(find((AM_Param_acc(:,:,4).*(~dTF_long_idx))==dTF_values(i) & (AM_Param_acc(:,:,2).*(~dTF_long_idx))==6));
end
%convert to 2D indexes
ReferSize = [size(AM_Param_acc,1),size(AM_Param_acc,2)];    
[I_10(:,1),I_10(:,2)] = ind2sub(ReferSize,TF_10_i); % from linear to 2D indexes
[I_05(:,1),I_05(:,2)] = ind2sub(ReferSize,TF_05_i); % from linear to 2D indexes


% plot only the perturbation responsive units
% make a 3-by-3 plots figure for each
% 1st row, stimulus example
% 2nd row, mean responses for positive, zero, and negative
% 3rd row, heatmap of responses for all trials, sorted by final TF 
PP_I = find(PertRespUnits_pos);
AM_UnitResponses_acc_s = AM_UnitResponses_smooth_acc(:,AM_UOI,:);
AM_Param_acc_s = AM_Param_acc(:,AM_UOI,:);
for i = 1:length(PP_I)
    
    % find all the TF's at the start of the trials
    StarTF_2dmask = AM_Param_acc_s(:,PP_I(i),2)
    start_TF_values = unique(StarTF_2dmask(~isnan(StarTF_2dmask(:))));
    % separate the responses for the different acceleration values
    dTF_2Dmask = squeeze(AM_Param_acc_s(:,PP_I(i),4));
    dTF_values = unique(dTF_2Dmask(~isnan(dTF_2Dmask(:))));
    Rec_OI = AM_Param_acc_s(1,PP_I(i),1);

    figure('Renderer', 'painters', 'Position', [10 10 1800 900]);
    title(['Unit ' num2str(PP_I(i))])
    %% plot the trial parameters
    for j = 1:length(start_TF_values)
        subplot(10,3,[j j+3]); % show trials with dTF duration of 0.5s
        clear t_id
        for k = 1:size(dTF_values,1)
            t_id(:,1) = ProjectData_acc.Trial_data{Rec_OI,1}.Starting_TF == start_TF_values(j);
            t_id(:,2) = ProjectData_acc.Trial_data{Rec_OI,1}.dTF == dTF_values(k);
            y_id_plot = min(find(t_id(:,1) & t_id(:,2)));
            plot(ProjectData_acc.Trial_data{Rec_OI,1}.time_in_trial{y_id_plot,1} , ...
                ProjectData_acc.Trial_data{Rec_OI,1}.TF_in_trial{y_id_plot,1} , ...
                'Color',[0.5 0.5 0.5]);
            hold on
        end
        dTF_dur = (ProjectData_acc.Trial_data{Rec_OI,1}.Ending_TF(1)-ProjectData_acc.Trial_data{Rec_OI,1}.Starting_TF(1))/ProjectData_acc.Trial_data{Rec_OI,1}.dTF(1);
        plot([0.5 0.5],[0 9.5],'g--','LineWidth',2); plot([0.5+dTF_dur 0.5+dTF_dur],[0 9.5],'g--','LineWidth',2); % dTF period
        plot([0 0],[0 9.5],'k:','LineWidth',4); plot([3 3],[0 9.5],'k:','LineWidth',4); % trial period
        set(gca,'TickDir','out','XTickLabel',[],'XTick',[0:0.5:3])
        box off
        ylabel('TF');
        ylim([-1 9.5]); xlim([-0.5 3.5]);
        % text(0.25,1,'example trial with starting TF=6Hz & dTF=-3:1:3')
        h = text(-0.1,2, 'trial start'); set(h,'Rotation',90);
        h = text(3.1,2, 'trial end'); set(h,'Rotation',90);
        h = text(0.5+0.05,2, 'dTF start'); set(h,'Rotation',90);
        h = text(0.5+dTF_dur+0.05,2, 'dTF stop'); set(h,'Rotation',90);
        subplot(10,3,[j j+3]);
        title(['trials with start_TF=' num2str(start_TF_values(j)) 'Hz & dTF=-3:1:3; dTF dur. ' num2str(dTF_dur) 's']);
    end
    %% plot the mean responses for final speeds below 3 Hz, at 3 Hz and above 3 Hz
    c = colormap('gray');
    for j = 1:length(start_TF_values)
        subplot(10,3,[j+6:3:j+6+6]); % show trials with dTF duration of 0.5s
        clear t_id
        t_id_tf = find(ProjectData_acc.Trial_data{Rec_OI,1}.Starting_TF == start_TF_values(j));
        t_id_dtf = ProjectData_acc.Trial_data{Rec_OI,1}.Ending_TF(t_id_tf) ;
        SingleUnitResponse = squeeze(AM_UnitResponses_acc_s(t_id_tf,PP_I(i),:));
        end_TF_values = unique(t_id_dtf);
        for g = 1:length(end_TF_values)
            shadedErrorBar( linspace(-1,4,size(SingleUnitResponse,2)), ...
                        nanmean(SingleUnitResponse(t_id_dtf==end_TF_values(g),:)), ...
                        nanstd(SingleUnitResponse(t_id_dtf==end_TF_values(g),:))/sqrt(size(SingleUnitResponse(t_id_dtf==end_TF_values(g),:),1)), ...
                        'lineprops',{'color',1-c(round(size(c,1)/(length(end_TF_values)+2)*g),:)} );
            hold on
        end
        plot([0.5 0.5],[0 9.5],'g--','LineWidth',2); plot([0.5+dTF_dur 0.5+dTF_dur],[0 9.5],'g--','LineWidth',2); % dTF period
        plot([0 0],[0 9.5],'k:','LineWidth',4); plot([3 3],[0 9.5],'k:','LineWidth',4); % trial period
        set(gca,'TickDir','out','XTickLabel',[],'XTick',[0:0.5:3])
        box off
        ylabel('Spikes/s');
        xlim([-0.5 3.5]);
        % text(0.25,1,'example trial with starting TF=6Hz & dTF=-3:1:3')
        h = text(-0.1,2, 'trial start'); set(h,'Rotation',90);
        h = text(3.1,2, 'trial end'); set(h,'Rotation',90);
        h = text(0.5+0.05,2, 'dTF start'); set(h,'Rotation',90);
        h = text(0.5+dTF_dur+0.05,2, 'dTF stop'); set(h,'Rotation',90);
        legend(num2str(end_TF_values),'Location','northwest','NumColumns',3)
    end
    %% plot responses for all trials in decreasing order of final speed value
    for j = 1:length(start_TF_values)
        subplot(10,3,[j+15:3:j+15+12]); % show t    rials with dTF duration of 0.5s
        clear t_id
        t_id_tf = find(ProjectData_acc.Trial_data{Rec_OI,1}.Starting_TF == start_TF_values(j));
        t_id_dtf = ProjectData_acc.Trial_data{Rec_OI,1}.Ending_TF(t_id_tf) ;
        SingleUnitResponse = squeeze(AM_UnitResponses_acc_s(t_id_tf,PP_I(i),:));
        % sort the responses
        [B I_inc] = sort(t_id_dtf,'descend');
        imagesc(linspace(-1,4,size(SingleUnitResponse,2)), 1:size(SingleUnitResponse(I_inc,:),1) , SingleUnitResponse(I_inc,:) );
        hold on
        c = gray;
        c = flipud(c);
        colormap(c);
        [A, I_fin_TF, n] = unique(B);
        for g = 1:length(I_fin_TF)
            plot([-0.5 3.5],[I_fin_TF(g) I_fin_TF(g)],'Color',[0 0 0],'LineWidth',1)
            text(0,I_fin_TF(g)+5,['TF = ' num2str(A(g)) ' Hz'])
        end
        set(gca, 'YDir','reverse')
        set(gca,'TickDir','out','XTick',[0:0.5:3])
        box off
        xlabel('seconds'); ylabel('trials');
        xlim([-0.5 3.5])         
    end
    
     if length(PP_I)==62
         saveas(gcf, ['E:\Dropbox\UCL\Tomaso-Aman\VisPerturbation\Fig_3\Acceleration\PertUnits_AccelResponses' filesep ...
                'AccelResponse_' ...
                'Rec#' num2str(Rec_OI) ...
                '_Unit#' num2str(PP_I(i))...
                '.png']);
     elseif length(PP_I)==47
        saveas(gcf, ['E:\Dropbox\UCL\Tomaso-Aman\VisPerturbation\Fig_3\Acceleration\PertUnits_AccelResponses_neg' filesep ...
                'AccelResponse_' ...
                'Rec#' num2str(Rec_OI) ...
                '_Unit#' num2str(PP_I(i))...
                '.png']);
    end        
            
            
    close all
    
end


%% plot responses irrespective of the initial TF
% plot heat map of all responses
% plot speed tuning curve
% save info of the perturbation responses
param_sel = AM_Param(:,SelectedCells,:);
 
PP_I = find(PertRespUnits_neg);
PP_I = 1:189;
% show only reliable units on both stimuli
PP_I = find(PertRespUnits_pos & I_bi'); % I_bi computed in line ~480
c = gray;   % colormap
for i = 1:length(PP_I)
    
    % find all the TF's at the start of the trials
    StarTF_2dmask = AM_Param_acc_s(:,PP_I(i),2)
    start_TF_values = unique(StarTF_2dmask(~isnan(StarTF_2dmask(:))));
    % separate the responses for the different acceleration values
    dTF_2Dmask = squeeze(AM_Param_acc_s(:,PP_I(i),4));
    dTF_values = unique(dTF_2Dmask(~isnan(dTF_2Dmask(:))));
    Rec_OI = AM_Param_acc_s(1,PP_I(i),1);
    Trial_Angle = squeeze(param_sel(:,PP_I(i),2));
    Trial_Pert = squeeze(param_sel(:,PP_I(i),3));
    TimeLine = linspace(-1,8.33,size(AM_UnitResponses_smooth__,3));
    Pert_lims =[param_sel(min(find(param_sel(:,PP_I(i),3)==1)),PP_I(i),4) param_sel(min(find(param_sel(:,PP_I(i),3)==1)),PP_I(i),5)];

    figure('Renderer', 'painters', 'Position', [10 10 1800 900]);
    title(['Unit ' num2str(PP_I(i))])
    %% plot the trial parameters
    subplot(10,4,[1:2 5:6]); % show trials with dTF duration of 0.5s
    clear t_id
    for k = 1:size(dTF_values,1)
        t_id(:,1) = ProjectData_acc.Trial_data{Rec_OI,1}.Starting_TF == 3;
        t_id(:,2) = ProjectData_acc.Trial_data{Rec_OI,1}.dTF == dTF_values(k);
        y_id_plot = min(find(t_id(:,1) & t_id(:,2)));
        plot(ProjectData_acc.Trial_data{Rec_OI,1}.time_in_trial{y_id_plot,1} , ...
            ProjectData_acc.Trial_data{Rec_OI,1}.TF_in_trial{y_id_plot,1} , ...
            'Color',[0.5 0.5 0.5]);
        hold on
    end
    dTF_dur = (ProjectData_acc.Trial_data{Rec_OI,1}.Ending_TF(1)-ProjectData_acc.Trial_data{Rec_OI,1}.Starting_TF(1))/ProjectData_acc.Trial_data{Rec_OI,1}.dTF(1);
    plot([0.5 0.5],[0 9.5],'g--','LineWidth',2); plot([0.5+dTF_dur 0.5+dTF_dur],[0 9.5],'g--','LineWidth',2); % dTF period
    plot([0 0],[0 9.5],'k:','LineWidth',4); plot([3.02 3.02],[0 9.5],'k:','LineWidth',4); % trial period
    plot([2 2],[0 9.5],'b:','LineWidth',4); plot([3 3],[0 9.5],'b:','LineWidth',4);% period for TF tuning curve
    set(gca,'TickDir','out','XTickLabel',[],'XTick',[0:0.5:3])
    box off
    ylabel('TF');
    ylim([-1 9.5]); xlim([-0.5 3.5]);
    % text(0.25,1,'example trial with starting TF=6Hz & dTF=-3:1:3')
    h = text(-0.1,2, 'trial start'); set(h,'Rotation',90);
    h = text(3.1,2, 'trial end'); set(h,'Rotation',90);
    h = text(0.5+0.05,2, 'dTF start'); set(h,'Rotation',90);
    h = text(0.5+dTF_dur+0.05,2, 'dTF stop'); set(h,'Rotation',90);
    h = text(2+0.05,8, 'Period TF tuning curve'); set(h,'Rotation',0);
    title(['trials with any starting TF & dTF=-3:1:3; dTF dur. ' num2str(dTF_dur) 's']);
    
    %% plot the mean responses for final speeds below 3 Hz, at 3 Hz and above 3 Hz
    subplot(10,4,[9:10 13:14]); % show trials with dTF duration of 0.5s
    clear t_id
    %t_id_tf = find(ProjectData_acc.Trial_data{Rec_OI,1}.Starting_TF == start_TF_values(j));
    t_id_dtf = ProjectData_acc.Trial_data{Rec_OI,1}.Ending_TF ; % final value of TF
    SingleUnitResponse = squeeze(AM_UnitResponses_acc_s(:,PP_I(i),:));
    end_TF_values = unique(t_id_dtf);
    for g = 1:length(end_TF_values)
        shadedErrorBar( linspace(-1,4,size(SingleUnitResponse,2)), ...
            nanmean(SingleUnitResponse(t_id_dtf==end_TF_values(g),:)), ...
            nanstd(SingleUnitResponse(t_id_dtf==end_TF_values(g),:))/sqrt(size(SingleUnitResponse(t_id_dtf==end_TF_values(g),:),1)), ...
            'lineprops',{'color',1-c(round(size(c,1)/(length(end_TF_values)+2)*g),:)} );
        hold on
        Y_max(g) = max(nanmean(SingleUnitResponse(t_id_dtf==end_TF_values(g),:)));
    end
    plot([0.5 0.5],[0 9.5],'g--','LineWidth',2); plot([0.5+dTF_dur 0.5+dTF_dur],[0 9.5],'g--','LineWidth',2); % dTF period
    plot([0 0],[0 9.5],'k:','LineWidth',4); plot([3 3],[0 9.5],'k:','LineWidth',4); % trial period
    set(gca,'TickDir','out','XTickLabel',[],'XTick',[0:0.5:3])
    box off
    ylabel('Spikes/s');
    xlim([-0.5 3.5]);
    % text(0.25,1,'example trial with starting TF=6Hz & dTF=-3:1:3')
    h = text(-0.1,2, 'trial start'); set(h,'Rotation',90);
    h = text(3.1,2, 'trial end'); set(h,'Rotation',90);
    h = text(0.5+0.05,2, 'dTF start'); set(h,'Rotation',90);
    h = text(0.5+dTF_dur+0.05,2, 'dTF stop'); set(h,'Rotation',90);
    h = text(0.05,max(Y_max), 'Lighter colour for higher TF'); set(h,'Rotation',0);
    %legend(num2str(end_TF_values))
    
    %% plot responses for all trials in decreasing order of final speed value
    subplot(10,4,[17:18 21:22 25:26 29:30 33:34 37:38]); % show t    rials with dTF duration of 0.5s
    clear t_id
    t_id_dtf = ProjectData_acc.Trial_data{Rec_OI,1}.Ending_TF ;
    t_id_dtfs = ProjectData_acc.Trial_data{Rec_OI,1}.Starting_TF ; 
    SingleUnitResponse = squeeze(AM_UnitResponses_acc_s(:,PP_I(i),:));
    % sort the responses
    [B I_inc] = sortrows([t_id_dtf t_id_dtfs],[-1 -2]);
    imagesc(linspace(-1,4,size(SingleUnitResponse,2)), 1:size(SingleUnitResponse(I_inc,:),1) , SingleUnitResponse(I_inc,:) );
    hold on
    c = gray;
    c = flipud(c);
    colormap(c);
    [A, I_fin_TF, n] = unique(B);
    for g = 1:length(I_fin_TF)
        plot([-0.5 3.5],[I_fin_TF(g) I_fin_TF(g)],'Color',[0 0 0],'LineWidth',1)
        text(-0.5,I_fin_TF(g)+5,[num2str(A(g)) ' Hz'])
    end
    set(gca, 'YDir','reverse')
    set(gca,'TickDir','out','XTick',[0:0.5:3])
    box off
    xlabel('seconds'); ylabel('trials');
    xlim([-0.5 3.5])
    
    %% plot speed tuning
    subplot(10,4, [3:4 7:8])
    errorbar(U_speed, Resp_mean(PP_I(i),:), Resp_sem(PP_I(i),:),'Color','k')
    set(gca,'box','off','TickDir','out'); box off
    xlabel('Temporal TF (cycles/s)'); ylabel('Firing rate (spikes/s)')
    xlim([0 max(U_speed)])
    
    %% plot responses to perturbation
    subplot(10,4, [15:16 19:20])
    errorbar(GratingDirections, GratingDirRespV(PP_I(i),:), GratingDirRespV_sem(PP_I(i),:), 'k');
    hold on
    errorbar(GratingDirections, GratingDirResp(PP_I(i),:), GratingDirResp_sem(PP_I(i),:), 'r');
    ylabel('spikes/s'); xlabel('Direction')
    ymax = nanmax(max(GratingDirRespV(PP_I(i),:))+max(GratingDirRespV_sem(PP_I(i),:)), max(GratingDirResp(PP_I(i),:))+max(GratingDirResp_sem(PP_I(i),:)));
    ymin = min(0,min(min(GratingDirRespV(PP_I(i),:))+min(GratingDirRespV_sem(PP_I(i),:)), min(GratingDirResp(PP_I(i),:))+min(GratingDirResp_sem(PP_I(i),:))));
    try
        ylim([ymin ymax]);
    end
    xlim([-45 360]);
    suptitle(['#' num2str(PP_I(i)) ...
              ', AUC=' num2str(AUC_shuffled(PP_I(i),1),2) ...
              ', MI=' num2str(DM(PP_I(i)),2) ...
              ', P' num2str(PertRespUnits_pos(PP_I(i))) ... 
              ', p_H_T = ' num2str(TuningPropsV.p_HT2(PP_I(i))<0.01,3) ...
              ', p_t_t = ' num2str(TuningPropsV.p_tt(PP_I(i))<0.01,3)])
    set(gca,'XTick',0:45:315,'box','off','TickDir','out')
    box off
    set(gca, 'XAxisLocation', 'top')
    
    %% plot perturbation responses
    PlotNR =[54 60; 42 48; 41 47; 40 46; 52 58; 64 70; 65 71; 66 72]; 
    UnitResponse = squeeze(UnitResponses(:,PP_I(i),:));
    for sp = 1 : size(PlotNR,1)
        subplot(12,6,PlotNR(sp,:))
        % visual
        % plot all trials
        Trial2plot = find(Trial_Angle==GratingDirections(sp));
        % plot mean
        [val FirstPart] = min(abs(TimeLine-Pert_lims(1)));
        shadedErrorBar(TimeLine(1:FirstPart-90),...
            mean(UnitResponse(Trial2plot,1:FirstPart-90)),...
            std(UnitResponse(Trial2plot,1:FirstPart-90))/sqrt(size(UnitResponse(Trial2plot,1:FirstPart-90),1)),...
            'lineprop',{'k-','markerfacecolor',[0.5 0.5 0.5]});
        hold on
        xlim([-0.8 max(TimeLine)-0.5]);
        set(gca,'box','off','TickDir','out')
        % perturbation
        Trial2plot_p = find(Trial_Angle==GratingDirections(sp) & Trial_Pert==1);
        Trial2plot_np= find(Trial_Angle==GratingDirections(sp) & Trial_Pert==0);
        % plot mean
        hold on
        shadedErrorBar(TimeLine(FirstPart-60:end),...
            mean(UnitResponse(Trial2plot_np,FirstPart-60:end)),...
            std(UnitResponse(Trial2plot_np,FirstPart-60:end))/sqrt(size(UnitResponse(Trial2plot_np,FirstPart-60:end),1)),...
            'lineprop',{'k-','markerfacecolor','k'});
        shadedErrorBar(TimeLine(FirstPart-60:end),...
            mean(UnitResponse(Trial2plot_p,FirstPart-60:end)),...
            std(UnitResponse(Trial2plot_p,FirstPart-60:end))/sqrt(size(UnitResponse(Trial2plot_p,FirstPart-60:end),1)),...
            'lineprop',{'r-','markerfacecolor','r'});
        
        Limit_Y_axis = max([ max(mean(UnitResponse(Trial2plot,1:FirstPart-90)))*1.1, ...
            max(mean(UnitResponse(Trial2plot_p,FirstPart-60:end)))*1.1, ...
            max(mean(UnitResponse(Trial2plot_np,FirstPart-60:end)))*1.1]);
        Trial_nr_p(i) = length(Trial2plot_p);
        Trial_nr_np(i) = length(Trial2plot_np);
        ylim([0 max(Limit_Y_axis)]);
        plot([0 0], [0 max(Limit_Y_axis)],'k:'); plot([max(TimeLine)-1 max(TimeLine)-1],[0 max(Limit_Y_axis)],'k:'); % trial start and stop
        plot([Pert_lims(1) Pert_lims(1)], [0,max(Limit_Y_axis)],'r:'); plot([Pert_lims(2) Pert_lims(2)], [0,max(Limit_Y_axis)],'r:'); % pert start and stop
        axis tight
        text(0, Limit_Y_axis*0.8, [num2str(GratingDirections(sp)) '^o'])
%         set(gca,'XTickLabel',[]);
%         %if i ~=5
%             set(gca,'YTickLabel',[]);
%             set(gca,'xcolor','none','ycolor','none')
        %end
        %legend({['pert. ON n=' num2str((Trial_nr_p(i)))],...
        %   ['pert. OFF n=' num2str((Trial_nr_np(i)))]});
    end
    
    %% plot mean response to all perturbation trials
     subplot(12,6,[53 59])
        hold on
        xlim([-0.8 max(TimeLine)-0.5]);
        set(gca,'TickDir','out')
        % perturbation
        Trial2plot_p = find(Trial_Pert==1);
        Trial2plot_np= find(Trial_Pert==0);
        % plot mean
        hold on
        shadedErrorBar(TimeLine,...
            mean(UnitResponse(Trial2plot_np,:)),...
            std(UnitResponse(Trial2plot_np,:))/sqrt(size(UnitResponse(Trial2plot_np,:),1)),...
            'lineprop',{'k-','markerfacecolor','k'});
        shadedErrorBar(TimeLine,...
            mean(UnitResponse(Trial2plot_p,:)),...
            std(UnitResponse(Trial2plot_p,:))/sqrt(size(UnitResponse(Trial2plot_p,:),1)),...
            'lineprop',{'r-','markerfacecolor','r'});
        
        Limit_Y_axis = max([ max(mean(UnitResponse(Trial2plot_np,:)))*1.1, ...
            max(mean(UnitResponse(Trial2plot_p,:)))*1.1]);
        Trial_nr_p(i) = length(Trial2plot_p);
        Trial_nr_np(i) = length(Trial2plot_np);
        ylim([0 max(Limit_Y_axis)]);
        plot([0 0], [0 max(Limit_Y_axis)],'k:'); plot([max(TimeLine)-1 max(TimeLine)-1],[0 max(Limit_Y_axis)],'k:'); % trial start and stop
        plot([Pert_lims(1) Pert_lims(1)], [0,max(Limit_Y_axis)],'r:'); plot([Pert_lims(2) Pert_lims(2)], [0,max(Limit_Y_axis)],'r:'); % pert start and stop
        axis tight
        text(0, Limit_Y_axis*0.8, 'All pert. trials')   
     
    %% save the figure
    if length(PP_I)==62
         saveas(gcf, ['E:\Dropbox\UCL\Tomaso-Aman\VisPerturbation\Fig_3\Acceleration\PertUnits_AccelResponses' filesep ...
                'AccelResponse_' ...
                'Rec#' num2str(Rec_OI) ...
                '_Unit#' num2str(PP_I(i))...
                'Pert.png']);
    elseif length(PP_I)==47
        saveas(gcf, ['E:\Dropbox\UCL\Tomaso-Aman\VisPerturbation\Fig_3\Acceleration\PertUnits_AccelResponses_neg' filesep ...
                'AccelResponse_' ...
                'Rec#' num2str(Rec_OI) ...
                '_Unit#' num2str(PP_I(i))...
                'Pert.png']);
    end
    close all
end



%% plot preferred TF VS MI of perturbation response (separate visually suppressed from excited units)

% additional selection criteria for reliable units
% apply 0.1 Hz threshold also for FR during acceleration stimulus
I_main = find(AM_UOI);
I_accel = find(AM_UOI_acc);
for i = 1:length(I_main)
    for j = 1:length(I_accel)
        if I_main(i) == I_accel(j)
            I_bi(i) = true;
        end
    end
end
% try to increase it to see what you get!
% manually discard those units?
I_bi = true(1,189);

% select all units responding to perturbation 
p_pert_th_01 = prctile(Sh_responses(:),99);
PertResp_units_01 = (AUC_shuffled(:,1)>p_pert_th_01);
p_pert_th_05 = prctile(Sh_responses(:),95);
PertResp_units_05 = (AUC_shuffled(:,1)>p_pert_th_05);

% select preferred TF frequency
Resp_mean_bi = Resp_mean(I_bi,:);
[Pref_TF TF_i] = max(Resp_mean_bi,[],2);

DM_bi = DM(I_bi);
PertResp_units_01 = PertResp_units_01(I_bi);
PertResp_units_05 = PertResp_units_05(I_bi);
figure
% p<0.01
subplot(2,2,1)
DM_2plot = DM_bi(PertResp_units_01); U_speed2plot = max(U_speed(TF_i(PertResp_units_01)),0);
plot(U_speed2plot(DM_2plot>0),DM_2plot(DM_2plot>0), '.r','MarkerSize',15)
hold on
plot(U_speed2plot(DM_2plot<0),DM_2plot(DM_2plot<0), '.k','MarkerSize',15)
ylabel('Perturbation MI_F_R')
ylim([-1 1])
set(gca,'box','off','TickDir','out','XTickLabel',[],'YTick',[-1 0 1])
grid on;
legend({['MI>0 n=' num2str(length(DM_2plot(DM_2plot>0)))],['MI<0 n=' num2str(length(DM_2plot(DM_2plot<0)))]})
subplot(2,2,1)
title('All pert. units (AUC - p<0.01)')
subplot(2,2,3)
histogram(U_speed2plot(DM_2plot>0),[0:0.5:9],'FaceColor',[1 0 0],'EdgeColor','none')
hold on
histogram(U_speed2plot(DM_2plot<0),[0:0.5:9],'FaceColor',[0 0 0],'EdgeColor','none')
xlabel('Preferred TF @ 0^o');
ylabel('Units');
set(gca,'box','off','TickDir','out')
subplot(2,2,3)
p = ranksum(U_speed2plot(DM_2plot>0),U_speed2plot(DM_2plot<0))
title(['ranksum test p=' num2str(p,2)]);
% p<0.05
subplot(2,2,2)
DM_2plot = DM_bi(PertResp_units_05); U_speed2plot = max(U_speed(TF_i(PertResp_units_05)),0);
semilogx(U_speed2plot(DM_2plot>0),DM_2plot(DM_2plot>0), '.r','MarkerSize',15)
hold on
semilogx(U_speed2plot(DM_2plot<0),DM_2plot(DM_2plot<0), '.k','MarkerSize',15)
ylabel('Perturbation MI_F_R')
ylim([-1 1])
set(gca,'box','off','TickDir','out','XTickLabel',[],'YTick',[-1 0 1])
grid on;
legend({['MI>0 n=' num2str(length(DM_2plot(DM_2plot>0)))],['MI<0 n=' num2str(length(DM_2plot(DM_2plot<0)))]})
subplot(2,2,2)
title('All pert. units (AUC - p<0.05)')
subplot(2,2,4)
histogram(U_speed2plot(DM_2plot>0),[0:0.5:9],'FaceColor',[1 0 0],'EdgeColor','none')
hold on
histogram(U_speed2plot(DM_2plot<0),[0:0.5:9],'FaceColor',[0 0 0],'EdgeColor','none')
xlabel('Preferred TF @ 0^o');
ylabel('Units');
set(gca,'box','off','TickDir','out')
subplot(2,2,4)
p = ranksum(U_speed2plot(DM_2plot>0),U_speed2plot(DM_2plot<0))
title(['ranksum test p=' num2str(p,2)]);


%% data from Allen Brain Institute - Data Portal
% http://observatory.brain-map.org/visualcoding
TF = [1 2 4 8 15];
Units_200_400 = [4609 1651 1225 880 613];
Units_550_570 = [52 44 64 66 35];

figure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,1,1)
I_bi = true(1,189);
% select all units responding to perturbation 
p_pert_th_05 = prctile(Sh_responses(:),95);
PertResp_units_05 = (AUC_shuffled(:,1)>p_pert_th_05);
% select preferred TF frequency
Resp_mean_bi = Resp_mean(I_bi,:);
[Pref_TF TF_i] = max(Resp_mean_bi,[],2);
DM_bi = DM(I_bi);
PertResp_units_01 = PertResp_units_01(I_bi);
PertResp_units_05 = PertResp_units_05(I_bi);
U_speed2plot = max(U_speed(TF_i(PertResp_units_01)),0);
DM_2plot = DM_bi(PertResp_units_05); U_speed2plot = max(U_speed(TF_i(PertResp_units_05)),0);
clear Units2Count
TFs = U_speed(U_speed>=0); 
Edges = [TFs-0.25; TFs(end)+0.25];
Edges = [-0.25 1.25 2.75 5.75 9.25]
Units2Count(:,1) = histcounts(U_speed2plot(DM_2plot>0),Edges); % positively modulated
Units2Count(:,2) = histcounts(U_speed2plot(DM_2plot<0),Edges); % negatively modulated
Units2Count(:,3) = histcounts(max(U_speed(TF_i(~PertResp_units_05)),0),Edges); % rest of units

xlabels = {'<=1', '1.5-2.5', '3-5.5', '6-9'}
H = bar(Units2Count/sum(Units2Count(:))*100,'stacked','EdgeColor','none')
set(gca,'XTick',1:length(Edges)-1,'XTickLabel',xlabels)
set(gca,'TickDir','out','box','off')
xlim([0.5 5.5])
myColors = [1 0 0; 0 0 1; 0 0 0];
for i = 1:3
    H(i).FaceColor = 'flat';
    H(i).CData = myColors(i,:);
end
legend({'MI>0 (n=62)','MI<0 (n=47)','rest (n=80)'})
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,1,2)
b = bar(Units_200_400/sum(Units_200_400)*100, 0.8,'FaceColor',[0 .5 .5],'EdgeColor','none')
ylabel('Units (%)')
set(gca,'XTickLabel',[],'TickDir','out','box','off')
title('Imaging depth: 200-400 \mum')
xlim([0.5 5.5])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,1,3)
bar(Units_550_570/sum(Units_550_570)*100, 0.8,'FaceColor',[0 .5 .5],'EdgeColor','none')
ylabel('Units (%)')
xlabel('Preferred temporal frequency (cycles/s)')
set(gca,'XTickLabel',TF,'TickDir','out','box','off')
title('Imaging depth: 550-570 \mum')
xlim([0.5 5.5])
log
