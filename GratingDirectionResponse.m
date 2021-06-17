%% Tomaso Muzzu - 21/10/2019 - UCL

%% response to grating direction %% CHECK CODE FOR BUGS!!! (NORMALIZATION)

% AM_Param :
% conditions = 1 --> nr of recording
% conditions = 2 --> grating direction [0:45:315]
% conditions = 3 --> 1 pert ON, 0 pert OFF
% conditions = 4 --> pert onset
% conditions = 5 --> pert offset

% find angle direction of each trial
Param = AM_Param(:,AM_UOI,2);
UnitResponses = AM_UnitResponses_smooth(:,AM_UOI,:)*60;
% for t = 1:size(UnitResponses,2)
%     OneUnitResponse = squeeze(UnitResponses(:,t,:));
%     UnitResponses_n(:,t,:) = OneUnitResponse./nanmax(OneUnitResponse,[],2);
% end
%UnitResponses = UnitResponses_n;
GratingDirections = unique(Param(~isnan(Param(:,1)),1));
clear Gr_Dir_2D t1 TrialsResp_OI
for i = 1:length(GratingDirections)
   Gr_Dir_2D(:,:,i) = Param==GratingDirections(i); 
   t1 = double(repmat(Gr_Dir_2D(:,:,i), 1, 1, size(UnitResponses,3))); % repeating the selection across time
   t1(t1(:)==0) = nan; 
   TrialsResp_OI{i} = t1.*UnitResponses;
end

%% compute significance of visual responses
% 1) compute MI for visual responses for different angles
ResponseTime = [-1 +4];
BonvisionFR = 60; %Hz
trialSide_samples = 60;
trialSide_seconds = 1; clear FR_vis MI_vis
for i = 1:length(TrialsResp_OI)
    for j = 1:length(TrialsResp_OI{i})
        SingleAngleMeanResponse = nanmean(squeeze(TrialsResp_OI{i}(~isnan(TrialsResp_OI{i}(:,j,1)),j,:)));
        FR_vis(j,i,1) = mean(SingleAngleMeanResponse(trialSide_samples:trialSide_samples+ResponseTime(2)*BonvisionFR)) - mean(SingleAngleMeanResponse(1:trialSide_samples));
    end
end
Param = AM_Param(:,AM_UOI,2);
RecordingsNR = unique(AM_Param(1,AM_UOI,1));
for rNR = 1:length(RecordingsNR)
    rNR_i(rNR) = min(find(AM_Param(1,AM_UOI,1)==RecordingsNR(rNR)))
end
tic
for sh = 2:1001 % shuffling directions for all trials
    Param_sh = Param;
    for rec = 1:length(rNR_i)-1
        Sh_Angle = Param(1:min(find(isnan(Param(:,rNR_i(end)))))-1,rNR_i(rec));
        Sh_Angle = Sh_Angle(randperm(length(Sh_Angle)));
        Param_sh(1:length(Sh_Angle),rNR_i(rec):rNR_i(rec+1)-1) = repmat(Sh_Angle, 1, rNR_i(rec+1)-rNR_i(rec));
    end
    Sh_Angle = Param(1:min(find(isnan(Param(:,rNR_i(end)))))-1,rNR_i(end));
    Sh_Angle = Sh_Angle(randperm(length(Sh_Angle)));
    Param_sh(1:length(Sh_Angle),rNR_i(end):size(Param_sh,2)) = repmat(Sh_Angle, 1, size(Param_sh,2)-rNR_i(end)+1);
   
   clear Gr_Dir_2D t1 TrialsResp_OI
   for i = 1:length(GratingDirections)
       Gr_Dir_2D(:,:,i) = Param_sh==GratingDirections(i);
       t1 = double(repmat(Gr_Dir_2D(:,:,i), 1, 1, size(UnitResponses,3))); % repeating the selection across time
       t1(t1(:)==0) = nan;
       TrialsResp_OI{i} = t1.*UnitResponses;
   end
   for i = 1:length(TrialsResp_OI)
       for j = 1:length(TrialsResp_OI{i})
           SingleAngleMeanResponse = nanmean(squeeze(TrialsResp_OI{i}(~isnan(TrialsResp_OI{i}(:,j,1)),j,:)));
           FR_vis(j,i,sh) = (mean(SingleAngleMeanResponse(trialSide_samples:trialSide_samples+ResponseTime(2)*BonvisionFR)) - mean(SingleAngleMeanResponse(1:trialSide_samples)));% / ...
               %(mean(SingleAngleMeanResponse(trialSide_samples:trialSide_samples+ResponseTime(2)*BonvisionFR)) + mean(SingleAngleMeanResponse(1:trialSide_samples)));
       end
   end
   sh
end
toc
for i = 1:size(FR_vis,1)
    for j = 1:size(FR_vis,2)
        FR_p(i,j) = sum(abs(FR_vis(i,j,1))>abs(squeeze(FR_vis(i,j,2:end))));
    end
end

% go through responses of single untis
ResponseTime = [-1 +4];
BonvisionFR = 60; %Hz
trialSide_samples = 60;
trialSide_seconds = 1;
clear GratingDirRespV
for k = 1:size(UnitResponses,2)
    clear UnitResp SingleUnitResponse_oneDir
    for i = 1:length(GratingDirections)
        SingleUnitResponse_oneDir = squeeze(TrialsResp_OI{i}(~isnan(TrialsResp_OI{i}(:,k,1)),k,1:trialSide_samples+ResponseTime(2)*BonvisionFR));
        %SingleUnitResponse_oneDir = SingleUnitResponse_oneDir./max(SingleUnitResponse_oneDir);
        BL = nanmean(SingleUnitResponse_oneDir(:,1:trialSide_samples),2);
        UnitResp(i) = nanmean(nanmean(SingleUnitResponse_oneDir(:,trialSide_samples+1:end)-BL));
    end
    UnitResp(isnan(UnitResp)) = 0;
    UnitResp = UnitResp+abs(min(UnitResp));
    GratingDirRespV(k,:) = UnitResp;
end

%% check whether the unit is negatively modulated by perturbation
% load info regarding modulation index of the visual stimulus
load('DM_visresp_shuffled.mat') 
% FR_UOI_VIS(1,p_i,:) % mean FR of first 4 seconds of stim
% FR_UOI_VIS(2,p_i,:) % mean FR of 1s before stim onset
% FR_UOI_VIS(3,p_i,:) % DM
% FR_UOI_VIS(4,p_i,:) % MI
% FR_UOI_VIS(5,p_i,:) % avg baseline FR
% select best units
SelectedCells = AM_UOI;
FR_UOI_VIS = FR_UOI_VIS(:,:,SelectedCells);
DM_vis = squeeze(FR_UOI_VIS(3,1,:));

% flip the tuning curve if modulation index DM is negative
for d = 1:size(GratingDirRespV,1)
    if DM_vis(d)<=0
       GratingDirRespV(d,:) = -(GratingDirRespV(d,:)-max(GratingDirRespV(d,:)));
    end
end

% Compute OSI and DSI for each unit 
clear TuningPropsV Pref_Dir OSI DSI L_ori L_dir
for k = 1:size(UnitResponses,2)
    [R_pref ind] = max(GratingDirRespV(k,:));
    R_null = GratingDirRespV(k,max(mod(ind+4,length(GratingDirections)),1));
    R_ortho1 = GratingDirRespV(k,max(mod(ind+2,length(GratingDirections)),1));
    R_ortho2 = GratingDirRespV(k,max(mod(ind+6,length(GratingDirections)),1));
    Pref_Dir(k,1) = GratingDirections(ind);% perferred direction; 
    % ORI = (R_pref+R_null)-(R_ortho1+R_ortho2)/(R_pref+R_null)+(R_ortho1+R_ortho2)
    OSI(k,1) = ((R_pref+R_null)-(R_ortho1+R_ortho2))/((R_pref+R_null)+(R_ortho1+R_ortho2)); 
    % DSI = (R_pref-R_null)/(R_pref+R_null)
    DSI(k,1) = (R_pref-R_null)/(R_pref+R_null);
    % circular variance or magnitude of the vector L = 1-CirVar, ORI space
    L_ori(k,1) = norm( GratingDirRespV(k,:)*exp(2i*degtorad(GratingDirections))/sum(GratingDirRespV(k,:)) ) ;
    P_ori(k,1) = rad2deg(angle(GratingDirRespV(k,:)*exp(2i*degtorad(GratingDirections))/sum(GratingDirRespV(k,:))))/2;
    if sign(P_ori(k,1))==-1
        P_ori(k,1) = P_ori(k,1)+180;
    end
    % circular variance or magnitude of the vector L = 1-CirVar, DIR space
    L_dir(k,1) = norm( GratingDirRespV(k,:)*exp(1i*degtorad(GratingDirections))/sum(GratingDirRespV(k,:))  );
    P_dir(k,1) = rad2deg(angle(GratingDirRespV(k,:)*exp(1i*degtorad(GratingDirections))/sum(GratingDirRespV(k,:))));
    if sign(P_dir(k,1))==-1
        P_dir(k,1) = 360+P_dir(k,1);
    end
end
TuningPropsV = table(Pref_Dir,OSI,DSI,L_ori,P_ori,L_dir,P_dir);

%% compute significance of response for Orientation tuning : Hotelling's T-squared test
clear GratingDirResp_trial
for k = 1:size(UnitResponses,2)
    clear UnitResp_all SingleUnitResponse_oneDir
    for i = 1:length(GratingDirections)
        SingleUnitResponse_oneDir = squeeze(TrialsResp_OI{i}(~isnan(TrialsResp_OI{i}(:,k,1)),k,1:trialSide_samples+ResponseTime(2)*BonvisionFR));
        %SingleUnitResponse_oneDir = SingleUnitResponse_oneDir./max(SingleUnitResponse_oneDir);
        BL = nanmean(SingleUnitResponse_oneDir(:,1:trialSide_samples),2);
        UnitResp_all{i} = nanmean(SingleUnitResponse_oneDir(:,trialSide_samples+1:end)-BL,2);
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
    p_HT2(i) = HotellingT2Test(Ori_vec,p_thres);
    TuningPropsV.p_HT2(i) = p_HT2(i);
end

%% compute significance of response for Direction tuning: direction dot product test
clear p_tt
for k = 1:size(GratingDirResp_trial,2)
    clear OriVec OriAxis Dir_vec ProjMagn
    % step 1 : calculate the average orientation vector onto the orientation
    % axis
    OriVec = (GratingDirRespV(k,:)*exp(2i*degtorad(GratingDirections)))/sum(GratingDirRespV(k,:));
    OriAxis = rad2deg(angle(GratingDirRespV(k,:)*exp(2i*degtorad(GratingDirections))))/2;
    if sign(OriAxis)==-1
        OriAxis = OriAxis+180;
    end
    OriAxis = deg2rad(OriAxis);
    % OriAxis = acos(real(OriVec)/imag(OriVec));
    % step 2 : calculate the magnitude of the projection of each direction
    % vector onto the orientation axis
    clear Dir_vec ProjMagn
    for j = 1:size(GratingDirResp_trial{k},2)
        Dir_vec(j,:) =  [real(GratingDirResp_trial{1,k}(:,j)'*exp(1i*degtorad(GratingDirections))) ...
            imag(GratingDirResp_trial{1,k}(:,j)'*exp(1i*degtorad(GratingDirections)))];
        ProjMagn(j) = dot(Dir_vec(j,:),[real(OriVec) imag(OriVec)]);
    end
    % step 3 : compute Student's T-test on the distribution of direction dot
    % products
    [h p_tt(k)] = ttest(ProjMagn); 
    TuningPropsV.p_tt(k) = p_tt(k);
end
p_thres = 0.001;
sum(p_HT2<p_thres)
sum(p_tt<p_thres)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look at tuning curves to vis. stim. onset
%% consider only units that perturbation responsive 
load('AUC_shuffled.mat')
Sh_responses = AUC_shuffled(:,2:end);
p_pert_th = prctile(Sh_responses(:),99);
Resp_PERT = AUC_shuffled(:,1)>p_pert_th;
Resp_ORI = TuningPropsV.p_HT2<p_thres & TuningPropsV.p_tt>p_thres; % units that are tuned for orientation but for direction
Resp_DIR = TuningPropsV.p_HT2>p_thres & TuningPropsV.p_tt<p_thres; % units that are tuned for direction but not orientation
Resp_BOTH = TuningPropsV.p_HT2<p_thres & TuningPropsV.p_tt<p_thres; % units that tuned for both 

Combo = Resp_PERT & Resp_ORI;
Combo = Resp_PERT & Resp_DIR;
Combo = Resp_PERT & Resp_BOTH;

GratingDirRespV_sel = GratingDirRespV(Combo,:);

plot_rows = 5; plot_cols = 10;
for u = 1:size(GratingDirRespV_sel,1)
    if mod(u,plot_rows*plot_cols)==1
        figure
        suptitle('visual stimulus response (norm.)')
    end
    nr_plot = mod(u,plot_rows*plot_cols);
    if nr_plot==0
        subplot(plot_rows,plot_cols,plot_rows*plot_cols)
    else
        subplot(plot_rows,plot_cols,nr_plot)
    end
    polarplot([deg2rad(GratingDirections); 0],[GratingDirRespV_sel(u,:) GratingDirRespV_sel(u,1)])
    thetaticklabels([])
end

Combo = Resp_PERT & Resp_ORI;
Combo = Resp_PERT & Resp_DIR;
Combo = Resp_PERT & Resp_BOTH;

GratingDirRespV_sel = GratingDirRespV(Combo,:);
GratingDirRespP_sel = GratingDirResp(Combo,:);
%%
Combo = Resp_PERT & Resp_ORI;
% UI_ori = [ 3 9 10 11 22 23 25 28 4 ]; FileSuffix = 'Ori';
Combo = Resp_PERT & Resp_DIR;
% UI_dir = [ 8 9 ]; FileSuffix = 'Dir';
Combo = Resp_PERT & Resp_BOTH;
% UI_both = [ 2 5 6 14 16 19 20 22 23 25 37 38 41 42 43 44 45 46 47 48 49 24 39 ]; FileSuffix = 'OriDir';
Combo = Resp_BOTH;
% everything --> apparently dir tuned #4
%% CIRCULAR PLOTS
Unit_Inds = find(Combo);
for UI = UI_ori %1:length(Unit_Inds)
    % UI = 3;
    UI_idx = Unit_Inds(UI);
    if DM(UI_idx)
        param_sel = AM_Param(:,SelectedCells,:);
        Pert_lims =[param_sel(min(find(param_sel(:,UI_idx,3)==1)),UI_idx,4) param_sel(min(find(param_sel(:,UI_idx,3)==1)),UI_idx,5)];
        TimeLine = linspace(-1,8.33,size(AM_UnitResponses_smooth,3));
        UnitResponse = squeeze(UnitResponses(:,UI_idx,:));
        Trial_Angle = squeeze(param_sel(:,UI_idx,2));
        Trial_Pert = squeeze(param_sel(:,UI_idx,3));
        figure
        set(gcf,'Position',[-1500 50 900 675])
        GratingDirectionsOrdered = [135 90 45 180 0 225 270 315];
        for i = 1 : length(GratingDirectionsOrdered)
            if i >= 5
                subplot(3,3,i+1)
            else
                subplot(3,3,i)
            end
            % visual
            % plot all trials
            Trial2plot = find(Trial_Angle==GratingDirectionsOrdered(i));
            %     for j = 1:length(Trial2plot)
            %         plot(TimeLine(1:FirstPart-90),UnitResponse(Trial2plot(j),TimeLine-Pert_lims(1)),'Color',[0.5 0.5 0.5]);
            %         hold on
            %     end
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
            Trial2plot_p = find(Trial_Angle==GratingDirectionsOrdered(i) & Trial_Pert==1);
            Trial2plot_np= find(Trial_Angle==GratingDirectionsOrdered(i) & Trial_Pert==0);
            % plot all trials
            %     for j = 1:length(Trial2plot)
            %         plot(TimeLine(FirstPart-60:end),UnitResponse(Trial2plot(j),FirstPart-60:end),'Color',[0.5 0.5 0.5]);
            %         hold on
            %     end
            % plot mean
            shadedErrorBar(TimeLine(FirstPart-60:end),...
                mean(UnitResponse(Trial2plot_p,FirstPart-60:end)),...
                std(UnitResponse(Trial2plot_p,FirstPart-60:end))/sqrt(size(UnitResponse(Trial2plot_p,FirstPart-60:end),1)),...
                'lineprop',{'r-','markerfacecolor','r'});
            hold on
            shadedErrorBar(TimeLine(FirstPart-60:end),... 
                mean(UnitResponse(Trial2plot_np,FirstPart-60:end)),...
                std(UnitResponse(Trial2plot_np,FirstPart-60:end))/sqrt(size(UnitResponse(Trial2plot_np,FirstPart-60:end),1)),...
                'lineprop',{'k-','markerfacecolor','k'});
            
            Limit_Y_axis(i) = max([ max(mean(UnitResponse(Trial2plot,1:FirstPart-90)))*1.1, ...
                max(mean(UnitResponse(Trial2plot_p,FirstPart-60:end)))*1.1, ...
                max(mean(UnitResponse(Trial2plot_np,FirstPart-60:end)))*1.1]);
            Trial_nr_p(i) = length(Trial2plot_p);
            Trial_nr_np(i) = length(Trial2plot_np);
        end
        for i = 1 : length(GratingDirectionsOrdered)
            if i >= 5
                subplot(3,3,i+1)
            else
                subplot(3,3,i)
            end
            ylim([0 max(Limit_Y_axis)]);
            plot([0 0], [0 max(Limit_Y_axis)],'k:'); plot([max(TimeLine)-1 max(TimeLine)-1],[0 max(Limit_Y_axis)],'k:'); % trial start and stop
            plot([Pert_lims(1) Pert_lims(1)], [0,max(Limit_Y_axis)],'r:'); plot([Pert_lims(2) Pert_lims(2)], [0,max(Limit_Y_axis)],'r:'); % pert start and stop
            legend({['pert. ON n=' num2str((Trial_nr_p(i)))],...
                ['pert. OFF n=' num2str((Trial_nr_np(i)))]});
        end
        subplot(3,3,5)
        polarplot([deg2rad(GratingDirections); 0],[GratingDirRespV(UI_idx,:) GratingDirRespV(UI_idx,1)],'k')
        % hold on
        % polarplot([deg2rad(GratingDirections); 0],[GratingDirResp(UI_idx,:) GratingDirResp(UI_idx,1)],'r')
        hold on
        polarplot([deg2rad(GratingDirections); 0],[GratingDirResp(UI_idx,:)+GratingDirResp_BL(UI_idx,:) GratingDirResp(UI_idx,1)+GratingDirResp_BL(UI_idx,1)],'r')
        suptitle(['Unit ' num2str(find(I==UI_idx),3) ', ' ...
            'AUC=' num2str(AUC_shuffled(UI_idx,1),2) ', ' ...
            'MI=' num2str(DM_sh(UI_idx,1),2) ', ' ...
            'Ori_a=' num2str(TuningPropsV.P_ori(UI_idx),3) ', ' ...
            'Ori_v='  num2str(TuningPropsV.L_ori(UI_idx),2) ', ' ...
            'Dir_a=' num2str(TuningPropsV.P_dir(UI_idx),3) ', '...
            'Dir_v=' num2str(TuningPropsV.L_dir(UI_idx),2) ]);
        %     saveas(gcf,['X:\DATA\PROJECTS\VisPerturbation\Figures\Fig_2\panels\ExampleUnits' filesep FileSuffix '_' num2str(UI) '.pdf']);
        %     saveas(gcf,['X:\DATA\PROJECTS\VisPerturbation\Figures\Fig_2\panels\ExampleUnits' filesep FileSuffix '_' num2str(UI) '.fig']);
        %
        % input('');
    end
end

%% LINEAR PLOTS
Unit_Inds = find(Combo);
for UI = 1:length(Unit_Inds)
    % UI = 3;
    UI_idx = Unit_Inds(UI);
    %if DM(UI_idx)<0
        param_sel = AM_Param(:,SelectedCells,:);
        Pert_lims =[param_sel(min(find(param_sel(:,UI_idx,3)==1)),UI_idx,4) param_sel(min(find(param_sel(:,UI_idx,3)==1)),UI_idx,5)];
        TimeLine = linspace(-1,8.33,size(AM_UnitResponses_smooth,3));
        UnitResponse = squeeze(UnitResponses(:,UI_idx,:));
        Trial_Angle = squeeze(param_sel(:,UI_idx,2));
        Trial_Pert = squeeze(param_sel(:,UI_idx,3));
        figure
        set(gcf,'Position',[-1500 50 900 675])
        GratingDirectionsOrdered = [135 90 45 180 0 225 270 315];
        for i = 1 : length(GratingDirectionsOrdered)
            if i >= 5
                subplot(3,3,i+1)
            else
                subplot(3,3,i)
            end
            % visual
            % plot all trials
            Trial2plot = find(Trial_Angle==GratingDirectionsOrdered(i));
            %     for j = 1:length(Trial2plot)
            %         plot(TimeLine(1:FirstPart-90),UnitResponse(Trial2plot(j),TimeLine-Pert_lims(1)),'Color',[0.5 0.5 0.5]);
            %         hold on
            %     end
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
            Trial2plot_p = find(Trial_Angle==GratingDirectionsOrdered(i) & Trial_Pert==1);
            Trial2plot_np= find(Trial_Angle==GratingDirectionsOrdered(i) & Trial_Pert==0);
            % plot all trials
            %     for j = 1:length(Trial2plot)
            %         plot(TimeLine(FirstPart-60:end),UnitResponse(Trial2plot(j),FirstPart-60:end),'Color',[0.5 0.5 0.5]);
            %         hold on
            %     end
            % plot mean
            shadedErrorBar(TimeLine(FirstPart-60:end),...
                mean(UnitResponse(Trial2plot_p,FirstPart-60:end)),...
                std(UnitResponse(Trial2plot_p,FirstPart-60:end))/sqrt(size(UnitResponse(Trial2plot_p,FirstPart-60:end),1)),...
                'lineprop',{'r-','markerfacecolor','r'});
            hold on
            shadedErrorBar(TimeLine(FirstPart-60:end),... 
                mean(UnitResponse(Trial2plot_np,FirstPart-60:end)),...
                std(UnitResponse(Trial2plot_np,FirstPart-60:end))/sqrt(size(UnitResponse(Trial2plot_np,FirstPart-60:end),1)),...
                'lineprop',{'k-','markerfacecolor','k'});
            
            Limit_Y_axis(i) = max([ max(mean(UnitResponse(Trial2plot,1:FirstPart-90)))*1.1, ...
                max(mean(UnitResponse(Trial2plot_p,FirstPart-60:end)))*1.1, ...
                max(mean(UnitResponse(Trial2plot_np,FirstPart-60:end)))*1.1]);
            Trial_nr_p(i) = length(Trial2plot_p);
            Trial_nr_np(i) = length(Trial2plot_np);
        end
        for i = 1 : length(GratingDirectionsOrdered)
            if i >= 5
                subplot(3,3,i+1)
            else
                subplot(3,3,i)
            end
            ylim([0 max(Limit_Y_axis)]);
            plot([0 0], [0 max(Limit_Y_axis)],'k:'); plot([max(TimeLine)-1 max(TimeLine)-1],[0 max(Limit_Y_axis)],'k:'); % trial start and stop
            plot([Pert_lims(1) Pert_lims(1)], [0,max(Limit_Y_axis)],'r:'); plot([Pert_lims(2) Pert_lims(2)], [0,max(Limit_Y_axis)],'r:'); % pert start and stop
            legend({['pert. ON n=' num2str((Trial_nr_p(i)))],...
                ['pert. OFF n=' num2str((Trial_nr_np(i)))]});
        end
        subplot(3,3,5)
        polarplot([deg2rad(GratingDirections); 0],[GratingDirRespV(UI_idx,:) GratingDirRespV(UI_idx,1)],'k')
        % hold on
        % polarplot([deg2rad(GratingDirections); 0],[GratingDirResp(UI_idx,:) GratingDirResp(UI_idx,1)],'r')
        hold on
        polarplot([deg2rad(GratingDirections); 0],[GratingDirResp(UI_idx,:)+GratingDirResp_BL(UI_idx,:) GratingDirResp(UI_idx,1)+GratingDirResp_BL(UI_idx,1)],'r')
        suptitle(['Unit ' num2str(find(I==UI_idx),3) ', ' ...
            'AUC=' num2str(AUC_shuffled(UI_idx,1),2) ', ' ...
            'MI=' num2str(DM_sh(UI_idx,1),2) ', ' ...
            'Ori_a=' num2str(TuningPropsV.P_ori(UI_idx),3) ', ' ...
            'Ori_v='  num2str(TuningPropsV.L_ori(UI_idx),2) ', ' ...
            'Dir_a=' num2str(TuningPropsV.P_dir(UI_idx),3) ', '...
            'Dir_v=' num2str(TuningPropsV.L_dir(UI_idx),2) ]);
        %     saveas(gcf,['X:\DATA\PROJECTS\VisPerturbation\Figures\Fig_2\panels\ExampleUnits' filesep FileSuffix '_' num2str(UI) '.pdf']);
        %     saveas(gcf,['X:\DATA\PROJECTS\VisPerturbation\Figures\Fig_2\panels\ExampleUnits' filesep FileSuffix '_' num2str(UI) '.fig']);
        %
        % input('');
    %end
end


%% sort the cells by their response to perturbation (reference for heat map)
[B,I,PertTrial_Ex] = PlotSortingOrder(ProjectData,SelectedResponses,AM_Param,AM_Speed,Units_Sel,0);
% I = descending order taking into account:[pert response; 2s vis stim respo; post-pert reponse]
I = I(~isnan(B(:,1)));
B = reshape(B(~isnan(B(:))),length(B(~isnan(B(:))))/3,3);

%%
figure
for r = 1:size(UnitResponse(Trial2plot,:),1)
    plot((UnitResponse(r,:)),'b')
    hold on
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% plot pref Ori visual VS pref Ori perturbation
figure
title('All units - Orientation')
plot(TuningPropsV.P_ori,TuningProps.P_ori,'.')
%load('AUC_shuffled.mat')
%Sh_responses = AUC_shuffled(:,2:end);
p_pert_th = prctile(Sh_responses(:),99);
PertResp_units = (AUC_shuffled(:,1)>p_pert_th);
hold on
plot(TuningPropsV.P_ori(PertResp_units),TuningProps.P_ori(PertResp_units),'or')
xlabel('Visual pref. orientation'); ylabel('Perturbation pref. orientation'); 
legend({'All units', 'Pert. responsive units'})

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
figure
p_thres = 0.001; PertResp_units = (AUC_shuffled(:,1)>p_pert_th);
set(gcf,'Position',[100 100 650 600])
subplot(4,4,[9 10 13 14])
title('Ori. tuned units (vis.)')
h1 = plot(TuningPropsV.P_ori(TuningPropsV.p_HT2<=p_thres),TuningProps.P_ori(TuningPropsV.p_HT2<=p_thres),'.','Color',[0.9 0.9 0.9],'MarkerSize',12)
hold on
plot(TuningPropsV.P_ori(PertResp_units & TuningPropsV.p_HT2<=p_thres),...
     TuningProps.P_ori(PertResp_units & TuningPropsV.p_HT2<=p_thres),'o','MarkerEdgeColor',[229 44 37]/255,'LineWidth',3)
 
scatterDiagHist(TuningPropsV.P_ori(PertResp_units & TuningPropsV.p_HT2<=p_thres),...
     TuningProps.P_ori(PertResp_units & TuningPropsV.p_HT2<=p_thres)) 
 
xlabel('Visual preferred orientation'); ylabel('Perturbation max resp. orientation'); 
legend({['Ori. tuned vis. units n=' num2str(length(TuningPropsV.P_ori(TuningPropsV.p_HT2<=p_thres)))], ...
        ['Pert. responsive units n=' num2str(length(TuningPropsV.P_ori(PertResp_units & TuningPropsV.p_HT2<=p_thres)))]})
set(gca,'YAxisLocation', 'left','XTick',0:45:180,'YTick',0:45:180,'TickDir','out'); box off;
%set(gca, 'XDir','reverse')
xlim([0 180]); ylim([0 180])

subplot(4,4,[1 2 5 6])
histogram(TuningPropsV.P_ori(TuningPropsV.p_HT2<=p_thres),0:10:180,'EdgeColor','none','FaceColor','k')
hold on
histogram(TuningPropsV.P_ori(PertResp_units & TuningPropsV.p_HT2<=p_thres),0:10:180,'EdgeColor','none','FaceColor',[229 44 37]/255)
set(gca,'YAxisLocation', 'left','XTick',0:45:180,'XTickLabel',[],'TickDir','out'); box off;
%set(gca, 'XDir','reverse')
ylabel('Units')
xlim([0 180])

subplot(4,4,[11 12 15 16])
histogram(TuningProps.P_ori(TuningPropsV.p_HT2<=p_thres),0:10:180,'EdgeColor','none','FaceColor','k')
hold on
histogram(TuningProps.P_ori(PertResp_units & TuningPropsV.p_HT2<=p_thres),0:10:180,'EdgeColor','none','FaceColor',[229 44 37]/255)
xlim([0 180])
set(gca,'XTick',0:45:180,'XTickLabel',[],'TickDir','out','XAxisLocation', 'bottom'); box off;
ylabel('Units');
set(gca, 'XDir','reverse')
view(90,90)


subplot(4,4,[3 4 7 8])
scatterDiagHist(TuningPropsV.P_ori(PertResp_units & TuningPropsV.p_HT2<=p_thres),...
     TuningProps.P_ori(PertResp_units & TuningPropsV.p_HT2<=p_thres))  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% plot pref Dir visual VS pref Dir perturbation
figure
plot(TuningPropsV.P_dir,TuningProps.P_dir,'.')
%load('AUC_shuffled.mat')
%Sh_responses = AUC_shuffled(:,2:end);
% p_pert_th = prctile(Sh_responses(:),99);
% PertResp_units = (AUC_shuffled(:,1)>p_pert_th);
hold on
plot(TuningPropsV.P_dir(PertResp_units),TuningProps.P_dir(PertResp_units),'or')
xlabel('Visual pref. orientation'); ylabel('Perturbation max resp. orientation'); 
legend({'All units', 'Pert. responsive units'})
title('All units - Direction')

figure
plot(TuningPropsV.P_dir(TuningPropsV.p_tt<=p_thres),TuningProps.P_dir(TuningPropsV.p_tt<=p_thres),'.')
hold on
plot(TuningPropsV.P_dir(PertResp_units & TuningPropsV.p_tt<=p_thres),...
     TuningProps.P_dir(PertResp_units & TuningPropsV.p_tt<=p_thres),'or')
xlabel('Visual pref. orientation'); ylabel('Perturbation pref. orientation'); 
legend({'Dir. tuned vis. units', 'Pert. responsive units'})
%title('Dir. tuned units (vis.)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot pref Dir visual VS pref Dir perturbation
%%
figure
set(gcf,'Position',[100 100 650 600])
subplot(4,4,[9 10 13 14])
plot(TuningPropsV.P_dir(TuningPropsV.p_tt<=p_thres),TuningProps.P_dir(TuningPropsV.p_tt<=p_thres),'.','Color',[0.9 0.9 0.9],'MarkerSize',12)
hold on
plot(TuningPropsV.P_dir(PertResp_units & TuningPropsV.p_tt<=p_thres),...
     TuningProps.P_dir(PertResp_units & TuningPropsV.p_tt<=p_thres),'o','MarkerEdgeColor',[229 44 37]/255,'LineWidth',3)
% figure
scatterDiagHist(TuningPropsV.P_dir(PertResp_units & TuningPropsV.p_tt<=p_thres),...
     TuningProps.P_dir(PertResp_units & TuningPropsV.p_tt<=p_thres))
xlabel('Visual preferred direction'); ylabel('Perturbation max resp. direction'); 
legend({['Dir. tuned vis. units n=' num2str(length(TuningPropsV.P_dir(TuningPropsV.p_tt<=p_thres)))],...
        ['Pert. responsive units n=' num2str(length(TuningPropsV.P_dir(PertResp_units & TuningPropsV.p_tt<=p_thres)))]})
set(gca,'XTick',0:90:360,'YTick',0:90:360,'TickDir','out'); box off;
xlim([0 360]); ylim([0 360])

subplot(4,4,[1 2 5 6])
histogram(TuningPropsV.P_dir(TuningPropsV.p_tt<=p_thres),0:20:360,'EdgeColor','none','FaceColor','k')
hold on
histogram(TuningPropsV.P_dir(PertResp_units & TuningPropsV.p_tt<=p_thres),0:20:360,'EdgeColor','none','FaceColor',[229 44 37]/255)
set(gca,'XTick',0:90:360,'XTickLabel',[],'TickDir','out'); box off;
ylabel('Units')
xlim([0 360])

subplot(4,4,[11 12 15 16])
histogram(TuningProps.P_dir(TuningPropsV.p_tt<=p_thres),0:20:360,'EdgeColor','none','FaceColor','k')
hold on
histogram(TuningProps.P_dir(PertResp_units & TuningPropsV.p_tt<=p_thres),0:20:360,'EdgeColor','none','FaceColor',[229 44 37]/255)
hold on
xlim([0 360])
set(gca,'XTick',0:90:360,'XTickLabel',[],'TickDir','out','XAxisLocation', 'bottom'); box off;
ylabel('Units');
set(gca, 'XDir','reverse')
view(90,90)

subplot(4,4,[3 4 7 8])
scatterDiagHist(TuningPropsV.P_dir(PertResp_units & TuningPropsV.p_tt<=p_thres),...
     TuningProps.P_dir(PertResp_units & TuningPropsV.p_tt<=p_thres))
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
