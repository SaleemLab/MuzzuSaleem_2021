%% Tomaso Muzzu - UCL - 30/04/2020  - Eye responses

%% load data
if ~exist('ProjectData','var')
    [ProjectData AM_UnitResponses AM_Param AM_Speed AM_UOI SelectedResponses AM_UnitResponses_smooth AM_EyeTracking] = LoadDataALL;
end

%% first 7 animalsz
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
    % select only perturbation responsive units
    load('AUC_shuffled_CTRL_1.mat')
    Sh_responses = AUC_shuffled(:,2:end);
    p_pert_th = prctile(Sh_responses(:),99);
    PertResp_units = (AUC_shuffled(:,1)>p_pert_th);
    % select only pos. modulated perturbation responsive units
    load('DM_pert_shuffled_CTRL.mat')
    DM = DM_sh(:,1);
    DM_sign_i(:,1) = DM>0;
    DM_sign_i(:,2) = DM<=0;
    % select only pos. modulated perturbation responsive units
    PertRespUnits_pos = PertResp_units & DM_sign_i(:,1);
    PertRespUnits_neg = PertResp_units & DM_sign_i(:,2);
end

%%
AM_EyeSel = AM_EyeTracking(:,AM_UOI,1:559,6); % pupil size
AM_EyeSel = AM_EyeTracking(:,AM_UOI,1:559,1:2); % x and y centre
AM_EyeSel = sqrt( (diff(AM_EyeSel(:,:,:,1),1,3).^2) + (diff(AM_EyeSel(:,:,:,2),1,3).^2) ); % eye movements

UnitResponses = AM_UnitResponses_smooth(:,AM_UOI,:);
clear UnitResponses_n
for t = 1:size(UnitResponses,2)
    OneUnitResponse = squeeze(UnitResponses(:,t,:));
    UnitResponses_n(:,t,:) = OneUnitResponse./nanmax(OneUnitResponse,[],2);
end
UnitResponses = UnitResponses_n;
% normalise eye movements
AM_EyeSel_n = AM_EyeSel;
for t = 1:length(AM_EyeSel_n)
    OneAM_EyeSel_n = squeeze(AM_EyeSel_n(:,t,:));
    AM_EyeSel_n(:,t,:) = OneAM_EyeSel_n./nanmax(OneAM_EyeSel_n,[],2);
end

% look at the speed profiles of each recording
AM_ParamSel = AM_Param(:,AM_UOI,:);
for i = 1:length(unique(AM_ParamSel(1,:,1)))
    AllRec_idx(i) = min(find(AM_ParamSel(1,:,1)==i));
end
AM_EyeSel_nSel = AM_EyeSel_n(:,AllRec_idx,:);
% find beginning and end of perturbation trials
% scan every recording and look at the mean speed responses during
% perturbation
Rec_OI = unique(AM_ParamSel(1,PertRespUnits_pos,1)); % recordings of interest

for i = 1:length(Rec_OI)
   RecIndex(i) = min(find((AM_ParamSel(1,:,1)==Rec_OI(i))))
end
clear MeanRecSpeed SemRecSpeed
for i = 1:length(RecIndex)
   % get perturbation trials
   PertTrials_1 = AM_ParamSel(:,RecIndex(i),3)==1; 
   PertTrials_0 = AM_ParamSel(:,RecIndex(i),3)==0;
    
   MeanRecEye(1,:,i) = nanmean(AM_EyeSel(PertTrials_1,RecIndex(i),:));
   MeanRecEye(2,:,i) = nanmean(AM_EyeSel(PertTrials_0,RecIndex(i),:));
   SemRecEye(1,:,i) = nanstd(AM_EyeSel(PertTrials_0,RecIndex(i),:))/sqrt(sum(PertTrials_1));
   SemRecEye(2,:,i) = nanstd(AM_EyeSel(PertTrials_0,RecIndex(i),:))/sqrt(sum(PertTrials_0));   
end

BonvisionFR = 60; %Hz  
trialSide_samples = 60; trialSide_seconds = 1; 
PertStartEnds = [round(nanmin(AM_ParamSel(:,RecIndex,4),[],1)*BonvisionFR+trialSide_samples) ;...
                round(nanmax(AM_ParamSel(:,RecIndex,5),[],1)*BonvisionFR+trialSide_samples)];
% extract predictors from speed signals
po_i = [min(PertStartEnds(1,:)) max(PertStartEnds(2,:))];
clear DiffResponsePercent_Pert
for i = 1:size(AM_EyeSel_nSel,2)
    
    UnitResponse = squeeze(AM_EyeSel_nSel(:,i,:));
    
    DiffResponsePercent_Pert(:,i,1) = (nanmean(UnitResponse(:,po_i(1):po_i(2)),2)./nanmean(UnitResponse(:,po_i(1)-trialSide_samples:po_i(1)),2)); % mean(pert)/mean(prepert)
    DiffResponsePercent_Pert(:,i,2) = (nanmean(UnitResponse(:,po_i(1):po_i(2)),2)-nanmean(UnitResponse(:,po_i(1)-trialSide_samples:po_i(1)),2)); % mean(pert)-mean(prepert)
    DiffResponsePercent_Pert(:,i,3) = (nansum(UnitResponse(:,po_i(1):po_i(2)),2)); % integral of mean(pert)
    DiffResponsePercent_Pert(:,i,4) = (nansum(UnitResponse(:,po_i(1)-trialSide_samples:po_i(2)),2)); % integral of mean(pert) and preceding 1s
    DiffResponsePercent_Pert(:,i,5) = (nanmean(UnitResponse(:,po_i(1):po_i(2)-trialSide_samples/2),2)-nanmean(UnitResponse(:,po_i(1)-trialSide_samples:po_i(1)),2))./...
                                      (nanmean(UnitResponse(:,po_i(1):po_i(2)-trialSide_samples/2),2)+nanmean(UnitResponse(:,po_i(1)-trialSide_samples:po_i(1)),2)); % mean(pert)-mean(prepert)/mean(pert+max(prepert) a.k.a depth of modulation
    DiffResponsePercent_Pert(:,i,6) = nanmean(UnitResponse(:,po_i(1):po_i(2)-trialSide_samples/2),2)./...
                                      (nanmean(UnitResponse(:,po_i(1):po_i(2)-trialSide_samples/2),2)+nanmean(UnitResponse(:,po_i(1)-trialSide_samples:po_i(1)),2)); % mean(pert)/mean(pert)+mean(prepert) a.k.a modulation index    
end

%% step 2: use these values with the actual response to fit a logistic model
Param = squeeze(AM_Param(:,AllRec_idx,3)); % save trial type 
clear SingleUnitResponse AUC_Eye
plotted=0;
for i = 1:size(DiffResponsePercent_Pert,2) 
    % perturbation period
    SingleUnitResponse = squeeze(DiffResponsePercent_Pert(:,i,:));
    % select trial responses
    Responses = Param(:,i);
    % scrub responses
    Responses = logical(Responses(~isnan(Responses)));
    % scrub predictor
    SingleUnitResponse = SingleUnitResponse(1:length(Responses),:);
    SingleUnitResponse(SingleUnitResponse==Inf,:) = 0;
    SingleUnitResponse(isnan(SingleUnitResponse(:))) = 0;
    % apply logistic model
    mdl = fitglm(SingleUnitResponse,Responses,'Distribution','binomial','Link','logit');
    scores = mdl.Fitted.Probability;
    % step 3: compute the ROC and evaluate the AUC for each unit
    [X,Y,T,AUC1] = perfcurve(Responses,scores,1);
    AUC_Eye(i,1) = AUC1;
end
tic
clear AUC_shuffled_Eye
parfor i = 1:size(DiffResponsePercent_Pert,2)
    SingleUnitResponse = squeeze(DiffResponsePercent_Pert(:,i,:));
    % select trial responses
    Responses = Param(:,i);
    % scrub responses
    Responses = logical(Responses(~isnan(Responses)));
    % scrub predictor
    SingleUnitResponse = SingleUnitResponse(1:length(Responses),:);
    SingleUnitResponse(SingleUnitResponse==Inf,:) = 0;
    SingleUnitResponse(isnan(SingleUnitResponse(:))) = 0;
    
    for sh_i = 1:1000
        % shuffle perturbation trals
        Responses_sh = Responses(randperm(length(Responses)));
        
        % apply logistic model
        mdl = fitglm(SingleUnitResponse,Responses_sh,'Distribution','binomial','Link','logit');
        scores = mdl.Fitted.Probability;
        % step 3: compute the ROC and evaluate the AUC for each unit
        [X,Y,T,AUC1] = perfcurve(Responses,scores,1);
        AUC_shuffled_Eye(i,sh_i) = AUC1;
    end
    i
end
toc

AUC_shuffled_Eye = cat(2,AUC_Eye,AUC_shuffled_Eye);
% save data from shuffling
save('AUC_shuffled_Eye.mat','AUC_shuffled_Eye')


%% show AUC of shuffled and original data
RecIndex = unique(AM_ParamSel(1,PertRespUnits_pos,1));
AUC_2_plot = AUC_shuffled_Eye(RecIndex,2:end);
AUC_Eye =  AUC_shuffled_Eye(RecIndex,1);
p_sh = prctile(AUC_2_plot(:),99);

figure
h = histogram(AUC_2_plot(:),[0:0.02:1],'FaceColor',[0.5 0.5 0.5],'Normalization','probability','EdgeColor','none')
hold on
h1 = histogram(AUC_Eye,[0:0.02:1],'FaceColor',[1 0 0],'Normalization','probability','EdgeColor','none')
%plot([p_sh(1) p_sh(1)], [0 max(h.BinCounts/length(AUC_shuffled_Sp(:)))*1.1],'k') 
plot([p_sh p_sh], [0 max(h1.BinCounts/length(AUC_Eye(AM_ParamSel(1,RecIndex,1))))*1.1],'k--','LineWidth',3)
xlabel('AUC of shuffled data'); ylabel('probability')
ll = legend({['AUC of shuffled data n=' num2str(size(AUC_2_plot,1)) 'k'],['AUC of sessions with resp units n=' num2str(size(AUC_2_plot,1))], '99% shuffled data'})
ll.Color = 'none'; ll.EdgeColor = 'none';
box off
set(gca,'FontSize',13,'TickDir','out');

AUC_Eye = AUC_shuffled_Eye(:,1);
AUC_Eye>p_sh

% how many units come from these 6 recordings?
R_psh = RecIndex((AUC_Eye(RecIndex)>p_sh));
for i =1:length(R_psh)
    nrUnitsXrec(i) = sum(AM_ParamSel(1,PertRespUnits_pos,1)==R_psh(i))
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% look at PSTH of eye movements at perturbation onsets
AM_ParamSel = AM_Param(:,AM_UOI,:);
% AM_EyeSel = AM_EyeTracking(:,AM_UOI,1:559,6);
AM_SpeedSel = AM_EyeSel_n;
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
clear Trials
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
suptitle('Eye movementes in sessions contataining pert. resp. units')
for i=1:length(RecIndex)
    if length(RecIndex)<10
        subplot(1,length(RecIndex),i)
    else
        subplot(4,6,i)
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
        MeanRecSpeed(2,PertStartEnds(1,i)-60:PertStartEnds(2,i)+60,i),...
        SemRecSpeed(2,PertStartEnds(1,i)-60:PertStartEnds(2,i)+60,i),...
        'lineprop',{'k-','markerfacecolor',[0.5 0.5 0.5]});
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
    xlabel('seconds'); ylabel('A.U.');
    box off
    axis tight
    if i == 1
        ll = legend({'All trials', 'Only trials with Pert.'},'fontsize',10);
        ll.Color = 'none'; ll.EdgeColor = 'none';
    end
    set(gca,'TickDir','out')
    %ylim([max(ylim_mm(1)-5,0) ylim_mm(2)+5]);
    ylim([0 1])
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% safety check to evaluate whether eye pupil correlate with speed
f = figure
set(f,'renderer','painters');
for i = 1:size(ProjectData,1)
    if i>1
        subplot(10,3,[i*3-2 i*3-2+1])
    else
        subplot(10,3,[1 2])
    end
    Time = ProjectData.Session_data{i,1}.Time{1,1}(1:length(ProjectData.EyeTracking{i,1}.Area));
    Speed = ProjectData.Session_data{i,1}.MouseSpeed{1,1}(1:length(ProjectData.EyeTracking{i,1}.Area));
    % smooth speed signal
    Width = 18; % 6 samples = 100 milliseconds
    Sigma = Width/3; % standard deviation in number of samples (converted from time in seconds)
    x_g = linspace(-Width/2, Width/2, Width);
    gaussFilter = exp(-x_g.^2/(2*Sigma^2));
    gaussFilter_ = gaussFilter / sum (gaussFilter); % normalize
    Speed_smooth = conv(squeeze(Speed), gaussFilter_, 'same');
    Eye_mov = ProjectData.EyeTracking{i,1}.Area;
%     figure
    SpeedNorm = Speed_smooth/max(Speed_smooth);
    %SpeedNorm = (SpeedNorm-mean(SpeedNorm))/std(SpeedNorm);
    plot(Time,smooth(SpeedNorm,5),'k')
    hold on
    EyeNorm = Eye_mov/max(Eye_mov);
    % EyeNorm = (EyeNorm-mean(EyeNorm))/std(EyeNorm);
    if i >= 6 
        hold off
        plot(Time(Time>200),smooth(SpeedNorm(Time>200),5),'k')
        hold on
        plot(Time(Time>200),EyeNorm(Time>200),'r')
    else
        plot(Time, smooth(EyeNorm,5),'r')
    end
    if i == size(ProjectData,1)
        set(gca,'TickDir','out','box','off')
        xlabel('Seconds')
    else
        set(gca,'TickDir','out','box','off','XTickLabel',[])
    end
    ylabel('Norm. units')
    if i == 1
        legend({'Speed','Pupil size'})
    end
    axis tight
    if i > 6
        [PearsonCoeff p] = corrcoef(smooth(SpeedNorm(~isnan(EyeNorm)),5),EyeNorm(~isnan(EyeNorm)));
    elseif i == 6
        [PearsonCoeff p] = corrcoef(smooth(SpeedNorm(Time>200),5),EyeNorm(Time>200));
        xlim([0 max(Time)]);
    else
        [PearsonCoeff p] = corrcoef(smooth(SpeedNorm,5),smooth(EyeNorm,5));
    end
    title(['Rec ' num2str(i) ', R=' num2str(PearsonCoeff(1,2),2) ', p=' num2str(p(1,2),2)]);
    
    % make tuning curve of pupil size wrt speed
    if i>1
        subplot(10,3,i*3)
    else
        subplot(10,3,3)
    end
    
    if i == 6
        Speed_smooth = Speed_smooth(Time>200);
        [SpeedNorm_sorted idx] = sort(Speed_smooth);
        EyeNorm_sorted = Eye_mov(Time>200);
        EyeNorm_sorted = EyeNorm_sorted(idx);
    else
        [SpeedNorm_sorted idx] = sort(Speed_smooth);
        EyeNorm_sorted = ProjectData.EyeTracking{i,1}.Area(idx);
    end
    clear TC_speed TC_eye
    TC_speed(1,1) = max(0, mean(Speed_smooth_sorted(SpeedNorm_sorted<=1))); % save value for first x data point
    TC_speed(1,2) = std(Speed_smooth_sorted(SpeedNorm_sorted<=1));% /sqrt(length(Speed_smooth_sorted(Speed_smooth_sorted<=1)));
    TC_eye(1,1) = mean(EyeNorm_sorted(SpeedNorm_sorted<=1)); % save value for first y data point
    TC_eye(1,2) = std(EyeNorm_sorted(SpeedNorm_sorted<=1)); %/sqrt(length(EyeNorm_sorted(Speed_smooth_sorted<=1)));
    idx_m = sum(SpeedNorm_sorted<=1)+1;
    spdp = 500;
    for j = 1:ceil((length(SpeedNorm_sorted)-idx_m)/spdp)
        TC_speed(j+1,1) = mean(SpeedNorm_sorted(idx_m+(j-1)*spdp:min(idx_m+j*spdp,length(SpeedNorm_sorted))));
        TC_speed(j+1,2) = std(SpeedNorm_sorted(idx_m+(j-1)*spdp:min(idx_m+j*spdp,length(SpeedNorm_sorted))));
        TC_eye(j+1,1) = mean(EyeNorm_sorted(idx_m+(j-1)*spdp:min(idx_m+j*spdp,length(EyeNorm_sorted))));
        TC_eye(j+1,2) = std(EyeNorm_sorted(idx_m+(j-1)*spdp:min(idx_m+j*spdp,length(EyeNorm_sorted))));
    end
    errorbar(TC_speed(:,1),TC_eye(:,1),TC_eye(:,2),TC_eye(:,2),TC_speed(:,2),TC_speed(:,2),'r:o')
    if mod(i,2)~=0
        ylabel('Pupil size (A.U.)');
    end
    set(gca,'TickDir','out','box','off')
    if i == 10
        xlabel('Speed (cm/s)')
    end
    i
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%


















