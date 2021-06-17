%% Tomaso Muzzu - UCL - 26/11/2019

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
load('AUC_shuffled.mat'); % load info of pert. responsive units
RespUnits = AUC_shuffled(:,1)>(prctile(reshape(AUC_shuffled(:,2:end),1,size(AUC_shuffled(:,2:end),1)*size(AUC_shuffled(:,2:end),2)),99));
Rec_OI = unique(AM_ParamSel(1,RespUnits,1)); % recordings of interest

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

%% plot the mean curves for each recordings
figure 
TimeLine = linspace(-trialSide_seconds,size(MeanRecSpeed,2)/BonvisionFR-trialSide_seconds,size(MeanRecSpeed,2));
suptitle('Avg. speed responses in recordings contataining pert. responsive units')
for i=1:length(RecIndex)
    subplot(4,7,i)
    shadedErrorBar( TimeLine,...
        MeanRecSpeed(1,:,i),...
        SemRecSpeed(1,:,i),...
        'lineprop',{'r-','markerfacecolor',[0.5 0.5 0.5]});
    hold on
    shadedErrorBar( TimeLine,...
        MeanRecSpeed(2,:,i),...
        SemRecSpeed(2,:,i),...
        'lineprop',{'k-','markerfacecolor',[0.5 0.5 0.5]}); 
    hold on
    plot([TimeLine(trialSide_samples) TimeLine(trialSide_samples)], [0 max(max(MeanRecSpeed(1,:,i)),max(MeanRecSpeed(2,:,i)))],'k:')
    plot([TimeLine(end-trialSide_samples) TimeLine(end-trialSide_samples)], [0 max(max(MeanRecSpeed(1,:,i)),max(MeanRecSpeed(2,:,i)))],'k:')
    plot([TimeLine(PertStartEnds(1,i)) TimeLine(PertStartEnds(1,i))], [0 max(max(MeanRecSpeed(1,:,i)),max(MeanRecSpeed(2,:,i)))],'r:')
    plot([TimeLine(PertStartEnds(2,i)) TimeLine(PertStartEnds(2,i))], [0 max(max(MeanRecSpeed(1,:,i)),max(MeanRecSpeed(2,:,i)))],'r:')
    xlabel('seconds'); ylabel('cm/s');
    box off
    axis tight
    if i == 1
        ll = legend({'Perturbation ON', 'Perturbation OFF'},'fontsize',10);
        ll.Color = 'none'; ll.EdgeColor = 'none';
    end
    set(gca,'TickDir','out')
end

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
figure
suptitle('Speed responses in recordings contataining pert. responsive units')
for i=1:length(RecIndex)
    
    subplot(4,7,i)
    imagesc(TimeLine,1:size(Trials{i},1),Trials{i})
    hold on
    plot([min(TimeLine) min(TimeLine)], [1 sum(AM_ParamSel(:,RecIndex(i),3)==1)], 'r', 'LineWidth',4)
    hold on
    plot([min(TimeLine) min(TimeLine)], [sum(AM_ParamSel(:,RecIndex(i),3)==1) ...
        sum(AM_ParamSel(:,RecIndex(i),3)==1)+sum(AM_ParamSel(:,RecIndex(i),3)==0)], 'k', 'LineWidth',4)
    hold on
    plot([TimeLine(PertStartEnds(1,i)) TimeLine(PertStartEnds(2,i)) TimeLine(PertStartEnds(2,i)) TimeLine(PertStartEnds(1,i)) TimeLine(PertStartEnds(1,i))], ...
            [0 0 sum(AM_ParamSel(:,RecIndex(i),3)==1) sum(AM_ParamSel(:,RecIndex(i),3)==1) 0],'r:')
    axis tight
    xlabel('seconds'); ylabel('trials')
    box off
    set(gca,'TickDir','out')
    colorbar
end

%% look at the correlation between firing rate and speed
% normalise firing rate
UnitResponses = AM_UnitResponses_smooth(:,AM_UOI,:);
clear UnitResponses_n
for t = 1:size(UnitResponses,2)
    OneUnitResponse = squeeze(UnitResponses(:,t,:));
    UnitResponses_n(:,t,:) = OneUnitResponse./nanmax(OneUnitResponse,[],2);
end
UnitResponses = UnitResponses_n;
% normalise speed
SpeedResponse_n = AM_SpeedSel;
for t = 1:length(SpeedResponse_n)
    OneSpeedResponse_n = squeeze(SpeedResponse_n(:,t,:));
    SpeedResponse_n(:,t,:) = OneSpeedResponse_n./nanmax(OneSpeedResponse_n,[],2);
end
%%
RespUnits_idx=find(RespUnits);
for i = 1:length(RecIndex)
    Unit_RoI = find(AM_ParamSel(1,RespUnits,1)==AM_ParamSel(1,RecIndex(i),1));
    Unit_idx = RespUnits_idx(Unit_RoI);
    PertTrials_1 = AM_ParamSel(:,RecIndex(i),3)==1;
    PertTrials_0 = AM_ParamSel(:,RecIndex(i),3)==0;
    for j = 1:length(Unit_idx)
        figure
        set(gcf, 'Units', 'centimeters', 'Position', [1, 2, 30, 15], 'PaperUnits', 'centimeters', 'PaperSize', [30, 20])
        % plot unit FR responses
%         if j>1
%             subplot(length(Unit_idx),3,j+2*j-2)
%         else
%             figure
%             subplot(length(Unit_idx),3,j)
%         end
        subplot(1,3,1)
        clear UnitResponseMean
        UnitResponseMean(1,:) = squeeze(nanmean(UnitResponses_n(PertTrials_1,Unit_idx(j),:)));
        UnitResponseMean(2,:) = squeeze(nanstd(UnitResponses_n(PertTrials_1,Unit_idx(j),:)))/sqrt(sum(PertTrials_1));
        shadedErrorBar(TimeLine,UnitResponseMean(1,:),UnitResponseMean(2,:),'lineprop',{'r-','markerfacecolor',[0.5 0.5 0.5]});
        hold on
        clear UnitResponseMean
        UnitResponseMean(1,:) = squeeze(nanmean(UnitResponses_n(PertTrials_0,Unit_idx(j),:)));
        UnitResponseMean(2,:) = squeeze(nanstd(UnitResponses_n(PertTrials_0,Unit_idx(j),:)))/sqrt(sum(PertTrials_1));
        shadedErrorBar(TimeLine,UnitResponseMean(1,:),UnitResponseMean(2,:),'lineprop',{'k-','markerfacecolor',[0.5 0.5 0.5]});
        hold on
        plot([TimeLine(trialSide_samples) TimeLine(trialSide_samples)], [0 1],'k:')
        plot([TimeLine(end-trialSide_samples) TimeLine(end-trialSide_samples)], [0 1],'k:')
        plot([TimeLine(PertStartEnds(1,i)) TimeLine(PertStartEnds(1,i))], [0 1],'r:')
        plot([TimeLine(PertStartEnds(2,i)) TimeLine(PertStartEnds(2,i))], [0 1],'r:')
        xlabel('seconds'); ylabel('norm FR');
        box off
        axis tight
        if j == 1
            ll = legend({'Perturbation ON', 'Perturbation OFF'},'fontsize',10);
            ll.Color = 'none'; ll.EdgeColor = 'none';
        end
        set(gca,'TickDir','out')   
    
        subplot(1,3,2)
        clear UnitResponseMean
        UnitResponseData = UnitResponses_n(~isnan(UnitResponses_n(:,Unit_idx(j),1)),Unit_idx(j),:);
        SpeedResponseData = SpeedResponse_n(~isnan(SpeedResponse_n(:,Unit_idx(j),1)),Unit_idx(j),:);
        scatter(SpeedResponseData(:),UnitResponseData(:),'.')
        fitobject = fit(SpeedResponseData(~isnan(SpeedResponseData(:))),UnitResponseData(~isnan(SpeedResponseData(:))),'poly1');
        hold on
        plot(fitobject,SpeedResponseData(:),UnitResponseData(:))
        xlabel('cm/s'); ylabel('norm FR');
        box off
        axis tight
    
        subplot(1,3,3)
        scatter(diff(SpeedResponseData(:),5),diff(UnitResponseData(:),5),'.')
        xlabel('speed change'); ylabel('FR change');
        box off
        axis tight
        
        saveas(gcf,['X:\DATA\PROJECTS\VisPerturbation\Figures\Fig_3\panels\PertRespUnits' filesep 'Rec_' num2str(i) '_Unit_' num2str(j) '.pdf']);
        close all
    end
end

%% look at trial by trial response to see if there is a correlation whatsoever
RespUnits_idx=find(RespUnits); k = 1; clear rho rho MI_corr MI_speed MI_neuro Speed_trial
for i = 1:length(RecIndex)
    Unit_RoI = find(AM_ParamSel(1,RespUnits,1)==AM_ParamSel(1,RecIndex(i),1));
    Unit_idx = RespUnits_idx(Unit_RoI);
    PertTrials_1 = AM_ParamSel(:,RecIndex(i),3)==1;
    PertTrials_0 = AM_ParamSel(:,RecIndex(i),3)==0;
    for j = 1:length(Unit_idx)
        clear UnitResponseData SpeedResponseData
        UnitResponseData = squeeze(UnitResponses_n(PertTrials_1,Unit_idx(j),PertStartEnds(1,i)-60:PertStartEnds(2,i)));
        SpeedResponseData = squeeze(SpeedResponse_n(PertTrials_1,Unit_idx(j),PertStartEnds(1,i)-60:PertStartEnds(2,i)));
        SpeedResponseData_real = squeeze(AM_SpeedSel(PertTrials_1,Unit_idx(j),PertStartEnds(1,i)-60:PertStartEnds(2,i)));
        rho_ = diag(corr(UnitResponseData',SpeedResponseData'));
        rho{k} = rho_(~isnan(rho_));
        % MI for speed and neural response
        Speed_trial{k} = mean(SpeedResponseData_real(:,1:60),2);
        MI_neuro = (mean(UnitResponseData(:,60:end),2)-mean(UnitResponseData(:,1:60),2))./...
                   (mean(UnitResponseData(:,60:end),2)+mean(UnitResponseData(:,1:60),2));
        MI_speed = (mean(SpeedResponseData(:,60:end),2)-mean(SpeedResponseData(:,1:60),2))./...
                   (mean(SpeedResponseData(:,60:end),2)+mean(SpeedResponseData(:,1:60),2));
        MI_corr{k} = [MI_neuro, MI_speed];      
        k = k + 1;
    end
end

figure
suptitle('Correlation between speed and FR')
for k = 1:length(rho)
    subplot(12,10,k)
    histogram(abs(rho{k}),0:0.05:1)
end

figure
for k = 1:length(rho)
    CorrSpeedPert(k) = median(abs(rho{k}));
end
histogram(CorrSpeedPert,0:0.05:1)
ylabel('distribution of median R')
title('Median correlation with pert period speed, n=116')

figure
%suptitle('MIneuro VS MIspeed')
for k = 1:length(MI_corr)
    subplot(12,10,k)
    %plot(MI_corr{k}(:,2),MI_corr{k}(:,1),'.');
    hold on
    x = MI_corr{k}(:,2);
    y = MI_corr{k}(:,1);
    fitobject = fit(x(~isnan(y)),y(~isnan(y)),'poly1');
    plot(fitobject,x,y)
    mdl = fitlm(x(~isnan(y)),y(~isnan(y)));
    p_MI_slope(k,:) = [mdl.Coefficients.Estimate(2) mdl.Coefficients.pValue(2)];   
end

%% show examples
figure
x = MI_corr{k}(:,2);
y = MI_corr{k}(:,1);
fitobject = fit(x(~isnan(y)),y(~isnan(y)),'poly1');
plot(fitobject,x,y)
mdl = fitlm(x(~isnan(y)),y(~isnan(y)));
xlim([-1 1]); ylim([-1 1]);
set(gca,'TickDir','out','box','off')
xlabel('MI speed'); ylabel('MI firing rate')
title(['unit #' num2str(k)])
ll = legend({'trials' ,['y=' num2str(mdl.Coefficients.Estimate(1),1) '+' num2str(mdl.Coefficients.Estimate(2),1) 'x']},'fontsize',10);
ll.Color = 'none'; ll.EdgeColor = 'none';
text(-0.8, 0.8, ['p=' num2str(p_MI_slope(k,2)),2],'fontsize',10)
   
sum(p_MI_slope(:,2)<0.01)
find(p_MI_slope(:,2)<0.01)

figure
histogram(p_MI_slope(:,1),-1:0.1:1,'FaceColor',[0.5 0.5 0.5],'EdgeColor','none') 
hold on
histogram(p_MI_slope(p_MI_slope(:,2)<0.01,1),-1:0.1:1,'FaceColor','r','EdgeColor','none') 
ylabel('units');
xlabel('1st coeff. of fit');
set(gca,'TickDir','out','box','off')

%% check modulation wrt speed pre perturbation

figure
%suptitle('MIneuro VS speed pre-pert')
clear p_slope
for k = 1:size(Speed_trial,2)
   % subplot(12,10,k)
    % fit a linear plot
    x = Speed_trial{k};
    y = MI_corr{k}(:,1);
    %fitobject = fit(x(~isnan(y)),y(~isnan(y)),'poly1');
	%hold on
    %plot(fitobject,x,y)
    %legend off
    mdl = fitlm(x(~isnan(y)),y(~isnan(y)));
    p_slope(k,:) =  [mdl.Coefficients.Estimate(2) mdl.Coefficients.pValue(2)];
end

figure
plot(p_slope,'.')
find(p_slope<0.01)

%% show examples
figure
x = Speed_trial{k};
y = MI_corr{k}(:,1);
fitobject = fit(x(~isnan(y)),y(~isnan(y)),'poly1');
plot(fitobject,x,y)
mdl = fitlm(x(~isnan(y)),y(~isnan(y)));
%xlim([0 max(x)]); ylim([-1 1]);
set(gca,'TickDir','out','box','off')
xlabel('Speed pre-pert (cm/s)'); ylabel('MI firing rate')
title(['unit #' num2str(k)])
ll = legend({'trials' ,['y=' num2str(mdl.Coefficients.Estimate(1),1) '+' num2str(mdl.Coefficients.Estimate(2),1) 'x']},'fontsize',10);
ll.Color = 'none'; ll.EdgeColor = 'none';
text(4,0.3, ['p=' num2str(mdl.Coefficients.pValue(2),2)],'fontsize',10)
   
sum(p_slope<0.01)
find(p_slope<0.01)

figure
histogram(p_slope(:,1),-0.1:0.01:0.1,'FaceColor',[0.5 0.5 0.5],'EdgeColor','none') 
hold on
histogram(p_slope(p_slope(:,2)<0.01,1),-0.1:0.01:0.1,'FaceColor','r','EdgeColor','none') 
ylabel('units');
xlabel('1st coeff. of fit');
set(gca,'TickDir','out','box','off')




%% rate of change correlation 
RespUnits_idx=find(RespUnits); k = 1; clear rho MI_corr MI_speed MI_neuro Speed_trial
for i = 1:length(RecIndex)
    Unit_RoI = find(AM_ParamSel(1,RespUnits,1)==AM_ParamSel(1,RecIndex(i),1));
    Unit_idx = RespUnits_idx(Unit_RoI);
    PertTrials_1 = AM_ParamSel(:,RecIndex(i),3)==1;
    PertTrials_0 = AM_ParamSel(:,RecIndex(i),3)==0;
    for j = 1:length(Unit_idx)
        clear UnitResponseData SpeedResponseData
        UnitResponseData = diff(squeeze(UnitResponses_n(PertTrials_1,Unit_idx(j),PertStartEnds(1,i)-60:PertStartEnds(2,i))),5);
        SpeedResponseData = diff(squeeze(SpeedResponse_n(PertTrials_1,Unit_idx(j),PertStartEnds(1,i)-60:PertStartEnds(2,i))),5);
        SpeedResponseData_real = diff(squeeze(AM_SpeedSel(PertTrials_1,Unit_idx(j),PertStartEnds(1,i)-60:PertStartEnds(2,i))),5);
        % Pearson's correlation coefficient
        rho_ = diag(corr(UnitResponseData',SpeedResponseData'));
        rho{k} = rho_(~isnan(rho_));
        % MI for speed and neural response
        Speed_trial{k} = mean(SpeedResponseData_real(:,1:60),2);
        MI_neuro = (mean(UnitResponseData(:,60:end),2)-mean(UnitResponseData(:,1:60),2))./...
                   (mean(UnitResponseData(:,60:end),2)+mean(UnitResponseData(:,1:60),2));
        MI_speed = (mean(SpeedResponseData(:,60:end),2)-mean(SpeedResponseData(:,1:60),2))./...
                   (mean(SpeedResponseData(:,60:end),2)+mean(SpeedResponseData(:,1:60),2));
        MI_corr{k} = [MI_neuro, MI_speed];      
        k = k + 1;
    end
end

figure
suptitle('Correlation between changes of speed and FR')
for k = 1:length(rho)
    subplot(12,10,k)
    histogram(abs(rho{k}),0:0.05:1)
end


%% fit linear regression model to see if beahviour can predict trial type
% look at the speed profiles of each recording
for i = 1:length(unique(AM_ParamSel(1,:,1)))
    AllRec_idx(i) = min(find(AM_ParamSel(1,:,1)==i));
end
SpeedResponse_nSel = SpeedResponse_n(:,AllRec_idx,:);
% extract predictors from speed signals
po_i = [min(PertStartEnds(1,:)) max(PertStartEnds(2,:))];
clear DiffResponsePercent_Pert
for i = 1:size(SpeedResponse_nSel,2)
    
    UnitResponse = squeeze(SpeedResponse_nSel(:,i,:));
    
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
clear SingleUnitResponse AUC_speed
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
    AUC_speed(i,1) = AUC1;
end

tic
clear AUC_shuffled_Sp
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
        AUC_shuffled_Sp(i,sh_i) = AUC1;
    end
    i
end
toc

% AUC_shuffled_Sp = cat(2,AUC_speed(:,1),AUC_shuffled_Sp);
save('AUC_shuffled_Speed.mat','AUC_shuffled_Sp')

figure
plot(AUC_shuffled_Sp,'.')
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
AUC_shuffled_Sp = cat(2, AUC_speed, AUC_shuffled_Sp);

% prctile(AUC_shuffled(:),95) = 0.5922  on 18/09/19 p_sh(1) = 0.5922;
% prctile(AUC_shuffled(:),99) = 0.6754  on 18/09/19 p_sh(2) = 0.6754;

if ~exist('AUC_shuffled_Sp','var')
    load('AUC_shuffled_Speed.mat')
    AUC_speed = AUC_shuffled_Sp(:,1);
    AUC_shuffled_Sp = AUC_shuffled_Sp(:,2:end);
end


p_sh = prctile(AUC_shuffled_Sp(:),99);
find(AUC_speed>p_sh)

RecIndex = unique(AM_ParamSel(1,UnitsPertPos,1));
AUC_2_plot = AUC_shuffled_Sp(RecIndex,:);
AUC_lin = AUC_2_plot(:,2:end);
p_sh = prctile(AUC_lin(:),99);
%find(AUC_speed(AM_ParamSel(1,RecIndex,1))>p_sh)

%% show AUC of shuffled and original data
figure
h = histogram(AUC_2_plot(:,2:end),[0:0.02:1],'FaceColor',[0.5 0.5 0.5],'Normalization','probability','EdgeColor','none')
hold on
h1 = histogram(AUC_2_plot(:,1),[0:0.02:1],'FaceColor',[1 0 0],'Normalization','probability','EdgeColor','none')
%plot([p_sh(1) p_sh(1)], [0 max(h.BinCounts/length(AUC_shuffled_Sp(:)))*1.1],'k') 
plot([p_sh p_sh], [0 max(h1.BinCounts/length(AUC_speed(AM_ParamSel(1,RecIndex,1))))*1.1],'k--','LineWidth',3)
xlabel('AUC of shuffled data'); ylabel('probability')
ll = legend({['AUC of shuffled data n=' num2str(size(AUC_2_plot,1)) 'k'],['AUC of sessions with resp units n=' num2str(size(AUC_2_plot,1))], '99% shuffled data'})
ll.Color = 'none'; ll.EdgeColor = 'none';
box off
set(gca,'FontSize',13,'TickDir','out');

sum(AUC_speed(AM_ParamSel(1,RecIndex,1))>p_sh)

%% compute the speed score (Pearson's correlation) between units and speed
% smoothing parameters
Width = round(0.3*BonvisionFR);
Sigma = Width/3; % standard deviation in number of samples (converted from time in seconds)
x_g = linspace(-Width/2, Width/2, Width);
gaussFilter = exp(-x_g.^2/(2*Sigma^2));
gaussFilter_ = gaussFilter / sum (gaussFilter); % normalize
k = 1;
clear FR_session CorrCoeffSpeed CorrCoeffSpeed_sh
for i = 1:size(ProjectData,1)
    SessionTime = ProjectData.Session_data{i,1}.Time{1,1}; 
    MouseSpeed = conv(ProjectData.Session_data{i,1}.MouseSpeed{1,1},gaussFilter_,'same');
    for j = 1:size(ProjectData.Units_Info{i,1},1)
        FR_session{i}(j,:) = conv(histcounts(ProjectData.Units_Info{i,1}.Spiketimes{j,1},SessionTime),gaussFilter_,'same');
        Rho = corrcoef(FR_session{i}(j,:),MouseSpeed);
        CorrCoeffSpeed(k) = Rho(1,2);
%         for ll = 1:1000
%             Lag = rand*20;
%             Sh_SC = ProjectData.Units_Info{i,1}.Spiketimes{j,1}+Lag;
%             Sh_SC(Sh_SC>max(SessionTime)) = Sh_SC(Sh_SC>max(SessionTime))-max(SessionTime);
%             Sh_SC = sort(Sh_SC);
%             FR_session_sh = conv(histcounts(Sh_SC,SessionTime),gaussFilter_,'same');
%             Rho = corrcoef(FR_session_sh,MouseSpeed);
%             CorrCoeffSpeed_sh(k,ll) = Rho(1,2);
%         end
        k = k + 1;
    end
end
save('CorrCoeffSpeed_sh.mat','CorrCoeffSpeed_sh');
% shuffle the firing rate and recompute the corr coeff. 1000 times if
% select only the units worth analysing
CorrCoeffSpeed_sh = CorrCoeffSpeed_sh(AM_UOI,:);

thres_1_99 = [prctile(CorrCoeffSpeed_sh(:),1) prctile(CorrCoeffSpeed_sh(:),99)];

for u = 1:length(CorrCoeffSpeed)
    p_speed(u) = sum(CorrCoeffSpeed(u)>thres_99(u))
end
CorrCoeffSpeed = CorrCoeffSpeed(AM_UOI);
figure
h = histogram(CorrCoeffSpeed_sh,[-1:0.02:1],'FaceColor',[0.5 0.5 0.5],'Normalization','probability','EdgeColor','none')
hold on
plot([thres_1_99(1) thres_1_99(1)],[0 max(h.BinCounts)/length(CorrCoeffSpeed_sh(:))], 'k:');
plot([thres_1_99(2) thres_1_99(2)],[0 max(h.BinCounts)/length(CorrCoeffSpeed_sh(:))], 'k:');
h1 = histogram(CorrCoeffSpeed,[-1:0.02:1],'FaceColor',[1 0 0],'Normalization','probability','EdgeColor','none')
xlabel('Speed Neuro Corr.'); ylabel('probability')
box off
set(gca,'TickDir','out')
title('Speed score')
% count how many units are modulated by running
sum(CorrCoeffSpeed>thres_1_99(2))
sum(CorrCoeffSpeed<thres_1_99(1))
Units_run = [find(CorrCoeffSpeed>thres_1_99(2)) find(CorrCoeffSpeed<thres_1_99(1))];
gg = 1; clear idx
for i = 1:length(Units_run)
    if ~isempty(find(Units_run(i)==RespUnits_idx))
        idx(gg) = find(Units_run(i)==RespUnits_idx);
        gg = gg + 1;
    end
end
% recordings that contain units modulated by speed
AM_ParamSel(1,RespUnits_idx(idx),1)

AM_ParamSel(1,RecIndex,1)

%% look at the modulation indexes of speed and perturbation responses.





