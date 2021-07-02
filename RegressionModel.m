%% Tomaso Muzzu - UCL -  20-05-2021

% Ridge regression model for spiking activity


% 1) load all data
% 2) scan through all experiments
% 3) prepare the data to be input to the model
     % use wheel speed, TF, Orientation, Direction, trial # and trial ID
     % use wheel speed as is
	 % use TF as is, or discretised between 0 and 1 (3 cycles/s)
	 % use Orientation and Direction as time-series
	 % save trial # and trial ID, I am not sure what they are useful
% 4) run the regressor and look at the performance
% 5) run the regressor multiple times removing the independent variables sequentially
% 6) save the performance and compare which indepedente variable is the best predictor


%% 1) load all data
%% Functions to plot Figure 2 - direction tuning
if ~exist('ProjectData','var')
    [ProjectData AM_UnitResponses AM_Param AM_Speed AM_UOI SelectedResponses AM_UnitResponses_smooth AM_EyeTracking] = LoadDataALL;
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


%% trial by trial

% initialise gaussian filter parameters for smoothing running speed
BonvisionFR = 60;
GaussFilter_W = 0.3; % seconds
Width = round(GaussFilter_W*BonvisionFR);
Sigma = Width/3; % standard deviation in number of samples (converted from time in seconds)
x_g = linspace(-Width/2, Width/2, Width);
gaussFilter = exp(-x_g.^2/(2*Sigma^2));
gaussFilter_ = gaussFilter / sum (gaussFilter); % normalize

Param = AM_Param(:,AM_UOI,:);
PertOnsets = Param(:,:,4); 
PertOffsets = Param(:,:,5);
pert_per(1) = round(max(PertOnsets(:))*60+60);
pert_per(2) = round(max(PertOffsets(:))*60+60);
AM_Speed_OI = AM_Speed(:,AM_UOI,:);
AM_UnitResponses_smooth_OI = AM_UnitResponses_smooth(:,AM_UOI,:);
clear R_sq Coeff_pvalue
parfor cell_oi = 1:length(DM)
    % Predictors
    Predictor_TF = squeeze(Param(~isnan(Param(:,cell_oi,1)),cell_oi,3));
    Predictor_Ori = mod(squeeze(Param(~isnan(Param(:,cell_oi,1)),cell_oi,2)),180);
    
    AM_Speed_OI_u = squeeze(AM_Speed_OI(~isnan(AM_Speed_OI(:,cell_oi,1)),cell_oi,:));
    %conv(AM_Speed_OI_u',gaussFilter_,'same');
    Predictor_speed = nanmean(AM_Speed_OI_u(:,pert_per(1):pert_per(2)),2)>2;
    
    % Predicted
    AM_UnitResponses_smooth_OI_u = squeeze(AM_UnitResponses_smooth_OI(~isnan(AM_UnitResponses_smooth_OI(:,cell_oi,1)),cell_oi,:));
    AM_UnitResponses_smooth_OI_u = AM_UnitResponses_smooth_OI_u/max(AM_UnitResponses_smooth_OI_u(:));
    Response_TBP = nansum(AM_UnitResponses_smooth_OI_u(:,pert_per(1):pert_per(2)),2);
    
    % select equal # of trials for pert-nopert conditions
    Trials_nopert = find(Predictor_TF==0);
    Trials_nopert_sel = Trials_nopert(randperm(length(Trials_nopert)));
    Trials_pert = find(Predictor_TF==1);
    trial_i = [Trials_pert; Trials_nopert_sel(1:length(Trials_pert))];
    %trial_i = trial_i(randperm(length(trial_i)));
    
    X = [Predictor_TF Predictor_Ori Predictor_speed Response_TBP];
    X = X(trial_i,:);
    
    for comb_i = 1:4
        % without interactions
        [trainedModel, validationRMSE] = trainRegressionModel(X,comb_i);
        R_sq(cell_oi,comb_i) = trainedModel.LinearModel.Rsquared.Adjusted;
        if comb_i == 4
            Coeff_pvalue(cell_oi,:) = trainedModel.LinearModel.Coefficients.pValue;  
        end
        %trainedModel.LinearModel.Rsquared.Adjusted
        % with interactions
%         [trainedModel1, validationRMSE1] = trainRegressionModel_comb(X,comb_i);
%         R_sq_c(cell_oi,comb_i) = trainedModel1.LinearModel.Rsquared.Adjusted;
%         if comb_i == 4
%             Coeff_pvalue_c(cell_oi,:) = trainedModel1.LinearModel.Coefficients.pValue;  
%         end
    end
    
    cell_oi
    
end

% plot ratios of explained variance for models using the following
% parameters:
% no combinations: no TF, no Ori, no Speed, altogether
% with combinations: no TF, no Ori, no Speed, altogether
% stat sign of parameters above: 
% no TF, no Ori, no Speed, no TF & no Ori, no TF & no Speed, no Ori & no Speed
% with combinations/interactions between the terms above
R_sq(R_sq<0) = nan;
%R_sq_c(R_sq_c<0) = nan;
R_sq_norm = R_sq./R_sq(:,4);
Coeff_pvalue_logic = Coeff_pvalue<0.05;
% R_sq_norm_c = R_sq_c./R_sq_c(:,4);
% Coeff_pvalue_c_logic = Coeff_pvalue_c<0.05;
datanames = {'no TF', 'no ori.', 'no speed', 'full model'};
hist_cols = {'r','g','b'};
%% look only at w/o interaction parameters
%% plot summary figure
figure
nr_col = 3;
for cp = 1:nr_col
    subplot(1,5,cp)
    histogram(log(R_sq_norm(:,cp)),-4:0.2:4,'FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.6)
    hold on
    histogram(log(R_sq_norm(Coeff_pvalue_logic(:,cp+1),cp)),-4:0.2:4,'FaceColor',hist_cols{cp},'FaceAlpha',0.6)
    legend({'all units n=1019', ['p<0.05 n=' num2str(sum(Coeff_pvalue_logic(:,cp+1)))]},'Location','northwest')
    legend boxoff 
    xlabel(['log (' datanames{cp} '/full)']); ylabel('Units')
    set(gca,'box','off','TickDir','out')
    title(datanames{cp})
end
% boxplots
subplot(1,5,4)
boxplot(log(R_sq_norm(:,1:end-1)),datanames(1:3),'Notch','on','orientation','vertical','outliersize',0.1)
h = findobj(gca,'Tag','Box');
bp_cols = flip(cat(2,hist_cols,'k'));
for j=1:3
    patch(get(h(j),'XData'),get(h(j),'YData'),'k','FaceAlpha',.6);
end
ylim([-3.5 2]); grid on
ylabel('log ratio of expl. var.'); title('w/o interactions (all units)')
set(gca,'box','off','TickDir','out')
% venn diagrams
subplot(1,5,5)
all_n = sum(cat(2,Coeff_pvalue_logic(:,2:end), ...
            Coeff_pvalue_logic(:,2) & Coeff_pvalue_logic(:,3), ...
            Coeff_pvalue_logic(:,2) & Coeff_pvalue_logic(:,4), ...
            Coeff_pvalue_logic(:,3) & Coeff_pvalue_logic(:,4), ...
            Coeff_pvalue_logic(:,2) & Coeff_pvalue_logic(:,3) & Coeff_pvalue_logic(:,4)));
[H,S] = venn(all_n, 'FaceColor',{'r','g','b'},'FaceAlpha',{0.6,0.6,0.6},'EdgeColor','none')
for i = 1:7
    text(S.ZoneCentroid(i,1), S.ZoneCentroid(i,2), num2str(all_n(i)))
end
legend(datanames(1:3),'Location','northwest'); legend boxoff
set(gca,'box','off','TickDir','out','visible','off')
axis equal, axis off

p_vals = [signrank(log(R_sq_norm(:,1))) signrank(log(R_sq_norm(:,2))) signrank(log(R_sq_norm(:,3)))]
p_vels = [ranksum((log(R_sq_norm(:,1))),(log(R_sq_norm(:,2)))) ...
          ranksum((log(R_sq_norm(:,1))),(log(R_sq_norm(:,3)))) ...
          ranksum((log(R_sq_norm(:,2))),(log(R_sq_norm(:,3))))];

% check how well the TF and running state can predict the variables

nanmedian(R_sq(Coeff_pvalue_logic(:,2),4)) % any TF
nanmedian(R_sq(Coeff_pvalue_logic(:,4),4)) % any running
nanmean(R_sq(Coeff_pvalue_logic(:,2) & ~Coeff_pvalue_logic(:,4),4)) % just TF & ori
nanmean(R_sq(Coeff_pvalue_logic(:,4) & ~Coeff_pvalue_logic(:,2),4)) % just running & ori
      
%% plot only the values of explained variance
figure
nr_col = 3;
for cp = 1:nr_col
    subplot(1,5,cp)
    histogram((R_sq(:,cp)),0:0.02:1,'FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.6)
    hold on
    histogram((R_sq(Coeff_pvalue_logic(:,cp+1),cp)),0:0.02:1,'FaceColor',hist_cols{cp},'FaceAlpha',0.6)
    legend({'all units n=1019', ['p<0.05 n=' num2str(sum(Coeff_pvalue_logic(:,cp+1)))]},'Location','northwest')
    legend boxoff 
    xlabel('Explained variance'); ylabel('Units')
    set(gca,'box','off','TickDir','out')
    title(datanames{cp})
end
cp = 4;
subplot(1,5,cp)
histogram((R_sq(:,cp)),0:0.02:1,'FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.6)
% hold on
% histogram((R_sq(Coeff_pvalue_logic(:,cp+1),cp)),0:0.02:1,'FaceColor',hist_cols{cp},'FaceAlpha',0.6)
legend({'all units n=1019'},'Location','northwest')
legend boxoff
xlabel('Explained variance'); ylabel('Units')
set(gca,'box','off','TickDir','out')
title(datanames{cp})
nanmean(R_sq(:,cp))
nanstd(R_sq(:,cp))/sqrt(length(R_sq(:,cp)))
nanmedian(R_sq(:,cp))



%% plot summary figure
figure
nr_col = 3;
for cp = 1:nr_col
    subplot(2,5,cp)
    histogram(R_sq_norm(:,cp),-0.5:0.05:2,'FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.6)
    hold on
    histogram(R_sq_norm(Coeff_pvalue_logic(:,cp+1),cp),-0.5:0.05:2,'FaceColor',hist_cols{cp},'FaceAlpha',0.6)
    legend({'all units n=1019', ['p<0.05 n=' num2str(sum(Coeff_pvalue_logic(:,cp+1)))]},'Location','northwest')
    legend boxoff 
    set(gca,'box','off','TickDir','out')
    title(datanames{cp})
    subplot(2,5,nr_col+2+cp)
    histogram(R_sq_norm_c(:,cp),-0.5:0.05:2,'FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.6)
    hold on
    histogram(R_sq_norm_c(Coeff_pvalue_c_logic(:,cp+1),cp),-0.5:0.05:2,'FaceColor',hist_cols{cp},'FaceAlpha',0.6)
    xlabel('ratio of expl. var.')
    legend({'all units n=1019', ['p<0.05 n=' num2str(sum(Coeff_pvalue_c_logic(:,cp+1)))]},'Location','northwest')
    legend boxoff 
    set(gca,'box','off','TickDir','out')
end
% boxplots
subplot(2,5,4)
boxplot(cat(2,R_sq_norm(:,1:end-1), R_sq(:,4)./R_sq_c(:,4)),cat(2,datanames(1:3),'without/with'),'Notch','on','orientation','vertical','outliersize',0.1)
h = findobj(gca,'Tag','Box');
bp_cols = flip(cat(2,hist_cols,'k'));
for j=1:4
    patch(get(h(j),'XData'),get(h(j),'YData'),'k','FaceAlpha',.6);
end
ylim([-1 2.5]); grid on
ylabel('ratio of expl. var.'); title('w/o interactions (all units)')
set(gca,'box','off','TickDir','out')
subplot(2,5,9)
boxplot(cat(2,R_sq_norm_c(:,1:end-1), R_sq_c(:,4)./R_sq(:,4)),cat(2,datanames(1:3),'with/without'),'Notch','on','orientation','vertical','outliersize',0.1)
h = findobj(gca,'Tag','Box');
for j=1:4
    patch(get(h(j),'XData'),get(h(j),'YData'),'k','FaceAlpha',.6);
end
ylim([-1 2.5]); grid on
ylabel('ratio of expl. var.'); title('with interactions (all units)')
set(gca,'box','off','TickDir','out')
% venn diagrams
subplot(2,5,5)
all_n = sum(cat(2,Coeff_pvalue_logic(:,2:end), ...
            Coeff_pvalue_logic(:,2) & Coeff_pvalue_logic(:,3), ...
            Coeff_pvalue_logic(:,2) & Coeff_pvalue_logic(:,4), ...
            Coeff_pvalue_logic(:,3) & Coeff_pvalue_logic(:,4), ...
            Coeff_pvalue_logic(:,2) & Coeff_pvalue_logic(:,3) & Coeff_pvalue_logic(:,4)));
[H,S] = venn(all_n, 'FaceColor',{'r','g','b'},'FaceAlpha',{0.6,0.6,0.6},'EdgeColor','none')
for i = 1:7
    text(S.ZoneCentroid(i,1), S.ZoneCentroid(i,2), num2str(all_n(i)))
end
legend(datanames(1:3),'Location','northwest'); legend boxoff
set(gca,'box','off','TickDir','out','visible','off')
axis equal, axis off
subplot(2,5,10)
all_n = sum(cat(2,Coeff_pvalue_c_logic(:,2:4), ...
            Coeff_pvalue_c_logic(:,2) & Coeff_pvalue_c_logic(:,3), ...
            Coeff_pvalue_c_logic(:,2) & Coeff_pvalue_c_logic(:,4), ...
            Coeff_pvalue_c_logic(:,3) & Coeff_pvalue_c_logic(:,4), ...
            Coeff_pvalue_c_logic(:,2) & Coeff_pvalue_c_logic(:,3) & Coeff_pvalue_c_logic(:,4)));
[H,S] = venn(all_n, 'FaceColor',{'r','g','b'},'FaceAlpha',{0.6,0.6,0.6},'EdgeColor','none');
for i = 1:7
    text(S.ZoneCentroid(i,1), S.ZoneCentroid(i,2), num2str(all_n(i)))
end
legend(datanames(1:3),'Location','northwest'); legend boxoff
set(gca,'box','off','TickDir','out','visible','off')
axis equal, axis off

%%


nanmean(log(R_sq(:,4)./R_sq_c(:,4)))


%%

R_sq_pos = R_sq(PertResp_units,:);
R_sq_pos = R_sq_pos./R_sq_pos(:,4);

R_sq_rest = R_sq(~PertResp_units,:);
R_sq_rest = R_sq_rest./R_sq_rest(:,4);

datanames = {'no TF', 'no direction', 'no speed'};
figure
for ct = 1:3
    subplot(3,1,ct)
    histogram(R_sq_pos(:,ct),-0.5:0.05:2)
    legend(datanames{ct})
end
figure
boxplot(R_sq_pos,'Notch','on')
ylim([-0.5 2])

R_sq_norm = R_sq./R_sq(:,4);
figure
for ct = 1:3
    subplot(3,1,ct)
    histogram(R_sq_norm(:,ct),-0.5:0.05:2)
    legend(datanames{ct})
end
figure
boxplot(R_sq_norm,'Notch','on')
ylim([-0.5 2])


% look at units whose TF coeff is significant
ci = 3;
%R_sq_sign = R_sq(Coeff_pvalue(:,ci+1)<0.05 & Coeff_pvalue(:,ci+2)<0.05,:)
R_sq_sign = R_sq(Coeff_pvalue(:,ci+1)<0.05,:);
R_sq_sign = R_sq_sign./R_sq_sign(:,4);

datanames = {'no TF', 'no orientation', 'no speed','full model'};
coeffnames = {'TF', 'Ori', 'Speed'};
figure
suptitle([coeffnames{ci} ' coeff p<0.05, n=' num2str(size(R_sq_sign,1))])
rows=3;
for ct = 1:rows
    subplot(rows,1,ct)
    histogram(R_sq_sign(:,ct),-0.5:0.05:2)
    legend(datanames{ct})
end

% plot explained variance
ci = 1;
R_sq_sign = R_sq(Coeff_pvalue(:,ci+1)<0.05,:);
figure
suptitle([coeffnames{ci} ' coeff p<0.05, n=' num2str(size(R_sq_sign,1))])
rows=4;
for ct = 1:rows
    subplot(rows,1,ct)
    histogram(R_sq_sign(:,ct),-0.5:0.02:1)
    legend(datanames{ct})
end



[m  cell_oi] = min(abs(R_sq(:,2)-0.7109))

AUC_shuffled(cell_oi,1)
DM(cell_oi)

%% example unit
SessionTrialType = Param(:,cell_oi,3); SessionTrialType = SessionTrialType(~isnan(SessionTrialType));
figure
subplot(4,1,1)
plot(linspace(-1,8.33,size(UnitResponses_smooth,3)), nanmean(squeeze(UnitResponses_smooth(SessionTrialType==1,cell_oi,:))) ,'r')
hold on
plot(linspace(-1,8.33,size(UnitResponses_smooth,3)), nanmean(squeeze(UnitResponses_smooth(SessionTrialType==0,cell_oi,:))) ,'k')
subplot(4,1,2:4)
[B,Ir] = sort(nanmean(squeeze(UnitResponses_smooth(1:length(SessionTrialType),cell_oi,pert_per(1):pert_per(2)))'),'descend');
imagesc(linspace(-1,8.33,size(UnitResponses_smooth,3)),1:length(SessionTrialType),squeeze(UnitResponses_smooth(Ir,cell_oi,:)))
hold on
imagesc(-1:0,1:length(SessionTrialType),repmat(SessionTrialType(Ir),1,2));






%% 2) scan through all experiments
% initialise gaussian filter parameters for smoothing running speed
BonvisionFR = 60;
GaussFilter_W = 0.3; % seconds
Width = round(GaussFilter_W*BonvisionFR);
Sigma = Width/3; % standard deviation in number of samples (converted from time in seconds)
x_g = linspace(-Width/2, Width/2, Width);
gaussFilter = exp(-x_g.^2/(2*Sigma^2));
gaussFilter_ = gaussFilter / sum (gaussFilter); % normalize

for rec = 1:size(ProjectData,1)
    %%%%%%%%%%%%%%%%%% 3) prepare data
    ES = ProjectData.Session_data{rec,1};
    % prepare FR's as array (# of neurons) X (# of independent variables)
    clear Response
    for neuron = 1:size(ES.SpikeInfo{1,1},2)-1
        Time_edges = ES.Time{1,1};
        [SpikeCount_temp, Edges] = histcounts(ES.SpikeInfo{1,1}{1,neuron+1},Time_edges);
        SpikeCount = SpikeCount_temp*BonvisionFR; % to get Hz
        UnitFR = [0 conv(SpikeCount, gaussFilter_, 'same')];
        
        Response(neuron,:) = UnitFR;
    end
    trialID_temp = zeros(length(ES.Time{1,1}),1);
    for tt_i = 1:length(ES.trialStartsEnds{1,1})
        trialID_temp(ES.trialStartsEnds{1,1}(tt_i,1):ES.trialStartsEnds{1,1}(tt_i,2)) = tt_i;
    end
    ES.trialID{1,1} = trialID_temp;
    ES.MouseSpeed{1,1} = conv(ES.MouseSpeed{1,1},gaussFilter_,'same');
    %TrialToUse = trialID_temp(ES.PertStartsEnds{1,1}(:,1)); % only perturbation trials
    TrialToUse = 1:length(ES.trialStartsEnds{1,1}); % all of the trials
    StimToUse = 4; % choose how many independent variables to use
    plotoption = 0;
    
    %%%%%%%%%%%%%%%%%% 4) run regression
    
%     for Stim_i = 1:max(StimToUse)
%         
%         [res, varList, meanSubResp_Used] = RidgeRegressionModelMain_SL(ES, Response(21,:),TrialToUse,Stim_i,plotoption); % optional: BestBeta
%         
%         RES_n{rec,Stim_i} = res;
%         
%     end
    
%     VarExplained = [[] ; VarExplained_temp];
    
end


% save('RM_res.mat','RES','-v7.3')

% load('RM_res.mat')
% contains (nr. recordings) X (nr. models)
% column 1: only temporal frequency
% column 2: only running speed
% column 3: only contrast
% column 4: only direction
% column 5: only orientation
% column 6: all variables together

rec
clear Expl_var
for i = 1:size(RES_n,2)
    Expl_var(:,i) = RES_n{rec,i}.bestPerformance;
end

GoodUnits = find(AM_UOI);
GoodUnits_rec = AM_Param(1,:,1);
GoodUnits_rec_ = find(AM_UOI(GoodUnits_rec==rec)); % only functioning units in recording indicated by rec

ExplVar_i = Expl_var(GoodUnits_rec_,:);

figure
for i = 1:size(ExplVar_i,1)
    subplot(1,size(ExplVar_i,1),i)
    plot(1:size(ExplVar_i(i,:),2),ExplVar_i(i,:))
    xlim([0 size(ExplVar_i(i,:),2)+1])
    set(gca,'xTick',1:5)
end

temp = AM_Param(1,AM_UOI,1) & PertRespUnits_pos';
GoodUnits_rec_ = find(temp(1:13))

figure
for i = 1:length(GoodUnits_rec_)
    subplot(length(GoodUnits_rec_),1,i)
    plot(1:2,res.bestBeta(1:2,GoodUnits_rec_(i)),'r')
    hold on
    plot(3,res.bestBeta(3,GoodUnits_rec(i)),'o')
    hold on
    plot(4,res.bestBeta(4,GoodUnits_rec(i)),'^')
    hold on
    plot(5:12,res.bestBeta(5:12,GoodUnits_rec(i)),'b')
    hold on
    plot(13:16,res.bestBeta(13:16,GoodUnits_rec(i)),'k')
    hold on
    xlim([0 16])
    set(gca,'xTick',[1 3 4 5 13 ],'xTicklabels',{'TF', 'run', 'contrast', 'dir', 'ori'});
    set(gca,'xTick',[1 3 4 5 13 ],'xTicklabels',{'TF', 'run', 'contrast', 'dir', 'ori'});
end 




tidx = ES.trialID{1}>0 & ismember(ES.trialID{1},TrialToUse) & ~isnan(varList(1).values) & ~isnan(varList(2).values);
for nn = 1:size(Response,1)
    CorrCoeff_temp = corrcoef(Response(nn,tidx),res.predictions.Pred(nn,:));
    CorrCoeff_FR(nn) = CorrCoeff_temp(1,2);
end
CorrCoeff_FR
input('')


figure
for gg =1:size(res.predictions.Pred,1)
    subplot(2,1,1)
    plot(meanSubResp_Used(gg,:))
    subplot(2,1,2)
    plot(res.predictions.Pred(gg,:))
    gg
    input('')
end
%
%     figure
%     for gg =1:size(res.predictions.Pred,1)
%         plot(res.bestBeta(:,gg))
%         gg
%         input('')
%     end

%% example cell
UnitResponses_smooth = AM_UnitResponses_smooth(:,AM_UOI,:);
Param = AM_Param(:,AM_UOI,:);
PertOnsets = Param(:,:,4); 
PertOffsets = Param(:,:,5);
pert_per(1) = round(max(PertOnsets(:))*60+60);
pert_per(2) = round(max(PertOffsets(:))*60+60);
SessionTrialType = Param(:,cell_oi,3); SessionTrialType = SessionTrialType(~isnan(SessionTrialType));
figure
subplot(4,1,1)
plot(linspace(-1,8.33,size(UnitResponses_smooth,3)), nanmean(squeeze(UnitResponses_smooth(SessionTrialType==1,cell_oi,:))) ,'r')
hold on
plot(linspace(-1,8.33,size(UnitResponses_smooth,3)), nanmean(squeeze(UnitResponses_smooth(SessionTrialType==0,cell_oi,:))) ,'k')
subplot(4,1,2:4)
[B,Ir] = sort(nanmean(squeeze(UnitResponses_smooth(1:length(SessionTrialType),cell_oi,pert_per(1):pert_per(2)))'),'descend');
imagesc(linspace(-1,8.33,size(UnitResponses_smooth,3)),1:length(SessionTrialType),squeeze(UnitResponses_smooth(Ir,cell_oi,:)))
hold on
imagesc(-1:0,1:length(SessionTrialType),repmat(SessionTrialType(Ir),1,2));

%% direction example

%% trial by trial

for cell_oi = 1:size()
% Predictors
Predictor_TF = squeeze(Param(~isnan(Param(:,cell_oi,1)),cell_oi,3));
Predictor_Dir = mod(squeeze(Param(~isnan(Param(:,cell_oi,1)),cell_oi,2)),180); Predictor_Dir(Predictor_Dir==0) = 180; Predictor_Dir = Predictor_Dir/max(Predictor_Dir);
AM_Speed_OI = AM_Speed(:,AM_UOI,:);
AM_Speed_OI = squeeze(AM_Speed_OI(~isnan(AM_Speed_OI(:,cell_oi,1)),cell_oi,:));
Predictor_speed = nanmean(AM_Speed_OI(:,pert_per(1):pert_per(2)),2)

% Predicted
AM_UnitResponses_smooth_OI = AM_UnitResponses_smooth(:,AM_UOI,:);
AM_UnitResponses_smooth_OI = squeeze(AM_UnitResponses_smooth_OI(~isnan(AM_UnitResponses_smooth_OI(:,cell_oi,1)),cell_oi,:));
Response_TBP = nansum(AM_UnitResponses_smooth_OI(:,pert_per(1):pert_per(2)),2);

X = [Predictor_TF Predictor_Dir Predictor_speed Response_TBP];

for comb_i = 1:4
    [trainedModel, validationRMSE] = trainRegressionModel(X,comb_i);
    R_sq(comb_i) = trainedModel.LinearModel.Rsquared.Adjusted
end


%%

% let's try the logistic classifier, instead of time, we use trials


