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

%% 2) scan through all experiments
% initialise gaussian filter parameters for firing rate
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
    
    for Stim_i = 1:max(StimToUse)
        
        [res, varList, meanSubResp_Used] = RidgeRegressionModelMain_SL(ES, Response(21,:),TrialToUse,Stim_i,plotoption); % optional: BestBeta
        
        RES_n{rec,Stim_i} = res;
        
    end
    
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
cell_oi = find(AUC_shuffled(:,1)==max(AUC_shuffled(:,1)));
SessionTrialType = Param(:,cell_oi,3); SessionTrialType = SessionTrialType(~isnan(SessionTrialType));
figure
subplot(4,1,1)
plot(linspace(-1,8.33,size(UnitResponses_smooth,3)), nanmean(squeeze(UnitResponses_smooth(SessionTrialType==1,cell_oi,:))) ,'r')
hold on
plot(linspace(-1,8.33,size(UnitResponses_smooth,3)), nanmean(squeeze(UnitResponses_smooth(SessionTrialType==0,cell_oi,:))) ,'k')
subplot(4,1,2:4)
[B,I] = sort(nanmean(squeeze(UnitResponses_smooth(1:length(SessionTrialType),cell_oi,pert_per(1):pert_per(2)))'),'descend');
imagesc(linspace(-1,8.33,size(UnitResponses_smooth,3)),1:length(SessionTrialType),squeeze(UnitResponses_smooth(I,cell_oi,:)))
hold on
imagesc(-1:0,1:length(SessionTrialType),repmat(SessionTrialType(I),1,2));
%%
ES.TF{1,1}
ES.MouseSpeed{1,1}
ES.Orientation{1,1}
ES.trialStartsEnds{1,1}
ES.Contrast{1,1}
trialID_temp = zeros(length(ES.Time{1,1}),1);
    for tt_i = 1:length(ES.trialStartsEnds{1,1})
        trialID_temp(ES.trialStartsEnds{1,1}(tt_i,1):ES.trialStartsEnds{1,1}(tt_i,2)) = tt_i;
    end

Predictor_TF = ES.TF{1,1}/max(ES.TF{1,1});
Predictor_speed = (ES.MouseSpeed{1,1}-min(ES.MouseSpeed{1,1}))/(max(ES.MouseSpeed{1,1})-min(ES.MouseSpeed{1,1}));
Predictor_Dir = ES.Orientation{1,1}/max(ES.Orientation{1,1});
Response_TBP = ((Response(21,:)-min(Response(21,:)))/(max(Response(21,:))-min(Response(21,:))))';

TrialToUse = trialID_temp(ES.PertStartsEnds{1,1}(:,1)); % only perturbation trials
TrialToUse = 1:length(ES.trialStartsEnds{1,1}); % all of the trials
tidx = ES.trialID{1}>0 & ismember(ES.trialID{1},TrialToUse) & ES.Contrast{1}>0.78;
    
% trying with binning
nr_bins = 10;
Predictor_speed = round(Predictor_speed*nr_bins)/nr_bins;
Response_TBP = round(Response_TBP*nr_bins)/nr_bins;


X = [Predictor_TF(tidx) Predictor_speed(tidx) Predictor_Dir(tidx) Response_TBP(tidx)];

    

% let's try the logistic classifier, instead of time, we use trials

Predictor_TF = ES.TF{1,1}/max(ES.TF{1,1});
Predictor_speed = (ES.MouseSpeed{1,1}-min(ES.MouseSpeed{1,1}))/(max(ES.MouseSpeed{1,1})-min(ES.MouseSpeed{1,1}));
Predictor_Dir = ES.Orientation{1,1}/max(ES.Orientation{1,1});
Response_TBP 

figure
histogram(Response_TBP(Predictor_TF==0),'Color,''r')
hold on
histogram(Response_TBP(Predictor_TF==1))
scatter(Predictor_TF(tidx),Response_TBP(tidx))



