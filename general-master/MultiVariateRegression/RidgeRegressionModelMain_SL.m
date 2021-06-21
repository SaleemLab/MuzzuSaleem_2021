function [res, varList,meanSubResp_Used] = RidgeRegressionModelMain_SL(ES, Response,TrialToUse,StimToUse,plotoption,BestBeta)
% ES = ES format
% Response= neuron response. same duration as ES 
% 

lambda = [0.01 0.1 1 5 10];
if nargin<6
    BestBeta=[];
end

% Setting up the shifted predictors
stepWidth = 0.3; %(in seconds) to shift continuous predictors
filtShift = 2; %1; %% in seconds, will shift forwards and/or backwards by this time window
Step = 1; %width of position filters in cm
BonvisionRate = 60; %(acquisition rate in Hz)
framesPerStep = floor(BonvisionRate * stepWidth);
ShiftVals = -floor(BonvisionRate * filtShift):framesPerStep:floor(BonvisionRate * filtShift);

%% listing all the predictors: specifying name, values, type of binning & number of bins
    i = 1; % temporal frequency
    varList(i).name      = 'TF';
    varList(i).values    = ES.TF{1,1};
    varList(i).type      = 'linearBinning'; %%%%%%%%%%%%%%%%
    varList(i).numBins   = 2;
    varList(i).numBins2  = 1;
    varList(i).shiftVals = [];

    i = 2; % running speed
    varList(i).name      = 'speed';
    varList(i).values    =  [0; ES.MouseSpeed{1,1}];
    varList(i).values(varList(i).values<0)    = 0;
    varList(i).type      = 'linear';
    varList(i).numBins   = 20;
    varList(i).numBins2  = 1;
    varList(i).shiftVals = [];

%     i = 3; % contrast
%     varList(i).name      = 'contrast';
%     varList(i).values    = ES.Contrast{1,1};
%     varList(i).type      = 'linear'; %%%%%%%%%%%%%%%%
%     varList(i).numBins   = 1;
%     varList(i).numBins2  = 5;
%     varList(i).shiftVals = [];
%     
    i = 3; % Direction in rad
    varList(i).name      = 'dir';
    varList(i).values    = ES.Orientation{1,1};
    varList(i).type      = 'linearBinning'; %%%%%%%%%%%%%%%%
    varList(i).numBins   = 8;
    varList(i).numBins2  = 5;
    varList(i).shiftVals = [];
%     
%     i = 5; % Orientation in rad
%     varList(i).name      = 'ori';
%     varList(i).values    = mod(ES.Orientation{1,1},pi);
%     varList(i).type      = 'linearBinning'; %%%%%%%%%%%%%%%%
%     varList(i).numBins   = 4;
%     varList(i).numBins2  = 5;
%     varList(i).shiftVals = [];

%% Processing the predicted variable (response)
meanSubResp = (Response - repmat(nanmean(Response,2),[1 size(Response,2)])); % used for model fitting

%% Processing each of the predictors
% tidx = ES.traj>0 & ES.contrast>0 & ~isnan(varList(1).values) & ~isnan(varList(2).values) & (varList(2).values>1) & ~isnan(meanSubResp');
tidx = ES.trialID{1}>0 & ismember(ES.trialID{1},TrialToUse) ...
       & ~isnan(varList(1).values)  ...
        & ~isnan(varList(2).values) ...
        & ~isnan(varList(2).values) ...
        & ES.Contrast{1}>0.78;
%tidx = true(length(ES.trialID{1}),1);

[stimMat, stimID] = prepRegStimMatrix(varList, tidx);

%% Run the regression
% PosPredIdx = 1:varList(1).numBins;
if StimToUse>3
    activeStim = stimID>=1; % This is basically considering all the stimuli(stimID==1) | (stimID==2);% This is position and speed; stimID~=3
else
    activeStim = stimID~=StimToUse; % it is actually stimulus to avoid
end
    
meanSubResp_Used=meanSubResp(:,tidx);
[res.bestBeta, res.bestPerformance, res.predictions, res.bestL] = runRidgeRegression_SL([stimMat(tidx,activeStim) stimMat(tidx,end)], meanSubResp(:,tidx), lambda, ES.trialID{1}(tidx),BestBeta);

%% Plot the betas (W hat)
if plotoption==1
figure
i = 1;
for iVar = 1:length(varList)
    idxs = find(stimID==iVar);
    activeStimID = stimID(activeStim);
    if sum(activeStimID==iVar)>0
        plot(idxs, res.bestBeta(activeStimID==iVar),'-o');
        hold on;
        legend_entry{i} = varList(iVar).name;
        i = i + 1;
    end
end
legend(legend_entry);
else
end
        
end