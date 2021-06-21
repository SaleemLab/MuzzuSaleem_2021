function [stimMat, stimID] = prepRegStimMatrix(varList, tidx)
% stimMat has a constant at the end.
% Format is varList example is: (Nav==es)
% i = 1; % position
% varList(i).name      = 'position';
% varList(i).values    = Nav.traj;
% varList(i).type      = 'linearBinning';
% varList(i).numBins   = 50;
% varList(i).shiftVals = [];
%
% i = 2; % speed
% varList(i).name      = 'speed';
% varList(i).values    = Nav.ballspeed;
% varList(i).values(varList(i).values<0)    = 0;
% varList(i).type      = 'linear';
% varList(i).numBins   = 5;
% varList(i).shiftVals = [];
%
% i = 3; % lick
% varList(i).name      = 'licks';
% varList(i).values    = Nav.lick;
% varList(i).type      = 'smoothedEvent';
% varList(i).numBins   = 1;
% varList(i).shiftVals = [];
%
% i = 4; % reward
% varList(i).name      = 'reward';
% varList(i).values    = double(~isnan(Nav.reward));
% varList(i).type      = 'smoothedEvent';
% varList(i).numBins   = 1;
% varList(i).shiftVals = ShiftVals;

for iVar = 1:length(varList)
    disp(['Processing the predictor: ' varList(iVar).name]);
    [predictors{iVar}, Nqout(iVar)] = prepPredictor4Regression(...
        varList(iVar).values, varList(iVar).type, tidx, varList(iVar).numBins, varList(iVar).numBins2);
    if ~isempty(varList(iVar).shiftVals)
        tempPredictor = predictors{iVar};
        predictors{iVar} = [];
        for iShift = 1:length(varList(iVar).shiftVals)
            predictors{iVar} = [predictors{iVar} circshift(tempPredictor,varList(iVar).shiftVals(iShift),1)];
        end
    end
end

% Making the stimulus matrix (inputs)
stimMat = [];
stimID = [];
for iVar = 1:length(predictors)
    stimMat = [stimMat predictors{iVar}];
    stimID  = [stimID  iVar*ones(1,size(predictors{iVar},2))];
end

stimMat = [stimMat ones(size(predictors{1},1),1)]; % adding a constant for the fitting

end




