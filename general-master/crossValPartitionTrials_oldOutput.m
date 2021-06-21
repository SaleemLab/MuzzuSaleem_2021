function CVO = crossValPartitionTrials_oldOutput(trialID, kfold, sepCV, alt)
% CVO = crossValPartitionTrials_oldOutput(trialID, kfold, sepCV, alt)
% Running the cross validation on trials basis, but getting the output in
% the old format, which is in time.

allTrials = unique(trialID(~isnan(trialID)));
if nargin==2
    CVO_Trials = crossValPartitionTrials(allTrials, kfold);
elseif nargin==3
    CVO_Trials = crossValPartitionTrials(allTrials, kfold, sepCV);
elseif nargin==4
    CVO_Trials = crossValPartitionTrials(allTrials, kfold, sepCV, alt);
end

% To convert to the old format
for ik = 1:kfold
    CVO.train{ik}= false(size(trialID));
    CVO.test{ik} = false(size(trialID));
    CVO.cv{ik}   = false(size(trialID));
    for iTrial = 1:length(CVO_Trials.train{ik})
        CVO.train{ik}(trialID==CVO_Trials.train{ik}(iTrial)) = true;
    end
    for iTrial = 1:length(CVO_Trials.test{ik})
        CVO.test{ik}(trialID==CVO_Trials.test{ik}(iTrial)) = true;
    end
    for iTrial = 1:length(CVO_Trials.cv{ik})
        CVO.cv{ik}(trialID==CVO_Trials.cv{ik}(iTrial)) = true;
    end
end
CVO.kfold = kfold;
CVO.separateCV = true;
end