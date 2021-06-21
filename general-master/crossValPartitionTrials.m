function CVO = crossValPartitionTrials(trialIDs, kfold, sepCV, alt)
% CVO = crossValPartitionTrials(trialIDs, kfold, sepCV, alt)
% trialIDs need to be unique numbers, not necessarily continous.
% default kfold is 10
% 'alt'   - is to have alternate trials in train & test/cv sets
% 'sepCV' - is whether there are separate test and cv sets
%
% Example Usage:
%         allTrials = unique(trialID(~isnan(trialID)));
%         CVO_Trials = crossValPartitionTrials(allTrials, kfold);
%         % To convert to the old format
%         for ik = 1:kfold
%             CVO.train{ik}= false(size(trialID));
%             CVO.test{ik} = false(size(trialID));
%             CVO.cv{ik}   = false(size(trialID));
%             for iTrial = 1:length(CVO_Trails.train{ik})
%                 CVO.train{ik}(trialID==CVO_Trails.train{ik}(iTrial)) = true;
%             end
%             for iTrial = 1:length(CVO_Trails.test{ik})
%                 CVO.test{ik}(trialID==CVO_Trails.test{ik}(iTrial)) = true;
%             end
%             for iTrial = 1:length(CVO_Trails.cv{ik})
%                 CVO.cv{ik}(trialID==CVO_Trails.cv{ik}(iTrial)) = true;
%             end
%         end
%         CVO.kfold = kfold;
%         CVO.separateCV = true;
% 
% Output:
% CVO - output stucture
% CVO.kfold = kfold;
% CVO.separateCV = true/false;
% CVO.train{n};
% CVO.test{n};
% CVO.cv{n};
%
% Aman Saleem - Sept 2019

if nargin<2
    kfold = 10;
end
if nargin<3
    sepCV = 1;
end
if nargin<4
    alt = 0;
end

if ~alt
    if length(trialIDs)<=(2*kfold)
        disp('WARNING: less than 2 trials in the test set')
        if length(trialIDs)<(kfold)
            errordlg('Less trials than kfold')
        end
    end
    
    rp = randperm(length(trialIDs));
    % spliting into k quantile bins
    binEnds = round(quantile(1:length(trialIDs), [(1/kfold):(1/kfold):1]));
    
    val = ones(1,length(trialIDs));
    for i = 2:length(binEnds)
        val(binEnds(i-1):binEnds(i)) = i;
    end
    
    shf_trialIDs(rp) = val;
    
    CVO.kfold = kfold;
    CVO.separateCV = sepCV;
    
    for ik = 1:kfold
        CVO.test{ik} = (shf_trialIDs==ik);
        if sepCV
            if ik<kfold
                CVO.cv{ik} = (shf_trialIDs==ik+1);
            else
                CVO.cv{ik} = (shf_trialIDs==1);
            end
        else
            CVO.cv{ik} = CVO.test{ik};
        end
        CVO.train{ik} = ~CVO.test{ik} & ~CVO.cv{ik};
        
        CVO.test{ik}    = trialIDs(CVO.test{ik});
        CVO.train{ik}   = trialIDs(CVO.train{ik});
        CVO.cv{ik}      = trialIDs(CVO.cv{ik});
    end
else
%     This is just to have alternate trials in the two sets
    CVO.kfold = 2;
    CVO.separateCV = false;
    CVO.test{1} = trialIDs(1:2:end);
    CVO.train{1} = trialIDs(2:2:end);
    CVO.cv{1} = trialIDs(1:2:end);
    
    CVO.test{2} = trialIDs(2:2:end);
    CVO.train{2} = trialIDs(1:2:end);
    CVO.cv{2} = trialIDs(2:2:end);
end