function [bestBeta, bestPerformance, predictions, bestL] = runRidgeRegression_SL(stimMat, Response, lambda_in, trialID,BestBeta)

%%%I have addedded the possibility of setting the beta values in the input
%%%and select the option of using them instead of hestimating them. To do
%%%so we go through the same initialization of the variables but then in
%%%the kfold loop we just calculate the Ylambda and the explained variance

if nargin<4
    trialID = [];
end

if ~isempty(BestBeta)
    useBeta=1;
else
    useBeta=0;
end

nCells  = size(Response,1);
nlam    = size(lambda_in,2);
nPreds  = size(stimMat,2);
nl_iter = 5; % number of iterations to refine lamda
% AllPredidx = 1:nPreds;

% best lambda
% if sum(stimMat(:))>0
% getting the cross-validataion set
kfold = 10; 
sepCV = 1;
if isempty(trialID)
    CVO = crossValPartition(ones(1,size(stimMat,1)), kfold,sepCV);
else
    if length(trialID)~=size(stimMat,1)
        errordlg('trialID length does not match stimMat');
        return
    end
    CVO = crossValPartitionTrials_oldOutput(trialID, kfold,sepCV);
end

% Initializing empty variables
Ypred_lambda = cell(nlam,1);
Ypred_lambda_nonPos = cell(nlam,1);
Ypred_lambda_Pos = cell(nlam,1);
betaHat = cell(nlam,1);
for il = 1:nlam
    Ypred_lambda{il} = NaN(size(Response));
    Ypred_lambda_nonPos{il} = NaN(size(Response));
    Ypred_lambda_Pos{il} = NaN(size(Response));
    betaHat{il} = NaN(kfold,nPreds,nCells);
end
varExpl = zeros(nCells, nlam);

%     Getting the training data
for iter = 1:kfold
    Xtest{iter} = stimMat(~CVO.train{iter} ,:);
    Ytest = Response(:,~CVO.train{iter})';
    Xtrain = stimMat(CVO.train{iter} ,:);
    Ytrain = Response(:,CVO.train{iter})';
    
    HessMat{iter} = Xtrain' * Xtrain;
    CorrMat{iter} = (Xtrain' * Ytrain);
end

lambda = repmat(lambda_in,[nCells, 1]);
%%%
%%% calculate EV for each lambda and each cell
%%%
%%% iterate over ranges of lamda
if useBeta==0
    for klambda = 1:nl_iter % not realted to nlam but just number of iterations to best lamda.
        if klambda == 1
            cellidx = {1:nCells};
            lambda_cellgroup = lambda(1,:);
        end
        
        for kcell = 1:numel(cellidx)
            %  train each cell across all folds
            for iter = 1:kfold
                for il = 1:nlam
                    beta = ...
                        (HessMat{iter} + eye(size(Xtrain,2)) * lambda_cellgroup(kcell,il)*(trace(HessMat{iter})/nPreds)) \ CorrMat{iter}(:,cellidx{kcell});
                    
                    
                    Ypred_lambda{il}(cellidx{kcell},~CVO.train{iter}) = (Xtest{iter}*beta)';
                    %                 if klambda == nlam
                    %                     if ~isempty(nonPosPredidx)
                    %                         Ypred_lambda_nonPos{il}(cellidx{kcell},~CVO.train{iter}) = (Xtest{iter}(:,nonPosPredidx)*beta(nonPosPredidx,:))';%calciumTraces(cellidx{kcell},~CVO.train{iter}) - (Xtest{iter}(:,PosPredidx)*beta(PosPredidx,:))';
                    %                     end
                    %                     if ~isempty(PosPredidx)
                    %                         Ypred_lambda_Pos{il}(cellidx{kcell},~CVO.train{iter}) = (Xtest{iter}(:,PosPredidx)*beta(PosPredidx,:))';%calciumTraces(cellidx{kcell},~CVO.train{iter}) - (Xtest{iter}(:,nonPosPredidx)*beta(nonPosPredidx,:))';
                    %                     end
                    %                 end
                    
                    betaHat{il}(iter,:,cellidx{kcell}) = beta;
                end
            end
            
            % test and calculate EV for each cell
            m_calciumTraces = repmat(mean(Response(cellidx{kcell},:), 2),[1 size(Response,2)]);
            for il = 1:nlam
                m_Ypred = repmat(mean(Ypred_lambda{il}(cellidx{kcell},:), 2),[1 size(Ypred_lambda{il},2)]);
                %             m_calciumTraces = repmat(mean(Response(cellidx{kcell},:), 2),[1 size(Response,2)]);
                varExpl(cellidx{kcell},il) = 1 - mean(((Ypred_lambda{il}(cellidx{kcell},:) - m_Ypred) - (Response(cellidx{kcell},:) - m_calciumTraces)).^2, 2)./mean((Response(cellidx{kcell},:) - m_calciumTraces).^2, 2);
            end
        end
        
        %%%% re-define the lamda values
        
        if klambda < nlam
            lambda_cellgroup = [];
            cellidx = [];
            for iCell  = 1:nCells
                [~, bestL] = max(varExpl(iCell,:));
                if bestL == 1
                    lambda(iCell,:) = [lambda(iCell,1)/100 lambda(iCell,1)/10 lambda(iCell,1) lambda(iCell,1) + (lambda(iCell,2)-lambda(iCell,1))/2 lambda(iCell,2)];
                elseif bestL == nlam
                    lambda(iCell,:) = [lambda(iCell,nlam-1) lambda(iCell,nlam) lambda(iCell,nlam)*2 lambda(iCell,nlam)*3 lambda(iCell,nlam)*4];
                else
                    lambda(iCell,:) = [lambda(iCell,bestL-1) lambda(iCell,bestL-1)+(lambda(iCell,bestL)-lambda(iCell,bestL-1))/2 lambda(iCell,bestL) lambda(iCell,bestL)+(lambda(iCell,bestL+1)-lambda(iCell,bestL))/2 lambda(iCell,bestL+1)];
                end
                Fnew = true;
                for k = 1:size(lambda_cellgroup,1)
                    if isequal(lambda(iCell,:),lambda_cellgroup(k,:))
                        Fnew = false;
                        cellidx{k} = [cellidx{k} iCell];
                    end
                end
                if Fnew
                    lambda_cellgroup = [lambda_cellgroup;lambda(iCell,:)];
                    cellidx{size(lambda_cellgroup,1)} = iCell;
                end
            end
        end
        %%% this part
    end
    
else
    cellidx = {1:nCells};
    for kcell = 1:numel(cellidx)
        for iter = 1:kfold
            Ypred_lambda{1}(cellidx{kcell},~CVO.train{iter}) = (Xtest{iter}*BestBeta)';
        end
        m_calciumTraces = repmat(mean(Response(cellidx{kcell},:), 2),[1 size(Response,2)]);
        m_Ypred = repmat(mean(Ypred_lambda{1}(cellidx{kcell},:), 2),[1 size(Ypred_lambda{1},2)]);
        varExpl(cellidx{kcell},1) = 1 - mean(((Ypred_lambda{1}(cellidx{kcell},:) - m_Ypred) - (Response(cellidx{kcell},:) - m_calciumTraces)).^2, 2)./mean((Response(cellidx{kcell},:) - m_calciumTraces).^2, 2);
    end
end
% end

bestBeta = zeros(nPreds, nCells);
bestL = zeros(nCells,1);
bestPerformance = zeros(nCells,1);
predictions.Pred = zeros(size(Response));
% predictions.PredPos = zeros(size(Response));
% predictions.PrednonPos = zeros(size(Response));
% if sum(stimMat(:))>0
for iCell  = 1:nCells
    if useBeta==0
    [bestPerformance(iCell), ibestL] = max(varExpl(iCell,:), [],2);
    bestL(iCell) = lambda(iCell,ibestL);
    bestBeta(:,iCell) = mean(betaHat{ibestL}(:,:,iCell), 1);
    predictions.Pred(iCell,:) = Ypred_lambda{ibestL}(iCell,:);
%     predictions.PredPos(iCell,:) = Ypred_lambda_Pos{ibestL}(iCell,:);
%     predictions.PrednonPos(iCell,:) = Ypred_lambda_nonPos{ibestL}(iCell,:);
    else
     bestPerformance(iCell)= varExpl(iCell,1);
     bestL(iCell)=NaN;
     bestBeta(:,iCell)=BestBeta(:,iCell);
     predictions.Pred(iCell,:) = Ypred_lambda{1}(iCell,:);
    end
end
% end