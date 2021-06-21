function [predictorOut, Nqout] = prepPredictor4Regression(predictor, type, tidx, Nbins, Nbins2)
% Function to take in a predictor for regression and output an
% appropriately processed version of it:
% [predictorOut, Nqout] = prepPredictor4Regression(predictor, type, tidx, Nbins)
% predictor - an array of the variable that is to be processed
% tidx      - this is the index of the time points that are being processed
%
% It basically does:
% For 'linear' type
% - centering & normalising
% For 'linearBinning' type
% - centering & normalising
% - binning the predictor into Nbins and outputing a series of numbers
% For 'smoothedEvent' type
% - This is for discreet events like licks, rewards.
% - smoothes the events with a gaussian of width <Nbins>
% For 'basisSet' type
%   [TODO: not made yet]
%
% Julien Fournier & Aman Saleem - Sept 2019

if nargin<2
    type = 'linearBinning';
end
if nargin<3
    tidx = ones(size(predictor));
end
if isempty(tidx)
    tidx = ones(size(predictor));
end
if nargin<4
    Nbins = 5;
end
if nargin<5
    Nbins2 = 1;
end

switch type
    case 'linear'
        % This is making it median zero and normalising by 95 percentile of the
        % data
        predictor = (predictor - median(predictor(tidx)))./quantile(predictor(tidx),0.95);
       
        predictorOut = predictor;
        Nqout = 0;
    case 'linearBinning'
        predictor = (predictor - median(predictor(tidx)))./quantile(predictor(tidx),0.95);
        
        [predictorBinned, Nqout] = linearBinning_local(predictor, Nbins, tidx);
        
        predictorBins = 0:max(predictorBinned); % range of filters
        predictorOut = zeros(length(predictorBinned),length(predictorBins)-1);
        for iStep = 1:length(predictorBins)-1
            predictorOut(:,iStep) = double(predictorBinned > predictorBins(iStep) & predictorBinned <= predictorBins(iStep+1));
        end
        %         predictorOut = posFilt;
    case 'smoothedEvent'
        predictorOut = smthEvents(predictor, Nbins);
        Nqout = 0;
    case 'basisSet'
        predictor = (predictor - median(predictor(tidx)))./quantile(predictor(tidx),0.95);
        [predictorBinned, Nqout] = linearBinning_local(predictor, Nbins, tidx);
        predictorBins = 0:max(predictorBinned); % range of filters
        predictorOut = zeros(length(predictorBinned),length(predictorBins)-1);
        for iStep = 1:length(predictorBins)-1
            predictorOut(:,iStep) = double(predictorBinned > predictorBins(iStep) & predictorBinned <= predictorBins(iStep+1));
            predictorOut(:,iStep) = smthEvents(predictorOut(:,iStep), Nbins2);
        end
    otherwise
        errordlg('Wrong <type> for the predictor')
end

    function [binnedVar, Nqout] = linearBinning_local(Var,Nq,tidxref)
        varbin_start = quantile(Var(tidxref),0.05);
        varbin_end = quantile(Var(tidxref),0.95);
        varbinedges = linspace(varbin_start,varbin_end,Nq+1);
        % The binnedVar is the bin index of each point
        binnedVar = NaN(size(Var));
        for ibin = 1:length(varbinedges)-1
            binnedVar(Var>=varbinedges(ibin) & Var<varbinedges(ibin+1)) = ibin;
        end
        binnedVar(Var < varbin_start) = 1;
        binnedVar(Var >= varbin_end) = Nq;
        
        bins = sort(unique(binnedVar(~isnan(binnedVar))),'ascend');
        Nqout = numel(bins);
        % This is centering the data again
        for iq = 1:Nqout
            binnedVar(binnedVar == bins(iq)) = iq + Nq;
        end
        binnedVar = binnedVar - Nq;
    end
    function [predictorOut, Nqout] = basisSetSeparation_local(Var, Nbins, tidxref)
        
%         for ipred = 1:size([predictorOut,2)
    end
    function out = smthEvents(in, win)
        sGrid = max([100 win*5]);
        s = (-sGrid:sGrid);
        
        sfilt = (1./(win*(sqrt(2*pi)))).*exp(-(s.^2)./(2*(win^2)));
        
        % padding the inputs
        pad = zeros(size(sfilt));
        pad_length = length(pad);
        if size(in,1)~=1
            flip_again = true;
            in = in';
            
        else
            flip_again = false;
        end
        pad_in = [pad in pad];
        pad_out         = conv(pad_in, sfilt, 'same');
        out = pad_out((pad_length+1) : (end-pad_length));
       
        if flip_again
            out = out';
        end
    end
end