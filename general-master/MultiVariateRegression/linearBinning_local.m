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