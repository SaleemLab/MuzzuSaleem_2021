%% Tomaso Muzzu - UCL - 16/05/2019

% function to plot neural responses of units to perturbation trials
% depending on grating direction

function ModulationIndex_AllCells(ProjectData,DiffDRPer_total,value2show,cond)

% Gaussian filter for smoothing out firing rate
GaussFilter_W = 0.3; % seconds
BonvisionFR = 60;
Width = round(GaussFilter_W*BonvisionFR);
Sigma = Width/3; % standard deviation in number of samples (converted from time in seconds)
x_g = linspace(-Width/2, Width/2, Width);
gaussFilter = exp(-x_g.^2/(2*Sigma^2));
gaussFilter_ = gaussFilter / sum (gaussFilter); % normalize

trialSide_seconds = 1; % take 60 samples before and after trial
trialSide_samples = 60; % take 60 samples before and after trial
cellcount = 1; 
clear PertRespMean_on PertRespMean_off VisRespMean
for rr = 1:size(ProjectData,1) % scan through sessions
    clear NormFR NormFR_sem NormFR_off NormFR_sem_off
    for i = 1:size(ProjectData.Units_Info{rr,1},1) % scan through units
        clear SingletrialSpikes SpikeCount_temp SpikeCount_smooth edges PertRespMeanCell
        trialcount = 1;
        for j = 1:size(ProjectData.Trial_data{rr,1},1) % scane through trials
            SingletrialSpikes = ProjectData.Trial_data{rr,1}.SpikesTrial{j,1}{1,i};
            TrialStart = ProjectData.Session_data{rr,1}.ACInfo{1,1}.trialStartsEnds(j,1);
            TrialEnd = ProjectData.Session_data{rr,1}.ACInfo{1,1}.trialStartsEnds(j,2);
            edges = linspace(TrialStart-trialSide_seconds,TrialEnd+trialSide_seconds,((TrialEnd+trialSide_seconds)-(TrialStart-trialSide_seconds))*60);
            [SpikeCount_temp, Edges] = histcounts(SingletrialSpikes,edges);
            SpikeCount_smooth(j,1:length(SpikeCount_temp)) = conv(SpikeCount_temp, gaussFilter_, 'same');
            if ProjectData.Trial_data{rr,1}.PerturbationON(j,1)
                Pert_idx_start(trialcount) = min(find(ProjectData.Trial_data{rr,1}.TF_in_trial{j,:}==0));
                Pert_idx_end(trialcount) = max(find(ProjectData.Trial_data{rr,1}.TF_in_trial{j,:}==0));
                PertRespMeanCell(trialcount,:,1) = SpikeCount_smooth(j,Pert_idx_start(trialcount):Pert_idx_start(trialcount)+60); % response during pert. for ranking all response
                PertRespMeanCell(trialcount,:,2) = SpikeCount_smooth(j,Pert_idx_end(trialcount)+1:Pert_idx_end(trialcount)+1+60); % response post perturbation
                PertRespMeanCell(trialcount,:,3) = SpikeCount_smooth(j,Pert_idx_start(trialcount)-60:Pert_idx_start(trialcount)); % response pre perturbation
                trialcount = trialcount + 1;
            end
        end 
        % perturbation trials means
        PertRespMean_on(cellcount,1) = mean(max(PertRespMeanCell(:,:,1),[],2)); % mean max response during perturbation
        PertRespMean_on(cellcount,2) = mean(max(PertRespMeanCell(:,:,2),[],2)); % mean max response post perturbation
        PertRespMean_on(cellcount,3) = mean(mean(PertRespMeanCell(:,:,3),2)); % mean response pre perturbation
        PertRespMean_on(cellcount,4) = mean(mean(PertRespMeanCell(:,:,1),2)); % mean response during perturbation
        PertRespMean_on(cellcount,5) = mean(mean(PertRespMeanCell(:,:,2),2)); % mean response post perturbation
        % non-perturbation trials means
        PertOFF_idx = find(ProjectData.Trial_data{rr,1}.PerturbationON==0);
        PertRespMean_off(cellcount,1) = mean(max(SpikeCount_smooth(PertOFF_idx,round(mean(Pert_idx_start)):round(mean(Pert_idx_end))),[],2));
        PertRespMean_off(cellcount,2) = mean(max(SpikeCount_smooth(PertOFF_idx,round(mean(Pert_idx_end))+1:round(mean(Pert_idx_end))+1+60),[],2));
        PertRespMean_off(cellcount,3) = mean(mean(SpikeCount_smooth(PertOFF_idx,round(mean(Pert_idx_start))-60:round(mean(Pert_idx_start))),2));
        PertRespMean_off(cellcount,4) = mean(mean(SpikeCount_smooth(PertOFF_idx,round(mean(Pert_idx_start)):round(mean(Pert_idx_end))),2));
        PertRespMean_off(cellcount,5) = mean(mean(SpikeCount_smooth(PertOFF_idx,round(mean(Pert_idx_end))+1:round(mean(Pert_idx_end))+1+60),2));
        
        % mean max response during perturbation period - not actual pert trial
        % mean max response post perturbation period - not actual pert trial
        % mean response pre perturbation period - not actual pert trial
        % means response during perturbation period - not actual pert trial
        % means response during posturbation period - not actual pert trial
        VisRespMean(cellcount,1) = mean(mean(SpikeCount_smooth(:,trialSide_samples+1:trialSide_samples+1+60),2));
        VisRespMean(cellcount,2) = mean(mean(SpikeCount_smooth(:,1:trialSide_samples),2));
        % mean max response to visual stimulus onset
        % mean response to 1s before visual stimulus onset
        cellcount = cellcount + 1;        
    end
end

ModI(:,1) = (PertRespMean_on(:,4)-PertRespMean_on(:,3))./(PertRespMean_on(:,4)+PertRespMean_on(:,3)); % modulation index to perturbation
ModI(:,2) = (PertRespMean_on(:,5)-PertRespMean_on(:,4))./(PertRespMean_on(:,5)+PertRespMean_on(:,4)); % modulation index to post perturbation
ModI(:,3) = (PertRespMean_off(:,4)-PertRespMean_off(:,3))./(PertRespMean_off(:,4)+PertRespMean_off(:,3)); % modulation index to perturbation period in non-pert trials
ModI(:,4) = (PertRespMean_off(:,5)-PertRespMean_off(:,4))./(PertRespMean_off(:,5)+PertRespMean_off(:,4)); % modulation index to post perturbation period in non-pert trials
ModI(:,5) = (VisRespMean(:,1)-VisRespMean(:,2))./(VisRespMean(:,1)+VisRespMean(:,2)); % modulation index to post perturbation

% compute hist
edges = 0:0.05:1;
Xaxis =0.05/2:0.05:0.05*20-0.05/2;
% Set up fittype and options.
ft = fittype( 'exp1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Normalize = 'on';
clear fitresult occ
for i = 1:size(ModI,2)
    [occ(:,i), edges] = histcounts(ModI(:,i),edges); % Pert period during pert trials
    [xData, yData] = prepareCurveData( Xaxis, occ(:,i)/sum(occ(:,i)) );
    % Fit model to data.
    [fitresult{i}, gof] = fit( xData, yData, ft, opts );
end

figure
set(gcf, 'Renderer', 'painters', 'Position',[-1500 500 1000 400]);
TitlesString = {'Perturbation period MI',...
                'Post-perturbation period MI',...
                'pert OFF Pert. MI',...
                'pert OFF post-pert. MI',...
                'Visual stimulus MI'};
subplot(1,2,1)
% Pert period during pert trials
b1 = bar(Xaxis, occ(:,1)/sum(occ(:,1)),0.8,'FaceColor','r')
b1.FaceAlpha = 0.5;
hold on
f1 = plot(fitresult{1},'r')
f1.LineWidth = 1;
hold on
% Pert period during non-pert trials
b2 = bar(Xaxis, occ(:,3)/sum(occ(:,3)),0.5,'FaceColor','b')
b2.FaceAlpha = 0.5;
hold on
f2 = plot(fitresult{3},'b')
f2.LineWidth = 1;
xlabel('Modulation Index','FontSize',12); ylabel('Fraction of units','FontSize',12);
set(gca,'TickDir','out')
box off
title(TitlesString{1},'FontSize',15); ylim([0 max(max(occ(:,1)/sum(occ(:,1))),max(occ(:,3)/sum(occ(:,3))))+0.01])
ll = legend([b1 b2],{'Pert. ON', 'Pert. OFF'},'FontSize',10)
ll.Color = 'none'; ll.EdgeColor = 'none'; 

subplot(1,2,2)
% Pert period during pert trials
b1 = bar(Xaxis, occ(:,2)/sum(occ(:,2)),0.8,'FaceColor','r')
b1.FaceAlpha = 0.5;
hold on
f1 = plot(fitresult{2},'r')
f1.LineWidth = 1;
hold on
% Pert period during non-pert trials
b2 = bar(Xaxis, occ(:,4)/sum(occ(:,4)),0.5,'FaceColor','b')
b2.FaceAlpha = 0.5;
hold on
f2 = plot(fitresult{4},'b')
f2.LineWidth = 1;
xlabel('Modulation Index','FontSize',12); ylabel('Fraction of units','FontSize',12);
set(gca,'TickDir','out')
box off
title(TitlesString{2},'FontSize',15); ylim([0 max(max(occ(:,2)/sum(occ(:,2))),max(occ(:,4)/sum(occ(:,4))))+0.01])
ll = legend([b1 b2],{'Pert. ON', 'Pert. OFF'},'FontSize',10)
ll.Color = 'none'; ll.EdgeColor = 'none'; 

figure
b1 = bar(Xaxis, occ(:,5)/sum(occ(:,5)),0.8,'FaceColor','c')
b1.FaceAlpha = 0.8;
xlabel('Modulation Index','FontSize',12); ylabel('Fraction of units','FontSize',12);
set(gca,'TickDir','out')
box off
title(TitlesString{5},'FontSize',15);


end



