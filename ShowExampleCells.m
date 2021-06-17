%% Tomaso Muzzu - UCL - 16/05/2019

% function to plot neural responses of units to perturbation trials

function ShowExampleCells(ProjectData,DiffDRPer_total,value2show,Units_idx,smoothing)

% Gaussian filter for smoothing out firing rate
GaussFilter_W = 0.6; % seconds
BonvisionFR = 60;
Width = round(GaussFilter_W*BonvisionFR);
Sigma = Width/3; % standard deviation in number of samples (converted from time in seconds)
x_g = linspace(-Width/2, Width/2, Width);
gaussFilter = exp(-x_g.^2/(2*Sigma^2));
gaussFilter_ = gaussFilter / sum (gaussFilter); % normalize

trialSide_seconds = 1; % take 60 samples before and after trial
trialSide_samples = 60; % take 60 samples before and after trial
cellcount = 1; 

for cc = 1:size(Units_idx,1)
    rr = Units_idx(cc,1);
    i = Units_idx(cc,2);
    trialcount = 1;
    clear SingletrialSpikes SpikeCount_temp SpikeCount_smooth edges PertRespMeanCell Pert_idx_start Pert_idx_end
    for j = 1:size(ProjectData.Trial_data{rr,1},1)
        SingletrialSpikes = ProjectData.Trial_data{rr,1}.SpikesTrial{j,1}{1,i};
        TrialStart = ProjectData.Session_data{rr,1}.ACInfo{1,1}.trialStartsEnds(j,1);
        TrialEnd = ProjectData.Session_data{rr,1}.ACInfo{1,1}.trialStartsEnds(j,2);
        edges = linspace(TrialStart-trialSide_seconds,TrialEnd+trialSide_seconds,((TrialEnd+trialSide_seconds)-(TrialStart-trialSide_seconds))*60);
        [SpikeCount_temp, Edges] = histcounts(SingletrialSpikes,edges);
        SpikeCount_smooth(j,1:length(SpikeCount_temp)) = SpikeCount_temp; %conv(SpikeCount_temp, gaussFilter_, 'same');
        if ProjectData.Trial_data{rr,1}.PerturbationON(j,1)
            Pert_idx_start(trialcount) = min(find(ProjectData.Trial_data{rr,1}.TF_in_trial{j,:}==0));
            Pert_idx_end(trialcount) = max(find(ProjectData.Trial_data{rr,1}.TF_in_trial{j,:}==0));
            PertRespMeanCell(trialcount,1) = mean(SpikeCount_smooth(j,Pert_idx_start(trialcount):Pert_idx_end(trialcount))); % mean response during pert. for ranking all response
            PertRespMeanCell(trialcount,2) = mean(SpikeCount_smooth(j,Pert_idx_end(trialcount)+1:Pert_idx_end(trialcount)+1+60)); % mean response post perturbation
            trialcount = trialcount + 1;
        end
    end
    PertRespMean(cellcount,1) = mean(PertRespMeanCell(:,1)); % mean response during pert. for ranking all response
    PertRespMean(cellcount,2) = mean(PertRespMeanCell(:,2)); % mean response post perturbation
    cellcount = cellcount + 1;
    
    clear FR
    FR(1,:) = sum(SpikeCount_smooth(ProjectData.Trial_data{rr,1}.PerturbationON==0,:)) / sum(ProjectData.Trial_data{rr,1}.PerturbationON==0);
    FR(2,:) = std(SpikeCount_smooth(ProjectData.Trial_data{rr,1}.PerturbationON==1,:)) / sqrt(length((ProjectData.Trial_data{rr,1}.PerturbationON==1)));
    FR(3,:) = sum(SpikeCount_smooth(ProjectData.Trial_data{rr,1}.PerturbationON==1,:)) / sum(ProjectData.Trial_data{rr,1}.PerturbationON==1);
    FR(4,:) = std(SpikeCount_smooth(ProjectData.Trial_data{rr,1}.PerturbationON==1,:)) / sqrt(length((ProjectData.Trial_data{rr,1}.PerturbationON==1)));
    
    figure
    TimeLine = linspace(ProjectData.Trial_data{1,1}.time_in_trial{1,1}(1)-trialSide_seconds,ProjectData.Trial_data{1,1}.time_in_trial{1,1}(end)-trialSide_seconds,size(FR,2));
    h1 = shadedErrorBar(TimeLine,smooth(FR(1,:)*60,smoothing),smooth(FR(2,:)*60,smoothing),'lineprops',{'Color',[0.5 0.5 0.5]})
    hold on
    h2 = shadedErrorBar(TimeLine,smooth(FR(3,:)*60,smoothing),smooth(FR(4,:)*60,smoothing),'lineprops',{'Color','k','markerfacecolor',[0.5 0.5 0.5]})
      
    hold on
    plot([0 0],...
        [0 max(max(FR*60))],'--r','LineWidth',0.75); % start of the trial
    hold on
    plot([ProjectData.Session_data{1,1}.ACInfo{1,1}.trialStartsEnds(1,2)-ProjectData.Session_data{1,1}.ACInfo{1,1}.trialStartsEnds(1,1)...
            ProjectData.Session_data{1,1}.ACInfo{1,1}.trialStartsEnds(1,2)-ProjectData.Session_data{1,1}.ACInfo{1,1}.trialStartsEnds(1,1)],...
            [0 max(max(FR*60))],'--r','LineWidth',0.75); % end of the trial
    hold on
    plot([TimeLine(round(mean(Pert_idx_start))) TimeLine(round(mean(Pert_idx_start)))],...
            [0 max(max(FR*60))],'--b','LineWidth',1); % start of the perturbation
    hold on
    plot([TimeLine(round(mean(Pert_idx_end))) TimeLine(round(mean(Pert_idx_end)))],...
            [0 max(max(FR*60))],'--b','LineWidth',1); % end of the perturbation
            
    xlim([-trialSide_seconds+0.3 max(TimeLine)-0.3]);
    ylim([0 max(max(FR*60))*1.2]);
    xlabel('seconds'); ylabel('FR (Hz)')
    box off
    set(gca,'TickDir','out')
    legend({['no pert (n=' num2str(sum(ProjectData.Trial_data{rr,1}.PerturbationON==0)) ')'],['with pert. (n=' num2str(sum(ProjectData.Trial_data{rr,1}.PerturbationON==1)) ')']},'Location','northwest')
    set(gca,'FontSize',15);
    DetailCell = {['Mouse: ' ProjectData.Mouse_name{rr,1}],...
                  ['recording: ' ProjectData.Recording_date{rr,1}]... 
                  ['unit ID: ' num2str(i)]};
    title(DetailCell);
end
       
    
end

