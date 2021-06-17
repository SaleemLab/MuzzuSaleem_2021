%% Tomaso Muzzu - UCL - 16/05/2019

% function to plot neural responses of units to perturbation trials

function ShowExampleCells_RunStill(ProjectData,DiffDRPer_total,value2show,Units_idx,smoothing)

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
    %% compute mean speed
    c1 = 1; clear PertStart_i SpikeCount_FR SpikeCount_FR_init SpikeCount_FR_pert SpikeCount_FR_postpert
    for ss = 1: size(SpikeCount_smooth,1)
        if ProjectData.Trial_data{rr,1}.PerturbationON(ss,1)==1
            PertStart_i(c1,1) = min(find(ProjectData.Trial_data{rr,1}.TF_in_trial{ss}==0)); % find index start of perturbation
            PertStart_i(c1,2) = max(find(ProjectData.Trial_data{rr,1}.TF_in_trial{ss}==0)); % find  duration of perturbation
            PertStart_i(c1,3) = ss; % perturbation trial ID
            PertStart_i(c1,4) = mean(ProjectData.Trial_data{rr,1}.Speed_in_trial{PertStart_i(c1,3)}(PertStart_i(c1,1):PertStart_i(c1,2))); % mean speed
            c1 = c1 + 1;
        end
    end
    % save mean speeds over the three periods of interest: initial 4
    % seconds, perturbation time window, 1s post perturbation
    OrigPertTrialOrder = ProjectData.Trial_data{rr,1}.PerturbationON;
    clear TrialSpeeds
    for ss = 1:size(SpikeCount_smooth,1)
        TrialSpeeds(ss,:) = [ mean(ProjectData.Trial_data{rr,1}.Speed_in_trial{ss,1}(round(trialSide_seconds*BonvisionFR):round(trialSide_seconds*BonvisionFR)+60*4)),... % initial 4 s
            mean(ProjectData.Trial_data{rr,1}.Speed_in_trial{ss,1}(round(mean(PertStart_i(:,1))):round(mean(PertStart_i(:,2))))),...
            mean(ProjectData.Trial_data{rr,1}.Speed_in_trial{ss,1}(round(mean(PertStart_i(:,2))):round(mean(PertStart_i(:,2))+60)))];
    end
    
    
    %%
    clear FR
    FR(1,:) = sum(SpikeCount_smooth(ProjectData.Trial_data{rr,1}.PerturbationON==0 & TrialSpeeds(:,2)<3,:))... % non pert, still MEAN
        / sum(ProjectData.Trial_data{rr,1}.PerturbationON==0 & TrialSpeeds(:,2)<3);
    FR(2,:) = std(SpikeCount_smooth(ProjectData.Trial_data{rr,1}.PerturbationON==0 & TrialSpeeds(:,2)<3,:))... % non pert, still STD
        / sqrt(length((ProjectData.Trial_data{rr,1}.PerturbationON==0 & TrialSpeeds(:,2)<3)));
    FR(3,:) = sum(SpikeCount_smooth(ProjectData.Trial_data{rr,1}.PerturbationON==1 & TrialSpeeds(:,2)<3,:))... % pert, still MEAN
        / sum(ProjectData.Trial_data{rr,1}.PerturbationON==1 & TrialSpeeds(:,2)<3);
    FR(4,:) = std(SpikeCount_smooth(ProjectData.Trial_data{rr,1}.PerturbationON==1 & TrialSpeeds(:,2)<3,:))... % pert, still STD
        / sqrt(length((ProjectData.Trial_data{rr,1}.PerturbationON==1 & TrialSpeeds(:,2)<3)));
    FR(5,:) = sum(SpikeCount_smooth(ProjectData.Trial_data{rr,1}.PerturbationON==0 & TrialSpeeds(:,2)>3,:))... % non pert, run MEAN
        / sum(ProjectData.Trial_data{rr,1}.PerturbationON==0 & TrialSpeeds(:,2)<3);
    FR(6,:) = std(SpikeCount_smooth(ProjectData.Trial_data{rr,1}.PerturbationON==0 & TrialSpeeds(:,2)>3,:))... % non pert, run STD
        / sqrt(length((ProjectData.Trial_data{rr,1}.PerturbationON==0 & TrialSpeeds(:,2)<3)));
    FR(7,:) = sum(SpikeCount_smooth(ProjectData.Trial_data{rr,1}.PerturbationON==1 & TrialSpeeds(:,2)>3,:))... % pert, run MEAN
        / sum(ProjectData.Trial_data{rr,1}.PerturbationON==1 & TrialSpeeds(:,2)<3);
    FR(8,:) = std(SpikeCount_smooth(ProjectData.Trial_data{rr,1}.PerturbationON==1 & TrialSpeeds(:,2)>3,:))... % pert, run STD
        / sqrt(length((ProjectData.Trial_data{rr,1}.PerturbationON==1 & TrialSpeeds(:,2)<3)));
    
    figure
    TimeLine = linspace(ProjectData.Trial_data{1,1}.time_in_trial{1,1}(1)-trialSide_seconds,ProjectData.Trial_data{1,1}.time_in_trial{1,1}(end)-trialSide_seconds,size(FR,2));
    h1 = shadedErrorBar(TimeLine,smooth(FR(7,:)*60,smoothing),smooth(FR(7,:)*60,smoothing),'lineprops',{'Color',[0.5 0.5 0.5]});
    hold on
    h2 = shadedErrorBar(TimeLine,smooth(FR(3,:)*60,smoothing),smooth(FR(4,:)*60,smoothing),'lineprops',{'Color','k','markerfacecolor',[0.5 0.5 0.5]});
      
    hold on
    plot([0 0],...
        [0 max(max(smooth(FR(3,:)*60,smoothing))*1.5,5)],'--r','LineWidth',0.75); % start of the trial
    hold on
    plot([ProjectData.Session_data{1,1}.ACInfo{1,1}.trialStartsEnds(1,2)-ProjectData.Session_data{1,1}.ACInfo{1,1}.trialStartsEnds(1,1)...
            ProjectData.Session_data{1,1}.ACInfo{1,1}.trialStartsEnds(1,2)-ProjectData.Session_data{1,1}.ACInfo{1,1}.trialStartsEnds(1,1)],...
            [0 max(max(smooth(FR(3,:)*60,smoothing))*1.5,5)],'--r','LineWidth',0.75); % end of the trial
    hold on
    plot([TimeLine(round(mean(Pert_idx_start))) TimeLine(round(mean(Pert_idx_start)))],...
            [0 max(max(smooth(FR(3,:)*60,smoothing))*1.5,5)],'--b','LineWidth',1); % start of the perturbation
    hold on
    plot([TimeLine(round(mean(Pert_idx_end))) TimeLine(round(mean(Pert_idx_end)))],...
            [0 max(max(smooth(FR(3,:)*60,smoothing))*1.5,5)],'--b','LineWidth',1); % end of the perturbation
            
    xlim([-trialSide_seconds+0.3 max(TimeLine)-0.3]);
    ylim([0 max(max(smooth(FR(3,:)*60,smoothing))*1.5,5)]);
    xlabel('seconds'); ylabel('FR (Hz)')
    box off
    set(gca,'TickDir','out')
    clear  n_PT_cond
    n_PT_cond(1) = sum(ProjectData.Trial_data{rr,1}.PerturbationON==1 & TrialSpeeds(:,2)>3);
    n_PT_cond(2) = sum(ProjectData.Trial_data{rr,1}.PerturbationON==1 & TrialSpeeds(:,2)<3);
    legend({['Running (n=' num2str(n_PT_cond(1)) ')'],['still (n=' num2str(n_PT_cond(2)) ')']},'Location','northwest')
    set(gca,'FontSize',15);
    DetailCell = {['Mouse: ' ProjectData.Mouse_name{rr,1}],...
                  ['recording: ' ProjectData.Recording_date{rr,1}]... 
                  ['unit ID: ' num2str(i)]};
    title(DetailCell);
end
       
    
end

