%% Tomaso Muzzu - UCL - 06/02/2019

% function to plot basic PSTH to visual responses looking only at a
% specific angle for all cells

function SingleCell_response_Perturbation(Array,SpikesTrial,ES,FileName,RecordingOI,Direction, RunningThreshold, MaxSpeed,StimulusParams)

trialSide_seconds = 1; % take 60 samples before and after trial
trialSide_samples = trialSide_seconds*round(length(ES.Time)/(max(ES.Time)-min(ES.Time))); % take 60 samples before and after trial
    
for i = 1:size(SpikesTrial,1)
    clear SingletrialSpikes SpikeCount SpeedTrial
    if mod(i,30)==1
        figure('Renderer', 'painters', 'Position', [10 10 2000 1000])
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for j = 1:size(SpikesTrial,2)   % scan through the trials
        SingletrialSpikes = SpikesTrial{i,j};
        edges = linspace(ES.ACInfo.trialStartsEnds(j,1)-trialSide_seconds,ES.ACInfo.trialStartsEnds(j,2)+trialSide_seconds,...
                        ((ES.ACInfo.trialStartsEnds(j,2)+trialSide_seconds)-(ES.ACInfo.trialStartsEnds(j,1)-trialSide_seconds))*60);
        [SpikeCount(j,1:length(edges)-1), Edges] = histcounts(SingletrialSpikes,edges);
        % save speed for every trial
        if ES.trialStartsEnds(j,1)-trialSide_samples>0
            SpeedTrial_temp = ES.MouseSpeed(ES.trialStartsEnds(j,1)-trialSide_samples:ES.trialStartsEnds(j,2)+trialSide_samples);
            SpeedTrial(j,:) = SpeedTrial_temp(1:min(ES.trialStartsEnds(:,2)-ES.trialStartsEnds(:,1))+60);
        end
    end
    
    FR(1,:) = sum(SpikeCount(Array.PerturbationON==0,:)) / sum(Array.PerturbationON==0);
    FR(2,:) = sum(SpikeCount(Array.PerturbationON==1,:)) / sum(Array.PerturbationON==1);
    
    if mod(i,30)==0
        subplot(5,6,30)
    else
        subplot(5,6,mod(i,30))
    end
    TimeLine = linspace(StimulusParams(1,1)-trialSide_seconds,max(StimulusParams(:,1))-trialSide_seconds,length(FR(1,:)));
    plot(TimeLine,smooth(FR(1,:)*60,9),'Color',[0.5 0.5 0.5])
    hold on
    plot(TimeLine,smooth(FR(2,:)*60,9),'k')
      
    plot([0 0],...
        [0 max(max(FR*60))],'r','LineWidth',0.75);
    hold on
    plot([ES.ACInfo.trialStartsEnds(1,2)-ES.ACInfo.trialStartsEnds(1,1) ES.ACInfo.trialStartsEnds(1,2)-ES.ACInfo.trialStartsEnds(1,1)],...
        [0 max(max(FR*60))],'r','LineWidth',0.75);
    hold on
    plot([ES.ACInfo.PertStartsEnds(end,2)-ES.ACInfo.trialStartsEnds(max(find(Array.PerturbationON)),1) ...
          ES.ACInfo.PertStartsEnds(end,2)-ES.ACInfo.trialStartsEnds(max(find(Array.PerturbationON)),1)],...
        [0 max(max(FR*60))],'b','LineWidth',0.75);
    hold on
    plot([ES.ACInfo.PertStartsEnds(end,1)-ES.ACInfo.trialStartsEnds(max(find(Array.PerturbationON)),1) ...
          ES.ACInfo.PertStartsEnds(end,1)-ES.ACInfo.trialStartsEnds(max(find(Array.PerturbationON)),1)],...
        [0 max(max(FR*60))],'b','LineWidth',0.75);
    
    xlim([-trialSide_seconds max(StimulusParams(:,1))-trialSide_seconds]);
    upperYlim = max(max(smooth(FR(1,:)*60,9)),max(smooth(FR(2,:)*60,9)));
    ylim([0 max(1.1*upperYlim,1)]);
    box off
    if i == 1
        legend({['no pert (n=' num2str(sum(Array.PerturbationON==0)) ')'],['with pert. (n=' num2str(sum(Array.PerturbationON==1)) ')']},'Location','northwest')
    end
    if i<25
        set(gca,'XTickLabel',[])
    end
    if mod(i,30)==0 || i==size(SpikesTrial,1)
        pattern = 'Matlab';
        StartLetter = strfind(FileName,pattern);
        nrTrialsPertDir = length(find(ES.PertTrial==1));
        nrTrialsNPertDir = length(find(ES.PertTrial==0));
        DetailCell = {[FileName(StartLetter+length(pattern)+1:StartLetter+length(pattern)+13) ', ' FileName(end-22:end-4)] , ...
            char(strcat(['All Cells, Direction=' num2str(Direction) ', n_p=' num2str(nrTrialsPertDir) ', n_n_p=' num2str(nrTrialsNPertDir)]))};
        DetailCell
        a = axes;
        t1 = title(DetailCell);
        a.Visible = 'off'; % set(a,'Visible','off');
        t1.Visible = 'on'; % set(t1,'Visible','on');
    end
end  




end