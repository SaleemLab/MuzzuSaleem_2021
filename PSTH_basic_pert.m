%% Tomaso Muzzu - UCL - 09/01/2018

% function to plot basic PSTH to visual responses 


function PSTH_basic_pert(i,SpikesTrial,ES,FileName,RecordingOI)

    clear SingletrialSpikes SpikeCount SpeedTrial
    
    figure('Renderer', 'painters', 'Position', [10 10 800 1000])
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(6,1,4:6) % plot single spikes
    hold off
    for j = 1:size(SpikesTrial,2)   % scan through the trials
        SingletrialSpikes = SpikesTrial{i,j};
        if ES.PertTrial(j) == 0
            plot(SingletrialSpikes-ES.ACInfo.trialStartsEnds(j,1),j*ones(length(SingletrialSpikes),1),'.','Color',[0.6 0.6 0.6]);
        else
            plot(SingletrialSpikes-ES.ACInfo.trialStartsEnds(j,1),j*ones(length(SingletrialSpikes),1),'.','Color',[0 0 0]);
            lastPertTrial = j;
        end
        hold on
        edges = linspace(ES.ACInfo.trialStartsEnds(j,1)-0.5,ES.ACInfo.trialStartsEnds(j,2)+0.5,((ES.ACInfo.trialStartsEnds(j,2)+0.5)-(ES.ACInfo.trialStartsEnds(j,1)-0.5))*60);
        [SpikeCount(j,1:length(edges)-1), Edges] = histcounts(SingletrialSpikes,edges); 
        % save speed for every trial
        SpeedTrial_temp = ES.MouseSpeed(ES.trialStartsEnds(j,1)-30:ES.trialStartsEnds(j,2)+30);
        SpeedTrial(j,:) = SpeedTrial_temp(1:min(ES.trialStartsEnds(:,2)-ES.trialStartsEnds(:,1))+60);
    end
    plot([ES.ACInfo.trialStartsEnds(j,2)-ES.ACInfo.trialStartsEnds(j,1) ES.ACInfo.trialStartsEnds(j,2)-ES.ACInfo.trialStartsEnds(j,1)],...
        [0 size(SpikesTrial,2)],'r','LineWidth',2.5);
    hold on
    plot([ES.ACInfo.trialStartsEnds(j,1)-ES.ACInfo.trialStartsEnds(j,1) ES.ACInfo.trialStartsEnds(j,1)-ES.ACInfo.trialStartsEnds(j,1)],...
        [0 size(SpikesTrial,2)],'r','LineWidth',2.5);
    hold on
    plot([ES.ACInfo.PertStartsEnds(end,2)-ES.ACInfo.trialStartsEnds(lastPertTrial,1) ES.ACInfo.PertStartsEnds(end,2)-ES.ACInfo.trialStartsEnds(lastPertTrial,1)],...
        [0 size(SpikesTrial,2)],'k','LineWidth',2);
    hold on
    plot([ES.ACInfo.PertStartsEnds(end,1)-ES.ACInfo.trialStartsEnds(lastPertTrial,1) ES.ACInfo.PertStartsEnds(end,1)-ES.ACInfo.trialStartsEnds(lastPertTrial,1)],...
        [0 size(SpikesTrial,2)],'k','LineWidth',2);
    xlim([-0.5 ES.ACInfo.trialStartsEnds(j,2)-ES.ACInfo.trialStartsEnds(j,1)+0.5])
    xlabel('seconds'); ylabel('trial'); set(gca,'TickDir','out')
    box off
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(6,1,1) % plot contrast and TF
    hold off
    plot([ES.ACInfo.trialStartsEnds(j,2)-ES.ACInfo.trialStartsEnds(j,1) ES.ACInfo.trialStartsEnds(j,2)-ES.ACInfo.trialStartsEnds(j,1)],...
        [0 3.2],'r','LineWidth',2.5,'HandleVisibility','off');
    hold on
    plot([ES.ACInfo.trialStartsEnds(j,1)-ES.ACInfo.trialStartsEnds(j,1) ES.ACInfo.trialStartsEnds(j,1)-ES.ACInfo.trialStartsEnds(j,1)],...
        [0 3.2],'r','LineWidth',2.5,'HandleVisibility','off');
    hold on % plot TF
    plot(ES.Time(ES.trialStartsEnds(lastPertTrial,1):ES.trialStartsEnds(lastPertTrial,2))-ES.Time(ES.trialStartsEnds(lastPertTrial,1)),...
         ES.TF(ES.trialStartsEnds(lastPertTrial,1):ES.trialStartsEnds(lastPertTrial,2)),...
         'k','LineWidth',1.5)
    hold on % plot contrast
    plot(ES.Time(ES.trialStartsEnds(j,1):ES.trialStartsEnds(j,2))-ES.Time(ES.trialStartsEnds(j,1)),...
         ES.Contrast(ES.trialStartsEnds(j,1):ES.trialStartsEnds(j,2)),...
         'b','LineWidth',1.5)
    xlim([-0.5 ES.ACInfo.trialStartsEnds(j,2)-ES.ACInfo.trialStartsEnds(j,1)+0.5])
    ylabel('contrast & TF'); set(gca,'TickDir','out','XTickLabel',[]);
    box off
    legend({'TF','Contrast'})
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(6,1,2) % plot mouse speed
    hold off
    plot([ES.ACInfo.trialStartsEnds(j,2)-ES.ACInfo.trialStartsEnds(j,1) ES.ACInfo.trialStartsEnds(j,2)-ES.ACInfo.trialStartsEnds(j,1)],...
        [0 ceil(max(mean(SpeedTrial(find(ES.PertTrial==1),:),1)))],'r','LineWidth',2.5,'HandleVisibility','off');
    hold on
    plot([ES.ACInfo.trialStartsEnds(j,1)-ES.ACInfo.trialStartsEnds(j,1) ES.ACInfo.trialStartsEnds(j,1)-ES.ACInfo.trialStartsEnds(j,1)],...
        [0 ceil(max(mean(SpeedTrial(find(ES.PertTrial==1),:),1)))],'r','LineWidth',2.5,'HandleVisibility','off');
    hold on
    plot([ES.ACInfo.PertStartsEnds(end,2)-ES.ACInfo.trialStartsEnds(lastPertTrial,1) ES.ACInfo.PertStartsEnds(end,2)-ES.ACInfo.trialStartsEnds(lastPertTrial,1)],...
        [0 ceil(max(mean(SpeedTrial(find(ES.PertTrial==1),:),1)))],'k','LineWidth',2,'HandleVisibility','off');
    hold on
    plot([ES.ACInfo.PertStartsEnds(end,1)-ES.ACInfo.trialStartsEnds(lastPertTrial,1) ES.ACInfo.PertStartsEnds(end,1)-ES.ACInfo.trialStartsEnds(lastPertTrial,1)],...
        [0 ceil(max(mean(SpeedTrial(find(ES.PertTrial==1),:),1)))],'k','LineWidth',2,'HandleVisibility','off');
    hold on % plot speed w/o perturbation
    plot(linspace(-0.5,(ES.ACInfo.trialStartsEnds(j,2)+0.5)-(ES.ACInfo.trialStartsEnds(j,1)-0.5),size(SpeedTrial,2)),...
         mean(SpeedTrial(find(ES.PertTrial==0),:),1),...
         'Color',[0.6 0.6 0.6],'LineWidth',1.5)
    hold on % plot speed w perturbation
    plot(linspace(-0.5,(ES.ACInfo.trialStartsEnds(j,2)+0.5)-(ES.ACInfo.trialStartsEnds(j,1)-0.5),size(SpeedTrial,2)),...
         mean(SpeedTrial(find(ES.PertTrial==1),:),1),...
         'Color',[0 0 0],'LineWidth',1.5)
    xlim([-0.5 ES.ACInfo.trialStartsEnds(j,2)-ES.ACInfo.trialStartsEnds(j,1)+0.5])
    ylabel('speed (cm/s)'); set(gca,'TickDir','out','XTickLabel',[]);
    box off
    legend({'no pert','with pert.'})
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(6,1,3) % plot avg firing rate
    hold off
    plot([ES.ACInfo.trialStartsEnds(j,2)-ES.ACInfo.trialStartsEnds(j,1) ES.ACInfo.trialStartsEnds(j,2)-ES.ACInfo.trialStartsEnds(j,1)],...
        [0 ceil(max(smooth(mean(SpikeCount,1)*60,9)))],'k','LineWidth',2,'HandleVisibility','off');
    hold on
    plot([ES.ACInfo.trialStartsEnds(j,1)-ES.ACInfo.trialStartsEnds(j,1) ES.ACInfo.trialStartsEnds(j,1)-ES.ACInfo.trialStartsEnds(j,1)],...
        [0 ceil(max(smooth(mean(SpikeCount,1)*60,9)))],'k','LineWidth',2,'HandleVisibility','off');
    hold on
    plot([ES.ACInfo.PertStartsEnds(end,2)-ES.ACInfo.trialStartsEnds(lastPertTrial,1) ES.ACInfo.PertStartsEnds(end,2)-ES.ACInfo.trialStartsEnds(lastPertTrial,1)],...
        [0 ceil(max(smooth(mean(SpikeCount,1)*60,9)))],'k','LineWidth',2,'HandleVisibility','off');
    hold on
    plot([ES.ACInfo.PertStartsEnds(end,1)-ES.ACInfo.trialStartsEnds(lastPertTrial,1) ES.ACInfo.PertStartsEnds(end,1)-ES.ACInfo.trialStartsEnds(lastPertTrial,1)],...
        [0 ceil(max(smooth(mean(SpikeCount,1)*60,9)))],'k','LineWidth',2,'HandleVisibility','off');
    hold on
    plot(linspace(-0.5,(ES.ACInfo.trialStartsEnds(j,2)+0.5)-(ES.ACInfo.trialStartsEnds(j,1)-0.5),size(SpikeCount,2)), ...
        smooth(mean(SpikeCount(find(ES.PertTrial==0),:),1)*60,9),'Color',[0.6 0.6 0.6],'LineWidth',2)
    hold on
    plot(linspace(-0.5,(ES.ACInfo.trialStartsEnds(j,2)+0.5)-(ES.ACInfo.trialStartsEnds(j,1)-0.5),size(SpikeCount,2)), ...
        smooth(mean(SpikeCount(find(ES.PertTrial==1),:),1)*60,9),'Color',[0 0 0],'LineWidth',2)
    xlim([-0.5 ES.ACInfo.trialStartsEnds(j,2)-ES.ACInfo.trialStartsEnds(j,1)+0.5])
    ylabel('FR (Hz)'); set(gca,'TickDir','out','XTickLabel',[])
    box off
    legend({'no pert','with pert.'})
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pattern = 'Matlab';
    StartLetter = strfind(FileName{1},pattern);
    DetailCell = {[FileName{1}(StartLetter+length(pattern)+1:StartLetter+length(pattern)+13) ', ' FileName{1}(end-22:end-4)] , ...
                   char(strcat('shank=', num2str(ES.SpikeInfo{3,i+1}), ', channel=', num2str(ES.SpikeInfo{4,i+1}), ... 
                   ', ID=', num2str(ES.SpikeInfo{5,i+1}), ', quality=', ES.SpikeInfo{2,i+1}))};
    a = axes;
    t1 = title(DetailCell);
    a.Visible = 'off'; % set(a,'Visible','off');
    t1.Visible = 'on'; % set(t1,'Visible','on');
    
    % CORRECT INFO FOR FIGURE TITLE AND MAKE IT UNIVERSAL
    % check folder if exists
    pattern = 'Matlab';
    StartLetter = strfind(FileName{1},pattern);
    if exist([FileName{1}(1:StartLetter-1) 'Figures/'], 'dir')
        saveas(gcf,[FileName{1}(1:StartLetter-1) 'Figures' filesep ...
                    'PSTH_Recording_' num2str(RecordingOI) '_cellID_' num2str(ES.SpikeInfo{5,i+1}) '_pert'],'png')
    else
        mkdir([FileName{1}(1:StartLetter-1) 'Figures'])
        % save plot in the relevant folder
        saveas(gcf,[FileName{1}(1:StartLetter-1) 'Figures' filesep ...
                    'PSTH_Recording_' num2str(RecordingOI) '_cellID_' num2str(ES.SpikeInfo{5,i+1}) '_pert'],'png')
    end
    close all

end