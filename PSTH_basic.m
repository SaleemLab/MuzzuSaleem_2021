%% Tomaso Muzzu - UCL - 09/01/2018

% function to plot basic PSTH to visual responses 

function PSTH_basic(i,SpikesTrial,ES,FileName,RecordingOI)

clear SingletrialSpikes SpikeCount
    
    figure('Renderer', 'painters', 'Position', [10 10 800 1000])
    subplot(5,1,3:5) % plot single spikes
    hold off
    for j = 1:size(SpikesTrial,2)   % scan through the trials
        SingletrialSpikes = SpikesTrial{i,j};
        plot(SingletrialSpikes-ES.ACInfo.trialStartsEnds(j,1),j*ones(length(SingletrialSpikes),1),'.','Color',[0.5 0.5 0.5]);
        hold on
        edges = linspace(ES.ACInfo.trialStartsEnds(j,1)-0.5,ES.ACInfo.trialStartsEnds(j,2)+0.5,((ES.ACInfo.trialStartsEnds(j,2)+0.5)-(ES.ACInfo.trialStartsEnds(j,1)-0.5))*60);
        [SpikeCount(j,1:length(edges)-1), Edges] = histcounts(SingletrialSpikes,edges); 
    end
    plot([ES.ACInfo.trialStartsEnds(j,2)-ES.ACInfo.trialStartsEnds(j,1) ES.ACInfo.trialStartsEnds(j,2)-ES.ACInfo.trialStartsEnds(j,1)],...
        [0 size(SpikesTrial,2)],'k','LineWidth',2);
    hold on
    plot([ES.ACInfo.trialStartsEnds(j,1)-ES.ACInfo.trialStartsEnds(j,1) ES.ACInfo.trialStartsEnds(j,1)-ES.ACInfo.trialStartsEnds(j,1)],...
        [0 size(SpikesTrial,2)],'k','LineWidth',2);
    hold on
    xlim([-0.5 ES.ACInfo.trialStartsEnds(j,2)-ES.ACInfo.trialStartsEnds(j,1)+0.5])
    xlabel('seconds'); ylabel('trial'); set(gca,'TickDir','out')
    box off
    
    subplot(5,1,1) % plot contrast
    hold off
    plot([ES.ACInfo.trialStartsEnds(j,2)-ES.ACInfo.trialStartsEnds(j,1) ES.ACInfo.trialStartsEnds(j,2)-ES.ACInfo.trialStartsEnds(j,1)],...
        [0 1],'k','LineWidth',2);
    hold on
    plot([ES.ACInfo.trialStartsEnds(j,1)-ES.ACInfo.trialStartsEnds(j,1) ES.ACInfo.trialStartsEnds(j,1)-ES.ACInfo.trialStartsEnds(j,1)],...
        [0 1],'k','LineWidth',2);
    hold on
    plot(ES.Time(ES.trialStartsEnds(j,1):ES.trialStartsEnds(j,2))-ES.Time(ES.trialStartsEnds(j,1)),...
         ES.Contrast(ES.trialStartsEnds(j,1):ES.trialStartsEnds(j,2)),...
         'r')
    xlim([-0.5 ES.ACInfo.trialStartsEnds(j,2)-ES.ACInfo.trialStartsEnds(j,1)+0.5])
    ylabel('contrast'); set(gca,'TickDir','out','XTickLabel',[])
    box off
    
    subplot(5,1,2) % plot avg firing rate
    hold off
    plot([ES.ACInfo.trialStartsEnds(j,2)-ES.ACInfo.trialStartsEnds(j,1) ES.ACInfo.trialStartsEnds(j,2)-ES.ACInfo.trialStartsEnds(j,1)],...
        [0 ceil(max(smooth(mean(SpikeCount,1)*60,9)))],'k','LineWidth',2);
    hold on
    plot([ES.ACInfo.trialStartsEnds(j,1)-ES.ACInfo.trialStartsEnds(j,1) ES.ACInfo.trialStartsEnds(j,1)-ES.ACInfo.trialStartsEnds(j,1)],...
        [0 ceil(max(smooth(mean(SpikeCount,1)*60,9)))],'k','LineWidth',2);
    hold on
    plot(linspace(-0.5,(ES.ACInfo.trialStartsEnds(j,2)+0.5)-(ES.ACInfo.trialStartsEnds(j,1)-0.5),size(SpikeCount,2)),smooth(mean(SpikeCount,1)*60,9),'Color',[0.6 0.6 0.6],'LineWidth',2)
    xlim([-0.5 ES.ACInfo.trialStartsEnds(j,2)-ES.ACInfo.trialStartsEnds(j,1)+0.5])
    ylabel('FR (Hz)'); set(gca,'TickDir','out','XTickLabel',[])
    box off
    
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
                    'PSTH_Recording_' num2str(RecordingOI) '_cellID_' num2str(ES.SpikeInfo{5,i+1})],'png')
    else
        mkdir([FileName{1}(1:StartLetter-1) 'Figures'])
        % save plot in the relevant folder
        saveas(gcf,[FileName{1}(1:StartLetter-1) 'Figures' filesep ...
                    ],'png')
    end
    close all

end