%% Tomaso Muzzu - UCL - 09/01/2018

% function to plot basic PSTH to visual responses looking only at a
% specific angle for all cells


function PSTH_AllCells_pert_angle(Array,SpikesTrial,ES,FileName,RecordingOI,Direction, RunningThreshold, MaxSpeed,StimulusParams)

for i = 1:size(SpikesTrial,1)
    clear SingletrialSpikes SpikeCount SpeedTrial
    if i == 1
        figure('Renderer', 'painters', 'Position', [10 10 500 1000])
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(6,1,4:6) % plot single spikes
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
        if ES.trialStartsEnds(j,1)-30>0
            SpeedTrial_temp = ES.MouseSpeed(ES.trialStartsEnds(j,1)-30:ES.trialStartsEnds(j,2)+30);
            SpeedTrial(j,:) = SpeedTrial_temp(1:min(ES.trialStartsEnds(:,2)-ES.trialStartsEnds(:,1))+60);
        end
    end
    AllSpikeCount{i} = SpikeCount;
    if i == size(SpikesTrial,1)
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
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if i == 1
        subplot(6,1,1) % plot contrast and TF
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
        legend({'TF','Contrast'},'Location','northwest')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        subplot(6,1,2) % plot mouse speed
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
        TimePoints  = linspace(StimulusParams(1,1),StimulusParams(end,1),length(Array.Speed_in_trial{1})); 
        Speeds = Array.Speed_in_trial(find(Array.PerturbationON==0 & Array.Direction==Direction & Array.Avg_Speed4Rank<MaxSpeed));
        N_traces(1) = length(Speeds);
        MeanSpeed(1,:) = mean(cell2mat(Speeds),1);
        plot(TimePoints, MeanSpeed, 'Color',[0.6 0.6 0.6],'LineWidth',1.5)
        hold on % plot speed w perturbation
        clear Speeds 
        Speeds = Array.Speed_in_trial(find(Array.PerturbationON==1 & Array.Direction==Direction & Array.Avg_Speed4Rank<MaxSpeed));
        N_traces(2) = length(Speeds);
        MeanSpeed(2,:) = mean(cell2mat(Speeds),1);
        plot(TimePoints, MeanSpeed(2,:), 'Color',[0 0 0],'LineWidth',1.5)
        xlim([-0.5 ES.ACInfo.trialStartsEnds(j,2)-ES.ACInfo.trialStartsEnds(j,1)+0.5])
        ylim = ([0 min(max(MeanSpeed(1),MeanSpeed(2))+5,40)]);
        ylabel('speed (cm/s)'); set(gca,'TickDir','out','XTickLabel',[]);
        box off
        legend({['no pert (n=' num2str(N_traces(1)) ')'],['with pert. (n=' num2str(N_traces(2)) ')']},'Location','northwest')
        ylim = ([0 min(max(MeanSpeed(1),MeanSpeed(2))+5,40)]);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if i == size(SpikesTrial,1)
        subplot(6,1,3) % plot avg firing rate
        hold off
        clear SpikeCount SpikeCount_Pert_0 SpikeCount_Pert_1
        m = 1; n = 1;
        for k = 1:length(AllSpikeCount)
            SpikeCount_temp =  AllSpikeCount{k};
            SpikeCount_Pert_0(m,:) = mean(SpikeCount_temp(find(ES.PertTrial==0 & ES.OrienTrial==Direction),:),1)*60;
            m = m + 1;
            SpikeCount_Pert_1(n,:) = mean(SpikeCount_temp(find(ES.PertTrial==1 & ES.OrienTrial==Direction),:),1)*60;
            n = n + 1;
        end
        plot(linspace(-0.5,(ES.ACInfo.trialStartsEnds(j,2)+0.5)-(ES.ACInfo.trialStartsEnds(j,1)-0.5),size(SpikeCount_Pert_0,2)), ...
            smooth(mean(SpikeCount_Pert_0(:,:),1),9),'Color',[0.6 0.6 0.6],'LineWidth',2)
        hold on
        plot(linspace(-0.5,(ES.ACInfo.trialStartsEnds(j,2)+0.5)-(ES.ACInfo.trialStartsEnds(j,1)-0.5),size(SpikeCount_Pert_1,2)), ...
            smooth(mean(SpikeCount_Pert_1(:,:),1),9),'Color',[0 0 0],'LineWidth',2)
        hold on
        plot([ES.ACInfo.trialStartsEnds(j,2)-ES.ACInfo.trialStartsEnds(j,1) ES.ACInfo.trialStartsEnds(j,2)-ES.ACInfo.trialStartsEnds(j,1)],...
            [0 ceil(max(smooth(mean(SpikeCount_Pert_1,1),9)))],'r','LineWidth',2,'HandleVisibility','off');
        hold on
        plot([ES.ACInfo.trialStartsEnds(j,1)-ES.ACInfo.trialStartsEnds(j,1) ES.ACInfo.trialStartsEnds(j,1)-ES.ACInfo.trialStartsEnds(j,1)],...
            [0 ceil(max(smooth(mean(SpikeCount_Pert_1,1),9)))],'r','LineWidth',2,'HandleVisibility','off');
        hold on
        plot([ES.ACInfo.PertStartsEnds(end,2)-ES.ACInfo.trialStartsEnds(lastPertTrial,1) ES.ACInfo.PertStartsEnds(end,2)-ES.ACInfo.trialStartsEnds(lastPertTrial,1)],...
            [0 ceil(max(smooth(mean(SpikeCount_Pert_1,1),9)))],'k','LineWidth',2,'HandleVisibility','off');
        hold on
        plot([ES.ACInfo.PertStartsEnds(end,1)-ES.ACInfo.trialStartsEnds(lastPertTrial,1) ES.ACInfo.PertStartsEnds(end,1)-ES.ACInfo.trialStartsEnds(lastPertTrial,1)],...
            [0 ceil(max(smooth(mean(SpikeCount_Pert_1,1),9)))],'k','LineWidth',2,'HandleVisibility','off');
        hold on
        xlim([-0.5 ES.ACInfo.trialStartsEnds(j,2)-ES.ACInfo.trialStartsEnds(j,1)+0.5])
        ylabel('FR (Hz)'); set(gca,'TickDir','out','XTickLabel',[])
        box off
        legend({'no pert','with pert.'},'Location','northwest')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        pattern = 'Matlab';
        StartLetter = strfind(FileName,pattern);
        DetailCell = {[FileName(StartLetter+length(pattern)+1:StartLetter+length(pattern)+13) ', ' FileName(end-22:end-4)] , ...
            char(strcat(['All Cells, Direction=' num2str(Direction)]))};
        a = axes;
        t1 = title(DetailCell);
        a.Visible = 'off'; % set(a,'Visible','off');
        t1.Visible = 'on'; % set(t1,'Visible','on');
        
        % CORRECT INFO FOR FIGURE TITLE AND MAKE IT UNIVERSAL
        % check folder if exists
        pattern = 'Matlab';
        %StartLetter = strfind(FileName{1},pattern);
%         if exist([FileName{1}(1:StartLetter-1) 'Figures/'], 'dir')
%             saveas(gcf,[FileName{1}(1:StartLetter-1) 'Figures' filesep ...
%                 'PSTH_Recording_' num2str(RecordingOI) '_AllCells_pert_' num2str(Direction)],'png')
%         else
%             mkdir([FileName{1}(1:StartLetter-1) 'Figures'])
%             % save plot in the relevant folder
%             saveas(gcf,[FileName{1}(1:StartLetter-1) 'Figures' filesep ...
%                 'PSTH_Recording_' num2str(RecordingOI) '_AllCells_pert_' num2str(Direction)],'png')
%         end
        %close all
    end
end

end