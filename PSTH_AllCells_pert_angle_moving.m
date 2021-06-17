%% Tomaso Muzzu - UCL - 10/01/2018

% function to plot basic PSTH to visual responses looking only at a
% specific angle for all cells differing between stationary and moving
% trials


function PSTH_AllCells_pert_angle_moving(Array,SpikesTrial,ES,FileName,RecordingOI,Direction,rth,MaxSpeed,StimulusParams)

for i = 1:size(SpikesTrial,1)
    clear SingletrialSpikes SpikeCount SpeedTrial
    if i == 1
        figure('Renderer', 'painters', 'Position', [10 10 500 1000])
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     subplot(6,1,4:6) % plot single spikes
    trialSide_seconds = 0.5; % take 60 samples before and after trial
    trialSide_samples = trialSide_seconds*round(length(ES.Time)/(max(ES.Time)-min(ES.Time))); % take 60 samples before and after trial
    for j = 1:size(SpikesTrial,2)   % scan through the trials
        SingletrialSpikes = SpikesTrial{i,j};
        if ES.PertTrial(j) == 0 && ES.OrienTrial(j) == Direction
%             plot(SingletrialSpikes-ES.ACInfo.trialStartsEnds(j,1),j*ones(length(SingletrialSpikes),1),'.','Color',[0.6 0.6 0.6]);
        elseif ES.PertTrial(j) == 1 && ES.OrienTrial(j) == Direction
%             plot(SingletrialSpikes-ES.ACInfo.trialStartsEnds(j,1),j*ones(length(SingletrialSpikes),1),'.','Color',[0 0 0]);
               lastPertTrial = j;
        end
        hold on
        edges = linspace(ES.ACInfo.trialStartsEnds(j,1)-trialSide_seconds,ES.ACInfo.trialStartsEnds(j,2)+trialSide_seconds,((ES.ACInfo.trialStartsEnds(j,2)+trialSide_seconds)-(ES.ACInfo.trialStartsEnds(j,1)-trialSide_seconds))*60);
        [SpikeCount(j,1:length(edges)-1), Edges] = histcounts(SingletrialSpikes,edges);
        % save speed for every trial
        SpeedTrial_temp = ES.MouseSpeed(ES.trialStartsEnds(j,1)-trialSide_samples:ES.trialStartsEnds(j,2)+trialSide_samples);
        SpeedTrial(j,:) = SpeedTrial_temp(1:min(ES.trialStartsEnds(:,2)-ES.trialStartsEnds(:,1))+60);
    end
    AllSpikeCount{i} = SpikeCount;
%     if i == 1
%         plot([ES.ACInfo.trialStartsEnds(j,2)-ES.ACInfo.trialStartsEnds(j,1) ES.ACInfo.trialStartsEnds(j,2)-ES.ACInfo.trialStartsEnds(j,1)],...
%             [0 size(SpikesTrial,2)],'r','LineWidth',2.5);
%         hold on
%         plot([ES.ACInfo.trialStartsEnds(j,1)-ES.ACInfo.trialStartsEnds(j,1) ES.ACInfo.trialStartsEnds(j,1)-ES.ACInfo.trialStartsEnds(j,1)],...
%             [0 size(SpikesTrial,2)],'r','LineWidth',2.5);
%         hold on
%         plot([ES.ACInfo.PertStartsEnds(end,2)-ES.ACInfo.trialStartsEnds(max(find(ES.PertTrial==1)),1) ES.ACInfo.PertStartsEnds(end,2)-ES.ACInfo.trialStartsEnds(max(find(ES.PertTrial==1)),1)],...
%             [0 size(SpikesTrial,2)],'k','LineWidth',2);
%         hold on
%         plot([ES.ACInfo.PertStartsEnds(end,1)-ES.ACInfo.trialStartsEnds(max(find(ES.PertTrial==1)),1) ES.ACInfo.PertStartsEnds(end,1)-ES.ACInfo.trialStartsEnds(max(find(ES.PertTrial==1)),1)],...
%             [0 size(SpikesTrial,2)],'k','LineWidth',2);
%         xlim([-0.5 ES.ACInfo.trialStartsEnds(j,2)-ES.ACInfo.trialStartsEnds(j,1)+0.5])
%         xlabel('seconds'); ylabel('trial'); set(gca,'TickDir','out')
%         box off
%     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if i == 1
        subplot(3,1,1) % plot contrast and TF
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
        legend({'TF','Contrast'},'Location','southwest')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        subplot(3,1,2) % plot mouse speed
        plot([ES.ACInfo.trialStartsEnds(j,2)-ES.ACInfo.trialStartsEnds(j,1) ES.ACInfo.trialStartsEnds(j,2)-ES.ACInfo.trialStartsEnds(j,1)],...
            [0 ceil(max(max(SpeedTrial)))],'r','LineWidth',2.5,'HandleVisibility','off');
        hold on
        plot([ES.ACInfo.trialStartsEnds(j,1)-ES.ACInfo.trialStartsEnds(j,1) ES.ACInfo.trialStartsEnds(j,1)-ES.ACInfo.trialStartsEnds(j,1)],...
            [0 ceil(max(max(SpeedTrial)))],'r','LineWidth',2.5,'HandleVisibility','off');
        hold on
        plot([ES.ACInfo.PertStartsEnds(end,2)-ES.ACInfo.trialStartsEnds(max(find(ES.PertTrial==1)),1) ES.ACInfo.PertStartsEnds(end,2)-ES.ACInfo.trialStartsEnds(max(find(ES.PertTrial==1)),1)],...
            [0 ceil(max(max(SpeedTrial)))],'k','LineWidth',2,'HandleVisibility','off');
        hold on
        plot([ES.ACInfo.PertStartsEnds(end,1)-ES.ACInfo.trialStartsEnds(max(find(ES.PertTrial==1)),1) ES.ACInfo.PertStartsEnds(end,1)-ES.ACInfo.trialStartsEnds(max(find(ES.PertTrial==1)),1)],...
            [0 ceil(max(max(SpeedTrial)))],'k','LineWidth',2,'HandleVisibility','off');
        hold on % plot speed w/o perturbation
        TimePoints  = linspace(StimulusParams(1,1),StimulusParams(end,1),length(Array.Speed_in_trial{1}));
        try
            clear Speeds
            Speeds=Array.Speed_in_trial(find(ES.PertTrial==0 & ES.OrienTrial==Direction & (Array.Avg_Speed4Rank>rth & Array.Avg_Speed4Rank<MaxSpeed)'));
            MeanSpeed = mean(cell2mat(Speeds),1);
            plot(TimePoints,MeanSpeed,'Color',[0.6 0.6 0.6],'LineWidth',1.5)
        catch
            plot(TimePoints,-1*ones(length(TimePoints),1),'Color',[0.6 0.6 0.6],'LineWidth',1.5)
        end
        try
            clear Speeds
            Speeds = Array.Speed_in_trial(find(ES.PertTrial==0 & ES.OrienTrial==Direction & (Array.Avg_Speed4Rank<rth & Array.Avg_Speed4Rank<MaxSpeed)'));
            MeanSpeed = mean(cell2mat(Speeds),1);
            plot(TimePoints,MeanSpeed,':','Color',[0.6 0.6 0.6],'LineWidth',1.5)
        catch
            plot(TimePoints,-1*ones(length(TimePoints),1),':','Color',[0.6 0.6 0.6],'LineWidth',1.5)
        end
        hold on % plot speed w perturbation
        try
            clear Speeds
            Speeds = Array.Speed_in_trial(find(ES.PertTrial==1 & ES.OrienTrial==Direction & (Array.Avg_Speed4Rank>rth & Array.Avg_Speed4Rank<MaxSpeed)'));
            MeanSpeed = mean(cell2mat(Speeds),1);
            plot(TimePoints,MeanSpeed,'Color',[0 0 0],'LineWidth',1.5)
        catch
            plot(TimePoints,-1*ones(length(TimePoints),1),'Color',[0.6 0.6 0.6],'LineWidth',1.5)
        end
        try
            clear Speeds
            Speeds = Array.Speed_in_trial(find(ES.PertTrial==1 & ES.OrienTrial==Direction & (Array.Avg_Speed4Rank<rth & Array.Avg_Speed4Rank<MaxSpeed)'));
            MeanSpeed = mean(cell2mat(Speeds),1);
            plot(TimePoints,MeanSpeed,':','Color',[0 0 0],'LineWidth',1.5)
        catch
            plot(TimePoints,-1*ones(length(TimePoints),1),':','Color',[0 0 0],'LineWidth',1.5)
        end
        xlim([-0.5 ES.ACInfo.trialStartsEnds(j,2)-ES.ACInfo.trialStartsEnds(j,1)+0.5])
        ylabel('speed (cm/s)'); set(gca,'TickDir','out','XTickLabel',[]);
        ylim([-0.5  max(max(cell2mat(Array.Speed_in_trial)))+3])
        box off
        N(1) = length(find(Array.PerturbationON==0 & Array.Direction==Direction & (Array.Avg_Speed4Rank>rth & Array.Avg_Speed4Rank<MaxSpeed))); % no perturbation running
        N(2) = length(find(Array.PerturbationON==0 & Array.Direction==Direction & (Array.Avg_Speed4Rank<rth & Array.Avg_Speed4Rank<MaxSpeed))); % no perturbation still
        N(3) = length(find(Array.PerturbationON==1 & Array.Direction==Direction & (Array.Avg_Speed4Rank>rth & Array.Avg_Speed4Rank<MaxSpeed))); % perturbation running
        N(4) = length(find(Array.PerturbationON==1 & Array.Direction==Direction & (Array.Avg_Speed4Rank<rth & Array.Avg_Speed4Rank<MaxSpeed))); % perturbation still
        legend({['no pert - moving (' num2str(N(1)) ')'], ...
                ['no pert - still (' num2str(N(2)) ')'], ...
                ['with pert. - moving (' num2str(N(3)) ')'], ...
                ['with pert. - still (' num2str(N(4)) ')']},...
                'Location','northwest')
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if i == size(SpikesTrial,1)
        subplot(3,1,3) % plot avg firing rate
        hold off
        clear SpikeCount SpikeCount_Pert_0_0 SpikeCount_Pert_0_1 SpikeCount_Pert_1_0 SpikeCount_Pert_1_1
        m1=1; m2=1; n1=1; n2=1;
        for k = 1:length(AllSpikeCount)
            SpikeCount_temp =  AllSpikeCount{k};
            SpikeCount_Pert_0_0(m1,:) = mean(SpikeCount_temp(find(ES.PertTrial==0 & ES.OrienTrial==Direction ...
                & (Array.Avg_Speed4Rank<rth & Array.Avg_Speed4Rank<MaxSpeed)'),:),1)*60; m1 = m1 + 1;
            SpikeCount_Pert_0_1(m2,:) = mean(SpikeCount_temp(find(ES.PertTrial==0 & ES.OrienTrial==Direction ...
                & (Array.Avg_Speed4Rank>rth & Array.Avg_Speed4Rank<MaxSpeed)'),:),1)*60; m2 = m2 + 1;
            SpikeCount_Pert_1_0(n1,:) = mean(SpikeCount_temp(find(ES.PertTrial==1 & ES.OrienTrial==Direction ...
                & (Array.Avg_Speed4Rank<rth & Array.Avg_Speed4Rank<MaxSpeed)'),:),1)*60; n1 = n1 + 1;
            SpikeCount_Pert_1_1(n2,:) = mean(SpikeCount_temp(find(ES.PertTrial==1 & ES.OrienTrial==Direction ...
                & (Array.Avg_Speed4Rank>rth & Array.Avg_Speed4Rank<MaxSpeed)'),:),1)*60; n2 = n2 + 1;
        end
        TimePoints_s = linspace(-0.5,(ES.ACInfo.trialStartsEnds(j,2)+0.5)-(ES.ACInfo.trialStartsEnds(j,1)-0.5),size(SpikeCount_Pert_0_0,2));
        plot(TimePoints_s, smooth(mean(SpikeCount_Pert_0_1(:,:),1),9),'Color',[0.6 0.6 0.6],'LineWidth',2)
        hold on
        plot(TimePoints_s, smooth(mean(SpikeCount_Pert_0_0(:,:),1),9),':','Color',[0.6 0.6 0.6],'LineWidth',2)
        hold on
        plot(TimePoints_s, smooth(mean(SpikeCount_Pert_1_1(:,:),1),9),'Color',[0 0 0],'LineWidth',2)
        hold on
        plot(TimePoints_s, smooth(mean(SpikeCount_Pert_1_0(:,:),1),9),':','Color',[0 0 0],'LineWidth',2)
        hold on
        plot([ES.ACInfo.trialStartsEnds(j,2)-ES.ACInfo.trialStartsEnds(j,1) ES.ACInfo.trialStartsEnds(j,2)-ES.ACInfo.trialStartsEnds(j,1)],...
            [0 ceil(max(smooth(mean(SpikeCount_Pert_0_0,1),9)))],'r','LineWidth',2,'HandleVisibility','off');
        hold on
        plot([ES.ACInfo.trialStartsEnds(j,1)-ES.ACInfo.trialStartsEnds(j,1) ES.ACInfo.trialStartsEnds(j,1)-ES.ACInfo.trialStartsEnds(j,1)],...
            [0 ceil(max(smooth(mean(SpikeCount_Pert_0_0,1),9)))],'r','LineWidth',2,'HandleVisibility','off');
        hold on
        plot([ES.ACInfo.PertStartsEnds(end,2)-ES.ACInfo.trialStartsEnds(max(find(ES.PertTrial==1)),1) ES.ACInfo.PertStartsEnds(end,2)-ES.ACInfo.trialStartsEnds(max(find(ES.PertTrial==1)),1)],...
            [0 ceil(max(smooth(mean(SpikeCount_Pert_0_0,1),9)))],'k','LineWidth',2,'HandleVisibility','off');
        hold on
        plot([ES.ACInfo.PertStartsEnds(end,1)-ES.ACInfo.trialStartsEnds(max(find(ES.PertTrial==1)),1) ES.ACInfo.PertStartsEnds(end,1)-ES.ACInfo.trialStartsEnds(max(find(ES.PertTrial==1)),1)],...
            [0 ceil(max(smooth(mean(SpikeCount_Pert_0_0,1),9)))],'k','LineWidth',2,'HandleVisibility','off');
        hold on
        xlim([-0.5 ES.ACInfo.trialStartsEnds(j,2)-ES.ACInfo.trialStartsEnds(j,1)+0.5])
        ylabel('FR (Hz)'); set(gca,'TickDir','out')
        xlabel('seconds');
        box off
%         legend({['no pert - moving (' num2str(N(1)) ')'], ...
%                 ['no pert - still (' num2str(N(2)) ')'], ...
%                 ['with pert. - moving (' num2str(N(3)) ')'], ...
%                 ['with pert. - still (' num2str(N(4)) ')']},...
%                 'Location','northwest')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        pattern = 'Matlab';
        StartLetter = strfind(FileName,pattern);
        nrTrialsPertDir = length(find(ES.PertTrial==1 & ES.OrienTrial==Direction));
        nrTrialsNPertDir = length(find(ES.PertTrial==0 & ES.OrienTrial==Direction));
        DetailCell = {[FileName(StartLetter+length(pattern)+1:StartLetter+length(pattern)+13) ', ' FileName(end-22:end-4)] , ...
            char(strcat(['All Cells, Direction=' num2str(Direction) ', n_p=' num2str(nrTrialsPertDir) ', n_n_p=' num2str(nrTrialsNPertDir)]))};
        a = axes;
        t1 = title(DetailCell);
        a.Visible = 'off'; % set(a,'Visible','off');
        t1.Visible = 'on'; % set(t1,'Visible','on');
        
        % CORRECT INFO FOR FIGURE TITLE AND MAKE IT UNIVERSAL
        % check folder if exists
%         pattern = 'Matlab';
%         StartLetter = strfind(FileName{1},pattern);
%         if exist([FileName{1}(1:StartLetter-1) 'Figures/'], 'dir')
%             saveas(gcf,[FileName{1}(1:StartLetter-1) 'Figures' filesep ...
%                 'PSTH_Recording_' num2str(RecordingOI) '_AllCells_pert_' num2str(Direction)],'png')
%         else
%             mkdir([FileName{1}(1:StartLetter-1) 'Figures'])
%             % save plot in the relevant folder
%             saveas(gcf,[FileName{1}(1:StartLetter-1) 'Figures' filesep ...
%                 'PSTH_Recording_' num2str(RecordingOI) '_AllCells_pert_' num2str(Direction)],'png')
%         end
%         close all
    end
end

end