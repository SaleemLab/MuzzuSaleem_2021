%% Tomaso Muzzu - UCL - 10/01/2018

% function to plot basic PSTH to visual responses looking only at a
% specific angle for all cells differing between stationary and moving
% trials


function SpeedPlot_pert_angle(SpikesTrial,ES,FileName,RecordingOI,Direction,SpeedThreshold,PertStatus)

clear SingletrialSpikes SpikeCount SpeedTrial
figure('Renderer', 'painters', 'Position', [-1930 0 950 1000])

% Gaussian filter
BonvisionFR = (ES.Time(end)-ES.Time(1))/length(ES.Time); 
Sigma = 0.2/BonvisionFR; % standard deviation in number of samples (converted from time in seconds)
Width = 3*Sigma; % convert size from seconds to number of samples
x_g = linspace(-Width/2, Width/2, Width);
gaussFilter = exp(-x_g.^2/(2*Sigma^2));
gaussFilter_ = gaussFilter / sum (gaussFilter); % normalize
SmoothedSpeed = conv(ES.MouseSpeed, gaussFilter_, 'same');

trialSide_seconds = 1; % take 60 samples before and after trial
trialSide_samples = trialSide_seconds*round(length(ES.Time)/(max(ES.Time)-min(ES.Time))); % take 60 samples before and after trial
for j = 1:size(ES.trialStartsEnds,1)   % scan through the trials
    % save speed for every trial
    SpeedTrial_temp = SmoothedSpeed(ES.trialStartsEnds(j,1)-trialSide_samples:ES.trialStartsEnds(j,2)+trialSide_samples);
    SpeedTrial(j,:) = SpeedTrial_temp(1:min(ES.trialStartsEnds(:,2)-ES.trialStartsEnds(:,1))+trialSide_samples);
end
% rank trials by avg running speed
[avg_speeds ind_] = sort(mean(SpeedTrial(:,:),2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(6,1,2) % plot contrast and TF
hold off
plot([ES.ACInfo.trialStartsEnds(j,2)-ES.ACInfo.trialStartsEnds(j,1) ES.ACInfo.trialStartsEnds(j,2)-ES.ACInfo.trialStartsEnds(j,1)],...
    [0 3.2],'r','LineWidth',2.5,'HandleVisibility','off');
hold on
plot([ES.ACInfo.trialStartsEnds(j,1)-ES.ACInfo.trialStartsEnds(j,1) ES.ACInfo.trialStartsEnds(j,1)-ES.ACInfo.trialStartsEnds(j,1)],...
    [0 3.2],'r','LineWidth',2.5,'HandleVisibility','off');
hold on % plot TF
plot(ES.Time(ES.trialStartsEnds(max(find(ES.PertTrial==1)),1):ES.trialStartsEnds(max(find(ES.PertTrial==1)),2))-ES.Time(ES.trialStartsEnds(max(find(ES.PertTrial==1)),1)),...
    ES.TF(ES.trialStartsEnds(max(find(ES.PertTrial==1)),1):ES.trialStartsEnds(max(find(ES.PertTrial==1)),2)),...
    'k','LineWidth',1.5)
hold on % plot contrast
plot(ES.Time(ES.trialStartsEnds(j,1):ES.trialStartsEnds(j,2))-ES.Time(ES.trialStartsEnds(j,1)),...
    ES.Contrast(ES.trialStartsEnds(j,1):ES.trialStartsEnds(j,2)),...
    'b','LineWidth',1.5)
plot([ES.Time(ES.PertStartsEnds(end,1))-ES.Time(ES.trialStartsEnds(max(find(ES.PertTrial==1)),1)) ES.Time(ES.PertStartsEnds(end,2))-ES.Time(ES.trialStartsEnds(max(find(ES.PertTrial==1)),1))],...
    [max(ES.TF) max(ES.TF)],':k','LineWidth',1.5,'HandleVisibility','off');
xlim([-trialSide_seconds ES.ACInfo.trialStartsEnds(j,2)-ES.ACInfo.trialStartsEnds(j,1)+trialSide_seconds]); ylim([0 3.2]);
ylabel('contrast & TF'); set(gca,'TickDir','out','YTick',[1:3],'YTickLabel',[1:3],'XTickLabels',[]);
box off
legend({'TF','Contrast'})


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(6,1,1) % plot speed entire session
plot(ES.Time(2:end),SmoothedSpeed)
xlabel('time'); ylabel('Speed');
if exist('SpeedThreshold','var')
    RunningThreshold = 1;
    PercRunning = length(SmoothedSpeed(find(SmoothedSpeed>RunningThreshold)))/length(SmoothedSpeed)*100;
    RunningThreshold = SpeedThreshold;
else
    RunningThreshold = 1;
    PercRunning = length(SmoothedSpeed(find(SmoothedSpeed>RunningThreshold)))/length(SmoothedSpeed)*100;
    RunningThreshold = -10; % dummy value for the other plot: show all trials.
end


ylabel('speed (cm/s)'); set(gca,'TickDir','out');
ylim([min(SmoothedSpeed) max(SmoothedSpeed)])
title(['avg speed=' num2str(mean(SmoothedSpeed),3) 'cm/s, ' num2str(PercRunning,3) '% running'])
box off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist('PertStatus','var')
    subplot(6,1,3:6) % plot speed of trials
    hold off
    runningtrials = 0; k = 1;
    clear TempAngleSpeed runningtrials
    if exist('Direction','var')
        TempDirecton = Direction;
    else
        TempDirecton = sort(unique(ES.OrienTrial));
    end
    for Dir = TempDirecton
        textON = 1;
        counter = 1;
        clear TempAngleSpeed
        for i = 1:size(SpeedTrial,1)
            if mean(SpeedTrial(ind_(i)))>RunningThreshold
                if ES.OrienTrial(ind_(i))==Dir && ES.PertTrial(ind_(i))==PertStatus
                    hold on
                    scatter(linspace(-trialSide_seconds,(ES.ACInfo.trialStartsEnds(j,2)+trialSide_seconds)-(ES.ACInfo.trialStartsEnds(j,1)-trialSide_seconds),size(SpeedTrial,2)),...
                        k*ones(length(SpeedTrial(ind_(i),:)),1), ...
                        [], SpeedTrial(ind_(i),:),'filled');
                    if textON == 1
                        text((ES.ACInfo.trialStartsEnds(j,2)+trialSide_seconds)-(ES.ACInfo.trialStartsEnds(j,1)-trialSide_seconds)+0.1,k,num2str(Dir))
                        textON = 0;
                    end
                    TempAngleSpeed(counter,:)=SpeedTrial(ind_(i),:);
                    counter = counter + 1;
                    runningtrials(k) = i;
                    k = k+1;                    
                end
            end
        end 
    end
    colorbar
    xlim([-trialSide_seconds ES.ACInfo.trialStartsEnds(j,2)-ES.ACInfo.trialStartsEnds(j,1)+trialSide_seconds+2])
    ylabel('trial'); set(gca,'TickDir','out');
    xlabel('seconds');
    if exist('PertStatus','var')
        title(['Trials perturbation status=' num2str(PertStatus) ', (n=' num2str(k) '), running thres=' num2str(RunningThreshold) ', direction=' num2str(TempDirecton)])
    elseif ~exist('PertStatus','var')
        title(['All trials, (n=' num2str(k) '), running thres=' num2str(RunningThreshold) ', direction=' num2str(TempDirecton)])
    end 
        plot([ES.ACInfo.trialStartsEnds(j,2)-ES.ACInfo.trialStartsEnds(j,1) ES.ACInfo.trialStartsEnds(j,2)-ES.ACInfo.trialStartsEnds(j,1)],...
            [min(min(SmoothedSpeed)) k],'r','LineWidth',2.5,'HandleVisibility','off');
        hold on
        plot([ES.ACInfo.trialStartsEnds(j,1)-ES.ACInfo.trialStartsEnds(j,1) ES.ACInfo.trialStartsEnds(j,1)-ES.ACInfo.trialStartsEnds(j,1)],...
            [min(min(SmoothedSpeed)) k],'r','LineWidth',2.5,'HandleVisibility','off');
        hold on
        plot([ES.Time(ES.PertStartsEnds(end,2))-ES.Time(ES.trialStartsEnds(max(find(ES.PertTrial(:)==1)))) ES.Time(ES.PertStartsEnds(end,2))-ES.Time(ES.trialStartsEnds(max(find(ES.PertTrial(:)==1))))],...
            [min(min(SmoothedSpeed)) k],'k','LineWidth',1.5,'HandleVisibility','off');
        hold on
        plot([ES.Time(ES.PertStartsEnds(end,1))-ES.Time(ES.trialStartsEnds(max(find(ES.PertTrial(:)==1)))) ES.Time(ES.PertStartsEnds(end,1))-ES.Time(ES.trialStartsEnds(max(find(ES.PertTrial(:)==1))))],...
             [min(min(SmoothedSpeed)) k],'k','LineWidth',1.5,'HandleVisibility','off');
        ylim([0 k]); xlim([-trialSide_seconds ES.ACInfo.trialStartsEnds(j,2)-ES.ACInfo.trialStartsEnds(j,1)+trialSide_seconds]);
        box off
else
    subplot(6,1,3:6) % plot speed of trials
    hold off
    runningtrials = 0; k = 1;
    clear TempAngleSpeed runningtrials
    if exist('Direction','var')
        TempDirecton = Direction;
    else
        TempDirecton = sort(unique(ES.OrienTrial));
    end
    for Dir = TempDirecton
        textON = 1;
        counter = 1;
        clear TempAngleSpeed
        for i = 1:size(SpeedTrial,1)
            if mean(SpeedTrial(ind_(i)))>RunningThreshold
                if ES.OrienTrial(ind_(i))==Dir
                    hold on
                    scatter(linspace(-trialSide_seconds,(ES.ACInfo.trialStartsEnds(j,2)+trialSide_seconds)-(ES.ACInfo.trialStartsEnds(j,1)-trialSide_seconds),size(SpeedTrial,2)),...
                        k*ones(length(SpeedTrial(ind_(i),:)),1), ...
                        [], SpeedTrial(ind_(i),:),'filled');
                    if textON == 1
                        text((ES.ACInfo.trialStartsEnds(j,2)+trialSide_seconds)-(ES.ACInfo.trialStartsEnds(j,1)-trialSide_seconds)+0.1,k,num2str(Dir))
                        textON = 0;
                    end
                    TempAngleSpeed(counter,:)=SpeedTrial(ind_(i),:);
                    counter = counter + 1;
                    runningtrials(k) = i;
                    k = k+1;
                end
            end
        end
    end
    colorbar
    xlim([-trialSide_seconds ES.ACInfo.trialStartsEnds(j,2)-ES.ACInfo.trialStartsEnds(j,1)+trialSide_seconds+2])
    ylabel('trial'); set(gca,'TickDir','out');
    xlabel('seconds');
    if exist('PertStatus','var')
        title(['Trials perturbation status=' num2str(PertStatus) ', (n=' num2str(length(runningtrials)) '), running thres=' num2str(RunningThreshold) ', direction=' num2str(TempDirecton)])
    elseif ~exist('PertStatus','var')
        title(['All trials, (n=' num2str(k) '), running thres=' num2str(RunningThreshold) ', direction=' num2str(TempDirecton)])
    else
        
        plot([ES.ACInfo.trialStartsEnds(j,2)-ES.ACInfo.trialStartsEnds(j,1) ES.ACInfo.trialStartsEnds(j,2)-ES.ACInfo.trialStartsEnds(j,1)],...
            [min(min(SmoothedSpeed)) k],'r','LineWidth',2.5,'HandleVisibility','off');
        hold on
        plot([ES.ACInfo.trialStartsEnds(j,1)-ES.ACInfo.trialStartsEnds(j,1) ES.ACInfo.trialStartsEnds(j,1)-ES.ACInfo.trialStartsEnds(j,1)],...
            [min(min(SmoothedSpeed)) k],'r','LineWidth',2.5,'HandleVisibility','off');
        ylim([0 k]); xlim([-trialSide_seconds ES.ACInfo.trialStartsEnds(j,2)-ES.ACInfo.trialStartsEnds(j,1)+trialSide_seconds]);
        box off
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pattern = 'Matlab';
StartLetter = strfind(FileName{1},pattern);
DetailCell = {[FileName{1}(StartLetter+length(pattern)+1:StartLetter+length(pattern)+13) ', ' FileName{1}(end-22:end-4)]};
% a = axes;
% t1 = title(DetailCell);
% a.Visible = 'off'; % set(a,'Visible','off');
% t1.Visible = 'on'; % set(t1,'Visible','on');

% % % CORRECT INFO FOR FIGURE TITLE AND MAKE IT UNIVERSAL
% check folder if exists
pattern = 'Matlab';
StartLetter = strfind(FileName{1},pattern);
DirInfo = {'AllDir', num2str(Direction)};
[DirSave ind_dir] = min([length(DirInfo{1}),length(DirInfo{2})]);
if exist([FileName{1}(1:StartLetter-1) 'Figures'], 'dir')
    if exist('PertStatus','var')
        saveas(gcf,[FileName{1}(1:StartLetter-1) 'Figures' filesep ...
            'Speed_' num2str(RecordingOI) '_Direction_' DirInfo{ind_dir} '_RunTH_' num2str(RunningThreshold) '_Pert_' num2str(PertStatus)],'png')
    else
        saveas(gcf,[FileName{1}(1:StartLetter-1) 'Figures' filesep ...
            'Speed_' num2str(RecordingOI) '_Direction_' DirInfo{ind_dir} '_RunTH_' num2str(RunningThreshold) '_'],'png')
    end
else
    mkdir([FileName{1}(1:StartLetter-1) 'Figures'])
    % save plot in the relevant folder
    if exist('PertStatus','var')
        saveas(gcf,[FileName{1}(1:StartLetter-1) 'Figures' filesep ...
            'Speed_' num2str(RecordingOI) '_Direction_' DirInfo{ind_dir} '_RunTH_' num2str(RunningThreshold) '_Pert_' num2str(PertStatus)],'png')
    else
        saveas(gcf,[FileName{1}(1:StartLetter-1) 'Figures' filesep ...
            'Speed_' num2str(RecordingOI) '_Direction_' DirInfo{ind_dir} '_RunTH_' num2str(RunningThreshold) '_'],'png')
    end
end
close all

%EOF
end
