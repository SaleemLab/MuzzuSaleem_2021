%% Tomaso Muzzu - UCL - 15/01/2018

% function to prepare data of speed profile of mice during visual stimuli
% presentaiton of grating with changing contrast and possible TF
% perturbation.

function [Array, StimulusParams] = GetSpeedData(ES,varargin)

%Speed_Plot(ES,FileName,RecordingOI,Direction,SpeedThreshold,PertStatus)

pnames = {'Direction','Running_Th','Pert_Status','Rank_mode','Norm_mode','log_speed','AvgSpeedPeriod','Smoothing'};

Dflt_Values = {0:45:315, 10, 2, 'AVGspeed', 'none', 0, [-1 0], 0.6};

[Direction, RunningThreshold, PerStatus, RankingBy, NormMode, logSpeed, AvgSpeedPeriod, GaussFilter_W] ...
    = internal.stats.parseArgs(pnames,Dflt_Values,varargin{:});

% Default inputs:
%%% Direction: 0:45:315; % plot all trials irrespective of the grating
% direction
%%% Running_Th: -10; % running threshold is set at -10 so all trials are
% considered
%%% Pert_Status = 2; % if 1 -> trials with perturbation; if 0 -> trials
% without perturbation
%%% Rank_mode = 'AVGspeed'; % rank the trials by the avg speed during the
% trials. Otherwise: 
% ------ 'TrialStart'= speed at -1s to +1s of trial onset. Period
% modifiable by user with input 'AvgSpeedPeriod'.
% ------ 'PertStart'= speed at -1s to 0s of perturbation onset. This only
% true if considering the trials with perturbation.  Period modifiable by 
% user at input.
%%% NormMode = 'none'; do not normalise the speed. Other options are:
% ------ 'MaxSpeedTrial'= normalise by max speed of trial
% ------ 'MaxSpeedSession'= normalise by max speed of session
% ------ 'z_score'= normalise by z score (X-mean)/std of session
%%% logSpeed = 1; % compute speed in log units. if 0 -> normal units.
%%% Smoothing = 0.2; % width of smoothing window in seconds

% Gaussian filter
BonvisionFR = (ES.Time(end)-ES.Time(1))/length(ES.Time); 
Width = round(GaussFilter_W/BonvisionFR);
Sigma = Width/3; % standard deviation in number of samples (converted from time in seconds)
x_g = linspace(-Width/2, Width/2, Width);
gaussFilter = exp(-x_g.^2/(2*Sigma^2));
gaussFilter_ = gaussFilter / sum (gaussFilter); % normalize
if sum(diff(ES.MousePosition))<0
    SmoothedSpeed = conv(-ES.MouseSpeed, gaussFilter_, 'same');
else 
    SmoothedSpeed = conv(ES.MouseSpeed, gaussFilter_, 'same');
end
SmoothedSpeed(SmoothedSpeed<=0) = 0.01;
% normalisation
switch NormMode 
    case 'none'
        SmoothedSpeed = SmoothedSpeed;
    case 'MaxSpeedSession'
        SmoothedSpeed = SmoothedSpeed/max(SmoothedSpeed);
    case 'z_score'
        SmoothedSpeed = (SmoothedSpeed - mean(SmoothedSpeed))/std(SmoothedSpeed);
    otherwise
        SmoothedSpeed = SmoothedSpeed;
end
% log scale
if logSpeed == 1
    SmoothedSpeed = log(SmoothedSpeed);
end
% save speed for every trial
trialSide_seconds = 1; % take 60 samples before and after trial
trialSide_samples = trialSide_seconds*round(length(ES.Time)/(max(ES.Time)-min(ES.Time))); % take 60 samples before and after trial
clear SpeedTrial_temp SpeedTrial SpeedTrial_PertStartEnds clear StimulusParams
r = 1;
for j = 1:size(ES.trialStartsEnds,1)   % scan through the trials
    %if ES.trialStartsEnds(j,1)-trialSide_samples>0
        SpeedTrial_temp = SmoothedSpeed(max(ES.trialStartsEnds(j,1)-trialSide_samples,1):min(length(SmoothedSpeed),ES.trialStartsEnds(j,2)+trialSide_samples));
        SpeedTrial_temp = interp1(1:length(SpeedTrial_temp),SpeedTrial_temp,1:560);
        if strcmp(NormMode, 'MaxSpeedTrial')
            SpeedTrial_temp = SpeedTrial_temp/max(SpeedTrial_temp);
        end
        SpeedTrial(j,:) = SpeedTrial_temp; clear SpeedTrial_temp
        %     if length(SpeedTrial_temp)>=(min(ES.trialStartsEnds(:,2)-ES.trialStartsEnds(:,1))+2*trialSide_samples)
        %         SpeedTrial(j,:) = SpeedTrial_temp(1:min(ES.trialStartsEnds(:,2)-ES.trialStartsEnds(:,1))+2*trialSide_samples); clear SpeedTrial_temp
        %     else
        %         SpeedTrial(j,:) = SpeedTrial_temp; clear SpeedTrial_temp
        %     end
        % save stimulus parameters: Time, contrast and TF
        Temp_t = ES.Time(max(ES.trialStartsEnds(j,1)-trialSide_samples,1):ES.trialStartsEnds(j,2)+trialSide_samples)-min(ES.Time(max(ES.trialStartsEnds(j,1)-trialSide_samples,1):ES.trialStartsEnds(j,2)+trialSide_samples));
        Temp_c = ES.Contrast(max(ES.trialStartsEnds(j,1)-trialSide_samples,1):ES.trialStartsEnds(j,2)+trialSide_samples);
        Temp_tf = ES.TF(max(ES.trialStartsEnds(j,1)-trialSide_samples,1):ES.trialStartsEnds(j,2)+trialSide_samples);
        StimulusParams{j,1}(:,1) = interp1(1:length(Temp_t), Temp_t, 1:560);
        StimulusParams{j,1}(:,2) = interp1(1:length(Temp_c), Temp_c, 1:560);
        StimulusParams{j,1}(:,3) = interp1(1:length(Temp_tf), Temp_tf, 1:560);
        if isfield(ES,'PertStartsEnds')
            if ES.PertTrial(j) == 1
                SpeedTrial_temp = SmoothedSpeed(ES.PertStartsEnds(r,1)-trialSide_samples:ES.PertStartsEnds(r,2)+trialSide_samples);
                if strcmp(NormMode, 'MaxSpeedTrial')
                    SpeedTrial_temp = SpeedTrial_temp/max(SpeedTrial_temp);
                end
                SpeedTrial_PertStartEnds(r,:) = SpeedTrial_temp(1:min(ES.PertStartsEnds(:,2)-ES.PertStartsEnds(:,1))+2*trialSide_samples)';
                r = r + 1;
                % save stimulus parameters: Time, contrast and TF
%                 Temp_t = ES.Time(ES.trialStartsEnds(j,1)-trialSide_samples:ES.trialStartsEnds(j,2)+trialSide_samples)-min(ES.Time(ES.trialStartsEnds(j,1)-trialSide_samples:ES.trialStartsEnds(j,2)+trialSide_samples));
%                 Temp_c = ES.Contrast(ES.trialStartsEnds(j,1)-trialSide_samples:ES.trialStartsEnds(j,2)+trialSide_samples);
%                 Temp_tf = ES.TF(ES.trialStartsEnds(j,1)-trialSide_samples:ES.trialStartsEnds(j,2)+trialSide_samples);
%                 StimulusParams(:,1) = Temp_t(1:min(ES.trialStartsEnds(:,2)-ES.trialStartsEnds(:,1))+2*trialSide_samples);
%                 StimulusParams(:,2) = Temp_c(1:min(ES.trialStartsEnds(:,2)-ES.trialStartsEnds(:,1))+2*trialSide_samples);
%                 StimulusParams(:,3) = Temp_tf(1:min(ES.trialStartsEnds(:,2)-ES.trialStartsEnds(:,1))+2*trialSide_samples);
            end
        else
            if j == size(ES.trialStartsEnds,1)-1
                % save stimulus parameters: Time, contrast and TF
%                 Temp_t = ES.Time(ES.trialStartsEnds(j,1)-trialSide_samples:ES.trialStartsEnds(j,2)+trialSide_samples)-min(ES.Time(ES.trialStartsEnds(j,1)-trialSide_samples:ES.trialStartsEnds(j,2)+trialSide_samples));
%                 Temp_c = ES.Contrast(ES.trialStartsEnds(j,1)-trialSide_samples:ES.trialStartsEnds(j,2)+trialSide_samples);
%                 Temp_tf = ES.TF(ES.trialStartsEnds(j,1)-trialSide_samples:ES.trialStartsEnds(j,2)+trialSide_samples);
%                 StimulusParams(:,1) = Temp_t(1:min(ES.trialStartsEnds(:,2)-ES.trialStartsEnds(:,1))+2*trialSide_samples);
%                 StimulusParams(:,2) = Temp_c(1:min(ES.trialStartsEnds(:,2)-ES.trialStartsEnds(:,1))+2*trialSide_samples);
%                 StimulusParams(:,3) = Temp_tf(1:min(ES.trialStartsEnds(:,2)-ES.trialStartsEnds(:,1))+2*trialSide_samples);
            end
        end
    %end
end
% ranking mode
clear avg_speeds ind
switch RankingBy 
    case 'AVGspeed' % rank the trials by the avg speed during the trial
        TimeRange = [1 , size(SpeedTrial(:,:),2)];     
    case 'TrialStart' % speed at -1s to 0s
        TimeRange = [max(1,1+trialSide_samples+round((AvgSpeedPeriod(1)/BonvisionFR))), ...
                                        1+trialSide_samples+round((AvgSpeedPeriod(2)/BonvisionFR))];
    case 'PertStart' % speed at -1s to 0s
        TimeRange = [1 , size(SpeedTrial(:,:),2)];
end
MeanSpeedInRange = mean(SpeedTrial(:,TimeRange(1):TimeRange(2)),2);
[sortedSpeed SpeedSort_idx] = sort(MeanSpeedInRange);

TimeRange_Pert =  [min(1,1+trialSide_samples-round((AvgSpeedPeriod(1)/BonvisionFR))), ...
                            1+trialSide_samples+round((AvgSpeedPeriod(2)/BonvisionFR))];
if isfield(ES,'PertStartsEnds')
    MeanSpeedinRange_pert = mean(SpeedTrial_PertStartEnds(:,TimeRange_Pert(1):TimeRange_Pert(2)),2);
    [sortedSpeed_pert SpeedSortPert_idx] = sort(MeanSpeedinRange_pert);
end
    

% save array with info necessary to plot everything
clear Array
cc = 1;
Array = table(1,1,1,1,{ones(1,length(SpeedTrial(1,:)))},1,{1},{1},{1},0,{0},...
        'VariableNames',{'Trial_ID';'Direction';'Trial_rank';'Avg_Speed4Rank';'Speed_in_trial';'PerturbationON'; 'time_in_trial'; 'contrast_in_trial'; 'TF_in_trial'; 'Pert_Rank';'Speed_in_pert'})
for i = 1:size(SpeedTrial,1)
    Array{i,1} = i; % order of display
    Array{i,2} = ES.OrienTrial(i); % grating direction
    Array{i,3} = SpeedSort_idx(i); % order of trials by ranking of choice
    Array{i,4} = MeanSpeedInRange(i); %SpeedSort_idx(i)); % avg speed used for ranking trials
    Array{i,5} = {SpeedTrial(i,:)}; % speed of every trial
    Array{i,7} = {StimulusParams{i,:}(:,1)}; % speed of every trial
    Array{i,8} = {StimulusParams{i,:}(:,2)}; % speed of every trial
    Array{i,9} = {StimulusParams{i,:}(:,3)}; % speed of every trial
    if isfield(ES,'PertStartsEnds')
        Array{i,6} = ES.PertTrial(i); % perturbation
        if ES.PertTrial(i)==1 && cc<=length(SpeedSortPert_idx)
            Array{i,10} = SpeedSortPert_idx(cc);
            Array{i,11} = {SpeedTrial_PertStartEnds(cc,:)};
            cc=cc+1;
        end
    end
end


% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subplot(6,1,2) % plot contrast and TF
% hold off
% plot([ES.ACInfo.trialStartsEnds(j,2)-ES.ACInfo.trialStartsEnds(j,1) ES.ACInfo.trialStartsEnds(j,2)-ES.ACInfo.trialStartsEnds(j,1)],...
%     [0 3.2],'r','LineWidth',2.5,'HandleVisibility','off');
% hold on
% plot([ES.ACInfo.trialStartsEnds(j,1)-ES.ACInfo.trialStartsEnds(j,1) ES.ACInfo.trialStartsEnds(j,1)-ES.ACInfo.trialStartsEnds(j,1)],...
%     [0 3.2],'r','LineWidth',2.5,'HandleVisibility','off');
% hold on % plot TF
% plot(ES.Time(ES.trialStartsEnds(max(find(ES.PertTrial==1)),1):ES.trialStartsEnds(max(find(ES.PertTrial==1)),2))-ES.Time(ES.trialStartsEnds(max(find(ES.PertTrial==1)),1)),...
%     ES.TF(ES.trialStartsEnds(max(find(ES.PertTrial==1)),1):ES.trialStartsEnds(max(find(ES.PertTrial==1)),2)),...
%     'k','LineWidth',1.5)
% hold on % plot contrast
% plot(ES.Time(ES.trialStartsEnds(j,1):ES.trialStartsEnds(j,2))-ES.Time(ES.trialStartsEnds(j,1)),...
%     ES.Contrast(ES.trialStartsEnds(j,1):ES.trialStartsEnds(j,2)),...
%     'b','LineWidth',1.5)
% plot([ES.Time(ES.PertStartsEnds(end,1))-ES.Time(ES.trialStartsEnds(max(find(ES.PertTrial==1)),1)) ES.Time(ES.PertStartsEnds(end,2))-ES.Time(ES.trialStartsEnds(max(find(ES.PertTrial==1)),1))],...
%     [max(ES.TF) max(ES.TF)],':k','LineWidth',1.5,'HandleVisibility','off');
% xlim([-trialSide_seconds ES.ACInfo.trialStartsEnds(j,2)-ES.ACInfo.trialStartsEnds(j,1)+trialSide_seconds]); ylim([0 3.2]);
% ylabel('contrast & TF'); set(gca,'TickDir','out','YTick',[1:3],'YTickLabel',[1:3],'XTickLabels',[]);
% box off
% legend({'TF','Contrast'})
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subplot(6,1,1) % plot speed entire session
% plot(ES.Time(2:end),SmoothedSpeed)
% xlabel('time'); ylabel('Speed');
% if exist('SpeedThreshold','var')
%     RunningThreshold = 1;
%     PercRunning = length(SmoothedSpeed(find(SmoothedSpeed>RunningThreshold)))/length(SmoothedSpeed)*100;
%     RunningThreshold = SpeedThreshold;
% else
%     RunningThreshold = 1;
%     PercRunning = length(SmoothedSpeed(find(SmoothedSpeed>RunningThreshold)))/length(SmoothedSpeed)*100;
%     RunningThreshold = -10; % dummy value for the other plot: show all trials.
% end
% 
% 
% ylabel('speed (cm/s)'); set(gca,'TickDir','out');
% ylim([min(SmoothedSpeed) max(SmoothedSpeed)])
% title(['avg speed=' num2str(mean(SmoothedSpeed),3) 'cm/s, ' num2str(PercRunning,3) '% running'])
% box off
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if exist('PertStatus','var')
%     subplot(6,1,3:6) % plot speed of trials
%     hold off
%     runningtrials = 0; k = 1;
%     clear TempAngleSpeed runningtrials
%     if exist('Direction','var')
%         TempDirecton = Direction;
%     else
%         TempDirecton = sort(unique(ES.OrienTrial));
%     end
%     for Dir = TempDirecton
%         textON = 1;
%         counter = 1;
%         clear TempAngleSpeed
%         for i = 1:size(SpeedTrial,1)
%             if mean(SpeedTrial(ind_(i)))>RunningThreshold
%                 if ES.OrienTrial(ind_(i))==Dir && ES.PertTrial(ind_(i))==PertStatus
%                     hold on
%                     scatter(linspace(-trialSide_seconds,(ES.ACInfo.trialStartsEnds(j,2)+trialSide_seconds)-(ES.ACInfo.trialStartsEnds(j,1)-trialSide_seconds),size(SpeedTrial,2)),...
%                         k*ones(length(SpeedTrial(ind_(i),:)),1), ...
%                         [], SpeedTrial(ind_(i),:),'filled');
%                     if textON == 1
%                         text((ES.ACInfo.trialStartsEnds(j,2)+trialSide_seconds)-(ES.ACInfo.trialStartsEnds(j,1)-trialSide_seconds)+0.1,k,num2str(Dir))
%                         textON = 0;
%                     end
%                     TempAngleSpeed(counter,:)=SpeedTrial(ind_(i),:);
%                     counter = counter + 1;
%                     runningtrials(k) = i;
%                     k = k+1;                    
%                 end
%             end
%         end 
%     end
%     colorbar
%     xlim([-trialSide_seconds ES.ACInfo.trialStartsEnds(j,2)-ES.ACInfo.trialStartsEnds(j,1)+trialSide_seconds+2])
%     ylabel('trial'); set(gca,'TickDir','out');
%     xlabel('seconds');
%     if exist('PertStatus','var')
%         title(['Trials perturbation status=' num2str(PertStatus) ', (n=' num2str(k) '), running thres=' num2str(RunningThreshold) ', direction=' num2str(TempDirecton)])
%     elseif ~exist('PertStatus','var')
%         title(['All trials, (n=' num2str(k) '), running thres=' num2str(RunningThreshold) ', direction=' num2str(TempDirecton)])
%     end 
%         plot([ES.ACInfo.trialStartsEnds(j,2)-ES.ACInfo.trialStartsEnds(j,1) ES.ACInfo.trialStartsEnds(j,2)-ES.ACInfo.trialStartsEnds(j,1)],...
%             [min(min(SmoothedSpeed)) k],'r','LineWidth',2.5,'HandleVisibility','off');
%         hold on
%         plot([ES.ACInfo.trialStartsEnds(j,1)-ES.ACInfo.trialStartsEnds(j,1) ES.ACInfo.trialStartsEnds(j,1)-ES.ACInfo.trialStartsEnds(j,1)],...
%             [min(min(SmoothedSpeed)) k],'r','LineWidth',2.5,'HandleVisibility','off');
%         hold on
%         plot([ES.Time(ES.PertStartsEnds(end,2))-ES.Time(ES.trialStartsEnds(max(find(ES.PertTrial(:)==1)))) ES.Time(ES.PertStartsEnds(end,2))-ES.Time(ES.trialStartsEnds(max(find(ES.PertTrial(:)==1))))],...
%             [min(min(SmoothedSpeed)) k],'k','LineWidth',1.5,'HandleVisibility','off');
%         hold on
%         plot([ES.Time(ES.PertStartsEnds(end,1))-ES.Time(ES.trialStartsEnds(max(find(ES.PertTrial(:)==1)))) ES.Time(ES.PertStartsEnds(end,1))-ES.Time(ES.trialStartsEnds(max(find(ES.PertTrial(:)==1))))],...
%              [min(min(SmoothedSpeed)) k],'k','LineWidth',1.5,'HandleVisibility','off');
%         ylim([0 k]); xlim([-trialSide_seconds ES.ACInfo.trialStartsEnds(j,2)-ES.ACInfo.trialStartsEnds(j,1)+trialSide_seconds]);
%         box off
% else
%     subplot(6,1,3:6) % plot speed of trials
%     hold off
%     runningtrials = 0; k = 1;
%     clear TempAngleSpeed runningtrials
%     if exist('Direction','var')
%         TempDirecton = Direction;
%     else
%         TempDirecton = sort(unique(ES.OrienTrial));
%     end
%     for Dir = TempDirecton
%         textON = 1;
%         counter = 1;
%         clear TempAngleSpeed
%         for i = 1:size(SpeedTrial,1)
%             if mean(SpeedTrial(ind_(i)))>RunningThreshold
%                 if ES.OrienTrial(ind_(i))==Dir
%                     hold on
%                     scatter(linspace(-trialSide_seconds,(ES.ACInfo.trialStartsEnds(j,2)+trialSide_seconds)-(ES.ACInfo.trialStartsEnds(j,1)-trialSide_seconds),size(SpeedTrial,2)),...
%                         k*ones(length(SpeedTrial(ind_(i),:)),1), ...
%                         [], SpeedTrial(ind_(i),:),'filled');
%                     if textON == 1
%                         text((ES.ACInfo.trialStartsEnds(j,2)+trialSide_seconds)-(ES.ACInfo.trialStartsEnds(j,1)-trialSide_seconds)+0.1,k,num2str(Dir))
%                         textON = 0;
%                     end
%                     TempAngleSpeed(counter,:)=SpeedTrial(ind_(i),:);
%                     counter = counter + 1;
%                     runningtrials(k) = i;
%                     k = k+1;
%                 end
%             end
%         end
%     end
%     colorbar
%     xlim([-trialSide_seconds ES.ACInfo.trialStartsEnds(j,2)-ES.ACInfo.trialStartsEnds(j,1)+trialSide_seconds+2])
%     ylabel('trial'); set(gca,'TickDir','out');
%     xlabel('seconds');
%     if exist('PertStatus','var')
%         title(['Trials perturbation status=' num2str(PertStatus) ', (n=' num2str(length(runningtrials)) '), running thres=' num2str(RunningThreshold) ', direction=' num2str(TempDirecton)])
%     elseif ~exist('PertStatus','var')
%         title(['All trials, (n=' num2str(k) '), running thres=' num2str(RunningThreshold) ', direction=' num2str(TempDirecton)])
%     else
%         
%         plot([ES.ACInfo.trialStartsEnds(j,2)-ES.ACInfo.trialStartsEnds(j,1) ES.ACInfo.trialStartsEnds(j,2)-ES.ACInfo.trialStartsEnds(j,1)],...
%             [min(min(SmoothedSpeed)) k],'r','LineWidth',2.5,'HandleVisibility','off');
%         hold on
%         plot([ES.ACInfo.trialStartsEnds(j,1)-ES.ACInfo.trialStartsEnds(j,1) ES.ACInfo.trialStartsEnds(j,1)-ES.ACInfo.trialStartsEnds(j,1)],...
%             [min(min(SmoothedSpeed)) k],'r','LineWidth',2.5,'HandleVisibility','off');
%         ylim([0 k]); xlim([-trialSide_seconds ES.ACInfo.trialStartsEnds(j,2)-ES.ACInfo.trialStartsEnds(j,1)+trialSide_seconds]);
%         box off
%     end
% end
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pattern = 'Matlab';
% StartLetter = strfind(FileName{1},pattern);
% DetailCell = {[FileName{1}(StartLetter+length(pattern)+1:StartLetter+length(pattern)+13) ', ' FileName{1}(end-22:end-4)]};
% % a = axes;
% % t1 = title(DetailCell);
% % a.Visible = 'off'; % set(a,'Visible','off');
% % t1.Visible = 'on'; % set(t1,'Visible','on');
% 
% % % % CORRECT INFO FOR FIGURE TITLE AND MAKE IT UNIVERSAL
% % check folder if exists
% pattern = 'Matlab';
% StartLetter = strfind(FileName{1},pattern);
% DirInfo = {'AllDir', num2str(Direction)};
% [DirSave ind_dir] = min([length(DirInfo{1}),length(DirInfo{2})]);
% if exist([FileName{1}(1:StartLetter-1) 'Figures'], 'dir')
%     if exist('PertStatus','var')
%         saveas(gcf,[FileName{1}(1:StartLetter-1) 'Figures' filesep ...
%             'Speed_' num2str(RecordingOI) '_Direction_' DirInfo{ind_dir} '_RunTH_' num2str(RunningThreshold) '_Pert_' num2str(PertStatus)],'png')
%     else
%         saveas(gcf,[FileName{1}(1:StartLetter-1) 'Figures' filesep ...
%             'Speed_' num2str(RecordingOI) '_Direction_' DirInfo{ind_dir} '_RunTH_' num2str(RunningThreshold) '_'],'png')
%     end
% else
%     mkdir([FileName{1}(1:StartLetter-1) 'Figures'])
%     % save plot in the relevant folder
%     if exist('PertStatus','var')
%         saveas(gcf,[FileName{1}(1:StartLetter-1) 'Figures' filesep ...
%             'Speed_' num2str(RecordingOI) '_Direction_' DirInfo{ind_dir} '_RunTH_' num2str(RunningThreshold) '_Pert_' num2str(PertStatus)],'png')
%     else
%         saveas(gcf,[FileName{1}(1:StartLetter-1) 'Figures' filesep ...
%             'Speed_' num2str(RecordingOI) '_Direction_' DirInfo{ind_dir} '_RunTH_' num2str(RunningThreshold) '_'],'png')
%     end
% end
% close all

%EOF
end
