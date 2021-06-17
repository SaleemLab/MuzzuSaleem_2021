%% Tomaso Muzzu - UCL - 06/02/2019

% function to plot basic PSTH to visual responses looking only at a
% specific angle for all cells


function SingleCell_response_PolarPlots(Array,SpikesTrial,ES,FileName,RecordingOI,Direction, RunningThreshold, MaxSpeed,StimulusParams,value2show)

trialSide_seconds = 1; % take 60 samples before and after trial
trialSide_samples = trialSide_seconds*round(length(ES.Time)/(max(ES.Time)-min(ES.Time))); % take 60 samples before and after trial

clear Tuning_dir Tuning_pert Tuning_postpert Tuning_dir_s Tuning_pert_s Tuning_postpert_s
for i = 1:size(SpikesTrial,1)
    clear SingletrialSpikes SpikeCount SpeedTrial
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
    
    %FR(1,:) = sum(SpikeCount(Array.PerturbationON==0,:)) / sum(Array.PerturbationON==0);
    %FR(2,:) = sum(SpikeCount(Array.PerturbationON==1,:)) / sum(Array.PerturbationON==1);
    clear FR_dir FR_dir_n FR_pert FR_pert_n FR_postpert FR_postpert_n
    % responses to different directions
    TrialInset = trialSide_samples:min(find(StimulusParams(:,2)==max(StimulusParams(:,2))))-trialSide_samples; % select first 4 seconds of trial
    for k = 1:length(Direction)
        FR_dir(k,:) = sum(SpikeCount(Array.Direction==Direction(k),TrialInset)) / sum(Array.Direction==Direction(k));
        FR_dir(k,:) = smooth(FR_dir(k,:),5);
        FR_dir_n(k) = sum(Array.Direction==Direction(k));
    end    
    % responses to different perturbation directions
    PertInset = min(find(StimulusParams(:,3)==min(StimulusParams(:,3)))):max(find(StimulusParams(:,3)==min(StimulusParams(:,3)))); % select first 4 seconds of trial
    for k = 1:length(Direction)
        FR_pert(k,:) = sum(SpikeCount(Array.Direction==Direction(k) & Array.PerturbationON==1,PertInset)) / sum(Array.Direction==Direction(k) & Array.PerturbationON==1);
        FR_pert(k,:) = smooth(FR_pert(k,:),5);
        FR_pert_n(k) = sum(Array.Direction==Direction(k)& Array.PerturbationON==1);
        FR_pert_ctrl(k,:) = sum(SpikeCount(Array.Direction==Direction(k) & Array.PerturbationON==0,PertInset)) / sum(Array.Direction==Direction(k) & Array.PerturbationON==0);
        FR_pert_ctrl(k,:) = smooth(FR_pert_ctrl(k,:),5);
        FR_pert_n_ctrl(k) = sum(Array.Direction==Direction(k)& Array.PerturbationON==0);
    end
    % responses to post perturbation periods
    PostPertInset = max(find(StimulusParams(:,3)==min(StimulusParams(:,3)))):length(StimulusParams(:,3))-trialSide_samples; % select first 4 seconds of trial
    for k = 1:length(Direction)
        FR_postpert(k,:) = sum(SpikeCount(Array.Direction==Direction(k) & Array.PerturbationON==1,PostPertInset)) / sum(Array.Direction==Direction(k) & Array.PerturbationON==1);
        FR_postpert(k,:) = smooth(FR_postpert(k,:),5);
        FR_postpert_n(k) = sum(Array.Direction==Direction(k) & Array.PerturbationON==1);
        FR_postpert_ctrl(k,:) = sum(SpikeCount(Array.Direction==Direction(k) & Array.PerturbationON==0,PostPertInset)) / sum(Array.Direction==Direction(k) & Array.PerturbationON==0);
        FR_postpert_ctrl(k,:) = smooth(FR_postpert_ctrl(k,:),5);
        FR_postpert_n_ctrl(k) = sum(Array.Direction==Direction(k) & Array.PerturbationON==0);
    end
    if isfield(ES,'Trials')
        FR_baseline(i) = mean(smooth(SpikeCount(:,1:trialSide_samples),5));
    else
        FR_baseline(i) = mean(smooth(SpikeCount(:,1:trialSide_samples),5));
    end
    
    Tuning_dir(i,:) = mean(FR_dir,2)' ;
    Tuning_dir_s(i,:) = std(FR_dir')./sqrt(FR_dir_n);
    Tuning_pert(i,:) = mean(FR_pert,2)' ;
    Tuning_pert_s(i,:) = std(FR_pert')./sqrt(FR_pert_n);
    Tuning_pert_ctrl(i,:) = mean(FR_pert_ctrl,2)' ;
    Tuning_pert_s_ctrl(i,:) = std(FR_pert_ctrl')./sqrt(FR_pert_n_ctrl);
    Tuning_postpert(i,:) = mean(FR_postpert,2)' ;
    Tuning_postpert_s(i,:) = std(FR_postpert')./sqrt(FR_postpert_n);
    Tuning_postpert_ctrl(i,:) = mean(FR_postpert_ctrl,2)' ;
    Tuning_postpert_s_ctrl(i,:) = std(FR_postpert_ctrl')./sqrt(FR_postpert_n_ctrl);
    
%     if value2show == 0
%         NormFR(k,:) = (FR(2,:)-FR(1,:))/max(max(FR));
%     else
%         NormFR(k,:) = (FR(2,:))/(max(FR(2,:)));
%     end
end

%% plot the direction tuning for each cell
clear SingletrialSpikes SpikeCount SpeedTrial
for i = 1:size(SpikesTrial,1) 
    if mod(i,30)==1
        figure('Renderer', 'painters', 'Position', [10 10 2000 1000])
    end
    if mod(i,30)==0
        subplot(5,6,30)
    else
        subplot(5,6,mod(i,30))
    end
    h1 = errorbar(Direction,Tuning_dir(i,:)*60,Tuning_dir_s(i,:)*60,'k');
    hold on
%     h11 = errorbar(Direction,(Tuning_dir(i,:)*60)-FR_baseline(i)*60,Tuning_dir_s(i,:)*60,'k:');
%     hold on
    h2 = errorbar(Direction,(Tuning_pert(i,:)-Tuning_pert_ctrl(i,:))*60,Tuning_pert_s(i,:)*60,'r');
    hold on
    h3 = errorbar(Direction,(Tuning_postpert(i,:)-Tuning_postpert_ctrl(i,:))*60,Tuning_postpert_s(i,:)*60','b');
    hold on
    h4 = plot(Direction, FR_baseline(i)*ones(1,length(Direction))*60,'Color',[0.5 0.5 0.5])
    set(gca,'XTick',Direction);
    xlim([min(Direction) max(Direction)])
    box off
    if i == 1
        legend({'Direction tuning' , ...
                '\Delta Perturbation response',...
                '\Delta Post-perturbation response',...
                'baseline'},'Location','northwest')
    end
    if i<25
        set(gca,'XTickLabel',[])
    end
    if mod(i,6)==1 || i==1
        ylabel('Hz');
    end
    if mod(i,30)==0 || i==size(SpikesTrial,1)
        pattern = 'Matlab';
        StartLetter = strfind(FileName,pattern);
        nrTrialsPertDir = length(find(ES.PertTrial==1));
        nrTrialsNPertDir = length(find(ES.PertTrial==0));
        DetailCell = {[FileName(StartLetter+length(pattern)+1:StartLetter+length(pattern)+13) ', ' FileName(end-22:end-4)] , ...
            char(strcat(['Directional responses for all directions, n_p=' num2str(nrTrialsPertDir) ', n_n_p=' num2str(nrTrialsNPertDir)])), ...
            char(strcat(['All Directions (' num2str(Direction') ', n_d=' num2str(FR_dir_n) ', n_p/n_p_p=' num2str(FR_pert_n)]))};
        DetailCell
        a = axes;
        t1 = title(DetailCell);
        a.Visible = 'off'; % set(a,'Visible','off');
        t1.Visible = 'on'; % set(t1,'Visible','on');
    end
end

% %% plot the perturbation direction tuning for each cell
% clear SingletrialSpikes SpikeCount SpeedTrial
% for i = 1:size(SpikesTrial,1) 
%     if mod(i,30)==1
%         figure('Renderer', 'painters', 'Position', [10 10 2000 1000])
%     end
%     if mod(i,30)==0
%         subplot(5,6,30)
%     else
%         subplot(5,6,mod(i,30))
%     end
%     errorbar(Direction,(Tuning_pert(i,:)-Tuning_pert_ctrl(i,:))*60,Tuning_pert_s(i,:)*60);
%     hold on
%     plot(Direction, FR_baseline(i)*ones(1,length(Direction))*60,'Color',[0.5 0.5 0.5])
%     set(gca,'XTick',Direction);
%     xlim([min(Direction) max(Direction)])
%     box off
%     if i == 1
%         legend({['no pert (n=' num2str(sum(Array.PerturbationON==0)) ')'],['with pert. (n=' num2str(sum(Array.PerturbationON==1)) ')']},'Location','northwest')
%     end
%     if i<25
%         set(gca,'XTickLabel',[])
%     end
%     if mod(i,6)==1 || i==1
%         ylabel('Hz');
%     end
%     if mod(i,30)==0 || i==size(SpikesTrial,1)
%         pattern = 'Matlab';
%         StartLetter = strfind(FileName,pattern);
%         nrTrialsPertDir = length(find(ES.PertTrial==1));
%         nrTrialsNPertDir = length(find(ES.PertTrial==0));
%         DetailCell = {[FileName(StartLetter+length(pattern)+1:StartLetter+length(pattern)+13) ', ' FileName(end-22:end-4)] , ...
%             char(strcat(['Pert. Resp - no pert. response, n_p=' num2str(nrTrialsPertDir) ', n_n_p=' num2str(nrTrialsNPertDir)]))};
%         DetailCell
%         a = axes;
%         t1 = title(DetailCell);
%         a.Visible = 'off'; % set(a,'Visible','off');
%         t1.Visible = 'on'; % set(t1,'Visible','on');
%     end
% end
% 
% %% plot the post perturbation direction tuning for each cell
% clear SingletrialSpikes SpikeCount SpeedTrial
% for i = 1:size(SpikesTrial,1) 
%     if mod(i,30)==1
%         figure('Renderer', 'painters', 'Position', [10 10 2000 1000])
%     end
%     if mod(i,30)==0
%         subplot(5,6,30)
%     else
%         subplot(5,6,mod(i,30))
%     end
%     errorbar(Direction,(Tuning_postpert(i,:)-Tuning_postpert_ctrl(i,:))*60,Tuning_postpert_s(i,:)*60);
%     hold on
%     plot(Direction, FR_baseline(i)*ones(1,length(Direction))*60,'Color',[0.5 0.5 0.5])
%     set(gca,'XTick',Direction);
%     xlim([min(Direction) max(Direction)])
%     box off
%     if i == 1
%         legend({['no pert (n=' num2str(sum(Array.PerturbationON==0)) ')'],['with pert. (n=' num2str(sum(Array.PerturbationON==1)) ')']},'Location','northwest')
%     end
%     if i<25
%         set(gca,'XTickLabel',[])
%     end
%     if mod(i,6)==1 || i==1
%         ylabel('Hz');
%     end
%     if mod(i,30)==0 || i==size(SpikesTrial,1)
%         pattern = 'Matlab';
%         StartLetter = strfind(FileName,pattern);
%         nrTrialsPertDir = length(find(ES.PertTrial==1));
%         nrTrialsNPertDir = length(find(ES.PertTrial==0));
%         DetailCell = {[FileName(StartLetter+length(pattern)+1:StartLetter+length(pattern)+13) ', ' FileName(end-22:end-4)] , ...
%             char(strcat(['PostPert. Resp - no pert. response, n_p=' num2str(nrTrialsPertDir) ', n_n_p=' num2str(nrTrialsNPertDir)]))};
%         DetailCell
%         a = axes;
%         t1 = title(DetailCell);
%         a.Visible = 'off'; % set(a,'Visible','off');
%         t1.Visible = 'on'; % set(t1,'Visible','on');
%     end
% end


end  





