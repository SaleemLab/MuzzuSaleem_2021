%% Tomaso Muzzu - UCL - 16/05/2019

% function to plot basic PSTH to visual responses looking only at a
% for all cells of all recordings


function AllCellsResponse(ProjectData,DiffDRPer_total,value2show)

% Gaussian filter for smoothing out firing rate
GaussFilter_W = 0.5; % seconds
BonvisionFR = 60;
Width = round(GaussFilter_W*BonvisionFR);
Sigma = Width/3; % standard deviation in number of samples (converted from time in seconds)
x_g = linspace(-Width/2, Width/2, Width);
gaussFilter = exp(-x_g.^2/(2*Sigma^2));
gaussFilter_ = gaussFilter / sum (gaussFilter); % normalize

trialSide_seconds = 1; % take 60 samples before and after trial
trialSide_samples = 60; % take 60 samples before and after trial
cellcount = 1; 
clear PertRespMean UnitSpikesDurSession MeanUnitFR
for rr = 1:size(ProjectData,1) 
    clear NormFR NormFR_sem NormFR_off NormFR_sem_off UnitSpikesDurSession
    for i = 1:size(ProjectData.Units_Info{rr,1},1)
        clear SingletrialSpikes SpikeCount_temp SpikeCount_smooth edges PertRespMeanCell
        trialcount = 1;
        for j = 1:size(ProjectData.Trial_data{rr,1},1)
            SingletrialSpikes = ProjectData.Trial_data{rr,1}.SpikesTrial{j,1}{1,i};
            TrialStart = ProjectData.Session_data{rr,1}.ACInfo{1,1}.trialStartsEnds(j,1);
            TrialEnd = ProjectData.Session_data{rr,1}.ACInfo{1,1}.trialStartsEnds(j,2);
            edges = linspace(TrialStart-trialSide_seconds,TrialEnd+trialSide_seconds,((TrialEnd+trialSide_seconds)-(TrialStart-trialSide_seconds))*60);
            [SpikeCount_temp, Edges] = histcounts(SingletrialSpikes,edges);
            SpikeCount_smooth(j,1:length(SpikeCount_temp)) = conv(SpikeCount_temp, gaussFilter_, 'same');
            if ProjectData.Trial_data{rr,1}.PerturbationON(j,1)
                Pert_idx_start = min(find(ProjectData.Trial_data{rr,1}.TF_in_trial{j,:}==0));
                Pert_idx_end = max(find(ProjectData.Trial_data{rr,1}.TF_in_trial{j,:}==0));
                PertRespMeanCell(trialcount,1) = mean(SpikeCount_smooth(j,Pert_idx_start:Pert_idx_end)); % mean response during pert. for ranking all response
                PertRespMeanCell(trialcount,2) = mean(SpikeCount_smooth(j,Pert_idx_end+1:Pert_idx_end+1+60)); % mean response post perturbation
                trialcount = trialcount + 1;
            end
        end 
        PertRespMean(cellcount,1) = mean(PertRespMeanCell(:,1)); % mean response during pert. for ranking all response
        PertRespMean(cellcount,2) = mean(PertRespMeanCell(:,2)); % mean response post perturbation
        %cellcount = cellcount + 1;
        
%         clear FR
%         FR(1,:) = sum(SpikeCount_smooth(ProjectData.Trial_data{rr,1}.PerturbationON==0,:)) / sum(ProjectData.Trial_data{rr,1}.PerturbationON==0);
%         FR(2,:) = sum(SpikeCount_smooth(ProjectData.Trial_data{rr,1}.PerturbationON==1,:)) / sum(ProjectData.Trial_data{rr,1}.PerturbationON==1);
        
        clear FR
        PertON_idx = find(ProjectData.Trial_data{rr,1}.PerturbationON==1);
        PertOFF_idx = find(ProjectData.Trial_data{rr,1}.PerturbationON==0);
        r_idx = randperm(length(PertOFF_idx));
        % take an equal nr of trials for both the pertON and pertOFF
        % conditions
        FR(1,:) = sum(SpikeCount_smooth(PertOFF_idx(r_idx(1:length(PertON_idx))),:),1) / length(PertOFF_idx(r_idx(1:length(PertON_idx)))); % mean FR when pert OFF
        FR(2,:) = sum(SpikeCount_smooth(PertON_idx,:),1) / length(PertON_idx); % mean FR when pert ON
        FR(3,:) = std(SpikeCount_smooth(PertOFF_idx(r_idx(1:length(PertON_idx))),:),1) / sqrt(length(PertOFF_idx(r_idx(1:length(PertON_idx))))); % SEM FR when pert OFF
        FR(4,:) = std(SpikeCount_smooth(PertON_idx,:),1) / sqrt(length(PertON_idx)); % SEM FR when pert ON
        
        % find mean FR during session of each unit
        SessionStart = ProjectData.Session_data{rr,1}.ACInfo{1,1}.trialStartsEnds(1,1);
        SessionEnd = ProjectData.Session_data{rr,1}.ACInfo{1,1}.trialStartsEnds(end,2);
        
        UnitSpikesDurSession{cellcount,1} = ProjectData.Units_Info{rr,1}.Spiketimes{i,1}(...
                    ProjectData.Units_Info{rr,1}.Spiketimes{i,1}>=SessionStart & ...
                    ProjectData.Units_Info{rr,1}.Spiketimes{i,1}<=SessionEnd);
        MeanUnitFR(cellcount,1) = length(UnitSpikesDurSession{cellcount,1})/(SessionEnd-SessionStart); % mean spiking rate
        % Reliability measure
        Time_1_3 = (SessionEnd-SessionStart)/3;
        UnitSpikesDurSession{cellcount,2} = ProjectData.Units_Info{rr,1}.Spiketimes{i,1}(...
                    ProjectData.Units_Info{rr,1}.Spiketimes{i,1}>=SessionStart & ...
                    ProjectData.Units_Info{rr,1}.Spiketimes{i,1}<=SessionStart+Time_1_3);
        MeanUnitFR(cellcount,2) = length(UnitSpikesDurSession{cellcount,2})/((SessionStart+Time_1_3)-SessionStart); % mean spiking rate
        UnitSpikesDurSession{cellcount,3} = ProjectData.Units_Info{rr,1}.Spiketimes{i,1}(...
                    ProjectData.Units_Info{rr,1}.Spiketimes{i,1}>=SessionEnd-Time_1_3 & ...
                    ProjectData.Units_Info{rr,1}.Spiketimes{i,1}<=SessionEnd);
        MeanUnitFR(cellcount,3) = length(UnitSpikesDurSession{cellcount,3})/(SessionEnd-(SessionEnd-Time_1_3)); % mean spiking rate
        
%         MeanUnitFR(cellcount,2) = 
        % find mean FR during gray screen trials or during inter-trial
        % intervals
        clear Gray_SC Gray_dur
        if isfield(ProjectData.Session_data{rr,1}.ACInfo{1,1},'trialZeroContrast')
            if ~isempty(ProjectData.Session_data{rr,1}.ACInfo{1,1}.trialZeroContrast)
                for g = 1:size(ProjectData.Session_data{rr,1}.ACInfo{1,1}.trialZeroContrast,1)
                    Gray_SC(g) = length(ProjectData.Units_Info{rr,1}.Spiketimes{i,1}( ...
                    ProjectData.Units_Info{rr,1}.Spiketimes{i,1}>ProjectData.Session_data{rr,1}.ACInfo{1,1}.trialZeroContrast(g,2) & ...
                    ProjectData.Units_Info{rr,1}.Spiketimes{i,1}<ProjectData.Session_data{rr,1}.ACInfo{1,1}.trialZeroContrast(g,1)));
                    Gray_dur(g) = ProjectData.Session_data{rr,1}.ACInfo{1,1}.trialZeroContrast(g,1)-ProjectData.Session_data{rr,1}.ACInfo{1,1}.trialZeroContrast(g,2);
                end
                Gray_FR(cellcount) = mean(Gray_SC./Gray_dur);
            end
        end
        
        cellcount = cellcount + 1;
        
        if value2show == 0
            NormFR(i,:) = (FR(2,:)-FR(1,:))/max(max(FR));
        else
            % pert ON
            NormFR(i,:) = FR(2,:)/max(max(FR(1:2,:)));
            NormFR_sem(i,:) = FR(4,:)/max(max(FR(1:2,:)));
            % pert OFF
            NormFR_off(i,:) = FR(1,:)/max(max(FR(1:2,:)));
            NormFR_sem_off(i,:) = FR(3,:)/max(max(FR(1:2,:)));
        end
    end
    if rr == 1
        % pert ON
        NormFR_allcells = NormFR(:,1:559);
        NormFR_allcells_sem = NormFR_sem(:,1:559);
        % pert OFF
        NormFR_allcells_off = NormFR_off(:,1:559);
        NormFR_allcells_off_sem = NormFR_sem_off(:,1:559);
    else
        % pert ON
        NormFR_allcells = [NormFR_allcells; NormFR(:,1:559)];
        NormFR_allcells_sem = [NormFR_allcells_sem; NormFR_sem(:,1:559)];
        % pert OFF
        NormFR_allcells_off = [NormFR_allcells_off; NormFR_off(:,1:559)];
        NormFR_allcells_off_sem = [NormFR_allcells_off_sem; NormFR_sem_off(:,1:559)];
    end
end

clear MeanUnitFR_th
FR_thres = 0.5;
MeanUnitFR_th=find(MeanUnitFR(:,1)>FR_thres);
NormFR_allcells_ = NormFR_allcells(MeanUnitFR_th,:);
NormFR_allcells_off_ = NormFR_allcells_off(MeanUnitFR_th,:);
PertRespMean_ = PertRespMean(MeanUnitFR_th,:);
% sort the cells by their response to perturbation 
SortingColumns = [PertRespMean_(:,1), ... % perturbation period
                  mean(NormFR_allcells_(:,trialSide_seconds:trialSide_seconds+60*2),2),... % first 2 seconds of trial
                  PertRespMean_(:,2)];   % second after perturbation
[B,I] = sortrows(SortingColumns,[-1 -2 -3]);

NormFR_allcells_ = NormFR_allcells_(I,:);
NormFR_allcells_off_ = NormFR_allcells_off_(I,:);

% make actual figure
TimeLine = linspace(ProjectData.Trial_data{1,1}.time_in_trial{1,1}(1)-trialSide_seconds,ProjectData.Trial_data{1,1}.time_in_trial{1,1}(end-1)-trialSide_seconds,size(NormFR_allcells,2));
figure
set(gcf, 'PaperUnits', 'centimeters ');
%set(gcf, 'Position', [0 0 20 30]);
set(gcf, 'Renderer', 'painters');
subplot(10,1,[5:10])
colormap(hot(256))
imagesc(TimeLine,linspace(1,size(NormFR_allcells_,1),size(NormFR_allcells_,1)), ...
    NormFR_allcells_);
colorbar
set(gca,'fontsize',13,'TickDir','out');
xlabel('seconds','fontsize',13); ylabel('units','fontsize',13);
TotNrPertTrials = 0;
xlim([-0.3 max(TimeLine)-0.3]);
box off

for rec = 1:size(DiffDRPer_total,1)
    TotNrPertTrials = TotNrPertTrials + DiffDRPer_total.TrialTypes{rec,1}(1);
end

% plot also contrast and TF
subplot(10,1,1)
PerturbationTrialExample = min(find(ProjectData.Trial_data{1,1}.PerturbationON));
h1 = plot(TimeLine,ProjectData.Trial_data{1,1}.contrast_in_trial{PerturbationTrialExample,1}(1:end-1),'k','LineWidth',2)
hold on
h2 = plot(TimeLine,ProjectData.Trial_data{1,1}.TF_in_trial{PerturbationTrialExample,1}(1:end-1),'b','LineWidth',2)
hold on
plot([0 0],...
    [0 max(ProjectData.Trial_data{1,1}.TF_in_trial{PerturbationTrialExample,1}(1:end-1))],'--r','LineWidth',0.75);
hold on
plot([ProjectData.Session_data{1,1}.ACInfo{1,1}.trialStartsEnds(1,2)-ProjectData.Session_data{1,1}.ACInfo{1,1}.trialStartsEnds(1,1)...
      ProjectData.Session_data{1,1}.ACInfo{1,1}.trialStartsEnds(1,2)-ProjectData.Session_data{1,1}.ACInfo{1,1}.trialStartsEnds(1,1)],...
    [0 max(ProjectData.Trial_data{1,1}.TF_in_trial{1,1}(1:end-1))],'--r','LineWidth',0.75);
hold on
plot([TimeLine(min(find(ProjectData.Trial_data{1,1}.TF_in_trial{PerturbationTrialExample,1}(1:end-1)==0))) ...
      TimeLine(max(find(ProjectData.Trial_data{1,1}.TF_in_trial{PerturbationTrialExample,1}(1:end-1)==0)))],...
    [max(ProjectData.Trial_data{1,1}.TF_in_trial{PerturbationTrialExample,1}(1:end-1)) max(ProjectData.Trial_data{1,1}.TF_in_trial{PerturbationTrialExample,1}(1:end-1))],'--b','LineWidth',1);

xlim([-0.3 max(TimeLine)-0.3]);
ylim([0 max(ProjectData.Trial_data{1,1}.TF_in_trial{PerturbationTrialExample,1}(1:end-1))+0.4]);
box off
set(gca,'TickDir','out','XTickLabel',[],'YTick',[0 0.8 3],'YTickLabel',[0 0.8 3])    
ll = legend({'Contrast', 'TF'},'location','northwest','fontsize',10)
ll.Color = 'none'; ll.EdgeColor = 'none'; 
set(gca,'FontSize',13);
colorbar
if value2show == 0
    title(['Diff response of pert. trials - no-pert. trials'])
else
    title(['Response during perturbation trials'])
end

subplot(10,1,[2:4])
colormap(jet(256))
h1= shadedErrorBar(TimeLine,...
                    mean(NormFR_allcells_(~isnan(sum(NormFR_allcells_,2)),:),1),...
                    std(NormFR_allcells_(~isnan(sum(NormFR_allcells_,2)),:),1)/sqrt(length(NormFR_allcells_(~isnan(sum(NormFR_allcells_,2)),:))), ...
                    'lineprops',{'r-','markerfacecolor','r'});
                % mean(NormFR_allcells_sem(~isnan(sum(NormFR_allcells,2)),:),1), ...
                %
hold on
h2 = shadedErrorBar(TimeLine,...
                    mean(NormFR_allcells_off_(~isnan(sum(NormFR_allcells_off_,2)),:),1),...
                    std(NormFR_allcells_off_(~isnan(sum(NormFR_allcells_off_,2)),:),1)/sqrt(length(NormFR_allcells_off_(~isnan(sum(NormFR_allcells_off_,2)),:))))
                    % mean(NormFR_allcells_off_sem(~isnan(sum(NormFR_allcells_off,2)),:) ,1))
hold on
plot([0 0],...
    [0 max(ProjectData.Trial_data{1,1}.TF_in_trial{1,1}(1:end-1))],'--r','LineWidth',0.75);
hold on
plot([ProjectData.Session_data{1,1}.ACInfo{1,1}.trialStartsEnds(1,2)-ProjectData.Session_data{1,1}.ACInfo{1,1}.trialStartsEnds(1,1)...
      ProjectData.Session_data{1,1}.ACInfo{1,1}.trialStartsEnds(1,2)-ProjectData.Session_data{1,1}.ACInfo{1,1}.trialStartsEnds(1,1)],...
    [0 max(ProjectData.Trial_data{1,1}.TF_in_trial{1,1}(1:end-1))],'--r','LineWidth',0.75);
hold on
plot([TimeLine(min(find(ProjectData.Trial_data{1,1}.TF_in_trial{PerturbationTrialExample,1}(1:end-1)==0))) ...
      TimeLine(min(find(ProjectData.Trial_data{1,1}.TF_in_trial{PerturbationTrialExample,1}(1:end-1)==0)))],...
    [0 1],'--b','LineWidth',1);
hold on
plot([TimeLine(max(find(ProjectData.Trial_data{1,1}.TF_in_trial{PerturbationTrialExample,1}(1:end-1)==0))) ...
      TimeLine(max(find(ProjectData.Trial_data{1,1}.TF_in_trial{PerturbationTrialExample,1}(1:end-1)==0)))],...
    [0 1],'--b','LineWidth',1);
colorbar
set(gca,'FontSize',13);
ylabel('norm. response','fontsize',13);
box off
set(gca,'TickDir','out','XTickLabel',[],'YTick',[0.4:0.1:0.7],'YTickLabel',[0.4:0.1:0.7]);
TotNrPertTrials = 0;
xlim([-0.3 max(TimeLine)-0.3]);
ylim([0.6 0.7])
ll = legend({'perturbation ON', 'perturbation OFF'},'fontsize',10);
ll.Color = 'none'; ll.EdgeColor = 'none'; 



%
% pattern = 'Matlab';
% StartLetter = strfind(FileName,pattern);
% nrTrialsPertDir = length(find(ES.PertTrial==1));
% nrTrialsNPertDir = length(find(ES.PertTrial==0));
% DetailCell = {[FileName(StartLetter+length(pattern)+1:StartLetter+length(pattern)+13) ', ' FileName(end-22:end-4)] , ...
%     char(strcat(['All Cells, All Directions, n_p=' num2str(nrTrialsPertDir) ', n_n_p=' num2str(nrTrialsNPertDir)]))};
% DetailCell;
% a = axes;
% t1 = title(DetailCell);
% a.Visible = 'off'; % set(a,'Visible','off');
% t1.Visible = 'on'; % set(t1,'Visible','on');

end




  





