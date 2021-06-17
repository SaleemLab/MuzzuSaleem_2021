%% Tomaso Muzzu - UCL - 16/05/2019

% function to plot basic PSTH to visual responses looking only at a
% for all cells of all recordings


function AllCellsResponse_Direction(ProjectData,DiffDRPer_total,value2show,DirectionOI)

% Gaussian filter for smoothing out firing rate
GaussFilter_W = 0.3; % seconds
BonvisionFR = 60;
Width = round(GaussFilter_W*BonvisionFR);
Sigma = Width/3; % standard deviation in number of samples (converted from time in seconds)
x_g = linspace(-Width/2, Width/2, Width);
gaussFilter = exp(-x_g.^2/(2*Sigma^2));
gaussFilter_ = gaussFilter / sum (gaussFilter); % normalize

trialSide_seconds = 1; % take 60 samples before and after trial
trialSide_samples = 60; % take 60 samples before and after trial
cellcount = 1; 
clear PertRespMean
for rr = 1:size(ProjectData,1) 
    clear NormFR NormFR_sem NormFR_off NormFR_sem_off
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
            if ProjectData.Trial_data{rr,1}.PerturbationON(j,1) && ProjectData.Trial_data{rr,1}.Direction(j,1)==DirectionOI
                Pert_idx_start = min(find(ProjectData.Trial_data{rr,1}.TF_in_trial{j,:}==0));
                Pert_idx_end = max(find(ProjectData.Trial_data{rr,1}.TF_in_trial{j,:}==0));
                PertRespMeanCell(trialcount,1) = mean(SpikeCount_smooth(j,Pert_idx_start:Pert_idx_end)); % mean response during pert. for ranking all response
                PertRespMeanCell(trialcount,2) = mean(SpikeCount_smooth(j,Pert_idx_end+1:Pert_idx_end+1+60)); % mean response post perturbation
                trialcount = trialcount + 1;
            end
        end 
        PertRespMean(cellcount,1) = mean(PertRespMeanCell(:,1)); % mean response during pert. for ranking all responses
        PertRespMean(cellcount,2) = mean(PertRespMeanCell(:,2)); % mean response post perturbation
        cellcount = cellcount + 1;
        
        clear FR
        PertON_idx = find(ProjectData.Trial_data{rr,1}.PerturbationON==1 & ProjectData.Trial_data{rr,1}.Direction==DirectionOI);
        PertOFF_idx = find(ProjectData.Trial_data{rr,1}.PerturbationON==0 & ProjectData.Trial_data{rr,1}.Direction==DirectionOI);
        r_idx = randperm(length(PertOFF_idx));
        FR(1,:) = sum(SpikeCount_smooth(PertOFF_idx(r_idx(1:length(PertON_idx))),:),1) ...
            / sum(PertOFF_idx(r_idx(1:length(PertON_idx))));
        FR(2,:) = sum(SpikeCount_smooth(PertON_idx,:),1) ...
            / sum(PertON_idx);
        FR(3,:) = std(SpikeCount_smooth(PertOFF_idx(r_idx(1:length(PertON_idx))),:),1) ...
            / sqrt(sum(PertOFF_idx(r_idx(1:length(PertON_idx)))));
        FR(4,:) = std(SpikeCount_smooth(PertON_idx,:),1) ...
            / sqrt(sum(PertON_idx));
        
        if value2show == 0
            NormFR(i,:) = (FR(2,:)-FR(1,:))/max(max(FR(1:2,:)));
            NormFR_sem(i,:) = ((FR(4,:)+FR(3,:))/2)/max(max(FR(3:4,:)));
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
            
% sort the cells by their response to perturbation 
SortingColumns = [PertRespMean(:,1), ... % perturbation period
                  mean(NormFR_allcells(:,trialSide_seconds:trialSide_seconds+60*2),2),... % first 2 seconds of trial
                  PertRespMean(:,2)];   % second after perturbation
[B,I] = sortrows(SortingColumns,[-1 -2 -3]);

%% make actual figure
TimeLine = linspace(ProjectData.Trial_data{1,1}.time_in_trial{1,1}(1)-trialSide_seconds,ProjectData.Trial_data{1,1}.time_in_trial{1,1}(end)-trialSide_seconds,size(NormFR_allcells,2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Renderer', 'painters', 'Position', [10 10 1000 1000])
subplot(10,1,[5:10])
colormap(parula(256))
imagesc(TimeLine,linspace(1,size(NormFR_allcells,1),size(NormFR_allcells,1)), ...
    NormFR_allcells(I,:));
colorbar
set(gca,'FontSize',15);
xlabel('seconds','fontsize',15); ylabel('units','fontsize',15);
TotNrPertTrials = 0;
xlim([-0.3 max(TimeLine)-0.3]);

for rec = 1:size(DiffDRPer_total,1)
    TotNrPertTrials = TotNrPertTrials + DiffDRPer_total.TrialTypes{rec,1}(1);
end

% plot also contrast and TF
subplot(10,1,[1 2])
h1 = plot(TimeLine,ProjectData.Trial_data{1,1}.contrast_in_trial{1,1}(1:end-1),'k','LineWidth',2)
hold on
h2 = plot(TimeLine,ProjectData.Trial_data{1,1}.TF_in_trial{1,1}(1:end-1),'b','LineWidth',2)
hold on
plot([0 0],...
    [0 max(ProjectData.Trial_data{1,1}.TF_in_trial{1,1}(1:end-1))],'--r','LineWidth',0.75);
hold on
plot([ProjectData.Session_data{1,1}.ACInfo{1,1}.trialStartsEnds(1,2)-ProjectData.Session_data{1,1}.ACInfo{1,1}.trialStartsEnds(1,1)...
      ProjectData.Session_data{1,1}.ACInfo{1,1}.trialStartsEnds(1,2)-ProjectData.Session_data{1,1}.ACInfo{1,1}.trialStartsEnds(1,1)],...
    [0 max(ProjectData.Trial_data{1,1}.TF_in_trial{1,1}(1:end-1))],'--r','LineWidth',0.75);
hold on
plot([TimeLine(min(find(ProjectData.Trial_data{1,1}.TF_in_trial{1,1}(1:end-1)==0))) ...
      TimeLine(max(find(ProjectData.Trial_data{1,1}.TF_in_trial{1,1}(1:end-1)==0)))],...
    [max(ProjectData.Trial_data{1,1}.TF_in_trial{1,1}(1:end-1)) max(ProjectData.Trial_data{1,1}.TF_in_trial{1,1}(1:end-1))],'--b','LineWidth',1);
xlim([-0.3 max(TimeLine)-0.3]);
ylim([0 max(ProjectData.Trial_data{1,1}.TF_in_trial{1,1}(1:end-1))+0.4]);
box off
set(gca,'TickDir','out','XTickLabel',[],'YTick',[0.8 3],'YTickLabel',[0.8 3])    
legend({'Contrast', 'TF'},'location','northwest')
set(gca,'FontSize',15);
colorbar
if value2show == 0
    title(['Diff response of pert. trials - no-pert. trials, ' num2str(DirectionOI) '^{o}'])
else
    title(['Response during perturbation trials, ' num2str(DirectionOI) '^{o}'])
end

subplot(10,1,[3 4])
colormap(parula(256))
h1= shadedErrorBar(TimeLine,...
                    mean(NormFR_allcells(~isnan(sum(NormFR_allcells,2)),:),1),...
                    mean(NormFR_allcells_sem(~isnan(sum(NormFR_allcells,2)),:) ,1),'lineprops',{'r-','markerfacecolor','r'});
                %std(NormFR_allcells(~isnan(sum(NormFR_allcells,2)),:),1)/sqrt(length(NormFR_allcells(~isnan(sum(NormFR_allcells,2)),:))))
hold on
h2 = shadedErrorBar(TimeLine,...
                    mean(NormFR_allcells_off(~isnan(sum(NormFR_allcells_off,2)),:),1),...
                    mean(NormFR_allcells_off_sem(~isnan(sum(NormFR_allcells_off,2)),:) ,1))
colorbar
set(gca,'FontSize',15);
ylabel('units','fontsize',15);
box off
set(gca,'TickDir','out','XTickLabel',[]);
TotNrPertTrials = 0;
xlim([-0.3 max(TimeLine)-0.3]);
ylim([0.2 0.4])

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




  





