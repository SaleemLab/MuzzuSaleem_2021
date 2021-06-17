%% Tomaso Muzzu - UCL - 16/05/2019

% function to plot basic PSTH to visual responses looking only at a
% for all cells of all recordings


function AllCellResponse_runstill(ProjectData,DiffDRPer_total,value2show,cond)

% Gaussian filter for smoothing out firing rate
GaussFilter_W = 0.6; % seconds
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
    clear NormFR_run NormFR_still NormFR_diff
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
        cellcount = cellcount + 1;
        
        %% compute mean speed
        if i == 1
            c1 = 1; clear PertStart_i SpikeCount_FR SpikeCount_FR_init SpikeCount_FR_pert SpikeCount_FR_postpert
            for ss = 1: size(SpikeCount_smooth,1)
                if ProjectData.Trial_data{rr,1}.PerturbationON(ss,1)==1
                    PertStart_i(c1,1) = min(find(ProjectData.Trial_data{rr,1}.TF_in_trial{ss}==0)); % find index start of perturbation
                    PertStart_i(c1,2) = max(find(ProjectData.Trial_data{rr,1}.TF_in_trial{ss}==0)); % find  duration of perturbation
                    PertStart_i(c1,3) = ss; % perturbation trial ID
                    PertStart_i(c1,4) = mean(ProjectData.Trial_data{rr,1}.Speed_in_trial{PertStart_i(c1,3)}(PertStart_i(c1,1):PertStart_i(c1,2))); % mean speed
                    c1 = c1 + 1;
                end
            end
            % save mean speeds over the three periods of interest: initial 4
            % seconds, perturbation time window, 1s post perturbation
            OrigPertTrialOrder = ProjectData.Trial_data{rr,1}.PerturbationON;
            clear TrialSpeeds
            for ss = 1:size(SpikeCount_smooth,1)
                TrialSpeeds(ss,:) = [ mean(ProjectData.Trial_data{rr,1}.Speed_in_trial{ss,1}(round(trialSide_seconds*BonvisionFR):round(trialSide_seconds*BonvisionFR)+60*4)),... % initial 4 s
                    mean(ProjectData.Trial_data{rr,1}.Speed_in_trial{ss,1}(round(mean(PertStart_i(:,1))):round(mean(PertStart_i(:,2))))),...
                    mean(ProjectData.Trial_data{rr,1}.Speed_in_trial{ss,1}(round(mean(PertStart_i(:,2))):round(mean(PertStart_i(:,2))+60)))];
            end
        end
        %%
        clear FR
        FR(1,:) = sum(SpikeCount_smooth(ProjectData.Trial_data{rr,1}.PerturbationON==0 & TrialSpeeds(:,2)>3,:),1) ...
                    / sum(ProjectData.Trial_data{rr,1}.PerturbationON==0 & TrialSpeeds(:,2)>3);
        FR(2,:) = sum(SpikeCount_smooth(ProjectData.Trial_data{rr,1}.PerturbationON==1 & TrialSpeeds(:,2)>3,:),1) ...
                    / sum(ProjectData.Trial_data{rr,1}.PerturbationON==1 & TrialSpeeds(:,2)>3);
        FR(3,:) = sum(SpikeCount_smooth(ProjectData.Trial_data{rr,1}.PerturbationON==0 & TrialSpeeds(:,2)<3,:),1) ...
                    / sum(ProjectData.Trial_data{rr,1}.PerturbationON==0 & TrialSpeeds(:,2)<3);
        FR(4,:) = sum(SpikeCount_smooth(ProjectData.Trial_data{rr,1}.PerturbationON==1 & TrialSpeeds(:,2)<3,:),1) ...
                    / sum(ProjectData.Trial_data{rr,1}.PerturbationON==1 & TrialSpeeds(:,2)<3);
        FR(5,:) = abs(FR(2,:)-FR(4,:));
        
        if value2show == 0
            NormFR_run(i,:) = (FR(2,:)-FR(1,:))/max(max(FR(1:2,:)));
            NormFR_still(i,:) = (FR(4,:)-FR(3,:))/max(max(FR(1:2,:)));
        else
            NormFR_run(i,:) = FR(2,:)/max(FR(2,:));
            NormFR_still(i,:) = FR(4,:)/max(FR(4,:));
            NormFR_diff(i,:) = FR(5,:);%/max(FR(5,:));
        end
    end
    n_PT_cond(rr,1) = sum(ProjectData.Trial_data{rr,1}.PerturbationON==1 & TrialSpeeds(:,2)>3);
    n_PT_cond(rr,2) = sum(ProjectData.Trial_data{rr,1}.PerturbationON==1 & TrialSpeeds(:,2)<3);
    n_PT_cond(rr,3) = size(ProjectData.Units_Info{rr,1},1);
    if rr == 1
        NormFR_allcells_run = NormFR_run(:,1:559);
        NormFR_allcells_still = NormFR_still(:,1:559);
        NormFR_allcells_diff = NormFR_diff(:,1:559);
    else
        NormFR_allcells_run = [NormFR_allcells_run; NormFR_run(:,1:559)];
        NormFR_allcells_still = [NormFR_allcells_still; NormFR_still(:,1:559)];
        NormFR_allcells_diff = [NormFR_allcells_diff; NormFR_diff(:,1:559)];
    end
    
end


%% sessions of interest
SofI = find(n_PT_cond(:,1)>10 & n_PT_cond(:,2)>10);
UnitsGroups = n_PT_cond(SofI,3);
n_PT_cond(:,4) = cumsum(n_PT_cond(:,3));
for f = 1:length(SofI)
    if SofI(f) == 1
        Units_ofI = 1:n_PT_cond(SofI(f),4);
    elseif f == 1 & SofI(f) ~= 1
        Units_ofI = n_PT_cond(SofI(f)-1,4):n_PT_cond(SofI(f),4);
    else
        Units_ofI = [Units_ofI n_PT_cond(SofI(f)-1,4):n_PT_cond(SofI(f),4)];
    end
end

%% sort the cells by their response to perturbation 

if cond == 1
    SortingColumns = [PertRespMean(Units_ofI,1), ...    % perturbation period
                      PertRespMean(Units_ofI,2), ...    % second after perturbation
                      mean(NormFR_allcells_run(Units_ofI,trialSide_seconds:trialSide_seconds+60*2),2)]; % first 2 seconds of trial
elseif cond == 0
    SortingColumns = [PertRespMean(Units_ofI,1), ...    % perturbation period
                      PertRespMean(Units_ofI,2), ...    % second after perturbation
                      mean(NormFR_allcells_still(Units_ofI,trialSide_seconds:trialSide_seconds+60*2),2)]; % first 2 seconds of trial
elseif cond == 2 
    %NormFR_allcells_diff = NormFR_allcells_run-NormFR_allcells_still;
    SortingColumns = [PertRespMean(Units_ofI,1), ...    % perturbation period
                      PertRespMean(Units_ofI,2), ...    % second after perturbation
                      mean(NormFR_allcells_diff(Units_ofI,trialSide_seconds:trialSide_seconds+60*2),2)]; % first 2 seconds of trial
    SortingColumns =  mean(NormFR_allcells_diff(Units_ofI,trialSide_seconds+60*4.5:trialSide_seconds+60*5.5),2); % trial for pert periods of selected cells           
                            
end
clear B I
[B,I] = sortrows(SortingColumns,[-1]);

%% make actual figure
TimeLine = linspace(ProjectData.Trial_data{1,1}.time_in_trial{1,1}(1)-trialSide_seconds,ProjectData.Trial_data{1,1}.time_in_trial{1,1}(end-1)-trialSide_seconds,size(NormFR_allcells_run,2));
figure('Renderer', 'painters', 'Position', [10 10 1000 1000])
subplot(6,1,[2:6])
colormap(parula(256))
if cond == 1
    Responses2plot = NormFR_allcells_run(Units_ofI,:,:);
    imagesc(TimeLine,linspace(1,size(Responses2plot,1),size(Responses2plot,1)), ...
        Responses2plot(I,:));
elseif cond == 0
    Responses2plot = NormFR_allcells_still(Units_ofI,:,:);
    imagesc(TimeLine,linspace(1,size(Responses2plot,1),size(Responses2plot,1)), ...
        Responses2plot(I,:));
elseif cond == 2
    Responses2plot = NormFR_allcells_diff(Units_ofI,:,:);
    imagesc(TimeLine,linspace(1,size(Responses2plot,1),size(Responses2plot,1)), ...
        Responses2plot(I,:)./max(Responses2plot(I,:),[],2));
end
colorbar
set(gca,'FontSize',15);
xlim([-0.5 max(TimeLine)-0.5]);
xlabel('seconds','fontsize',15); ylabel('units','fontsize',15);
TotNrPertTrials = 0;
for rec = 1:size(DiffDRPer_total,1)
    TotNrPertTrials = TotNrPertTrials + DiffDRPer_total.TrialTypes{rec,1}(1);
end

% plot also contrast and TF
subplot(6,1,1)
hold off
PerturbationTrialExample = min(find(ProjectData.Trial_data{1,1}.PerturbationON));
h1 = plot(TimeLine,ProjectData.Trial_data{1,1}.contrast_in_trial{PerturbationTrialExample,1}(1:end-1),'k','LineWidth',2)
hold on
h2 = plot(TimeLine,ProjectData.Trial_data{1,1}.TF_in_trial{PerturbationTrialExample,1}(1:end-1),'b','LineWidth',2)
hold on
plot([0 0],...
    [0 max(ProjectData.Trial_data{1,1}.TF_in_trial{1,1}(1:end-1))],'--r','LineWidth',0.75);
hold on
plot([ProjectData.Session_data{1,1}.ACInfo{1,1}.trialStartsEnds(1,2)-ProjectData.Session_data{1,1}.ACInfo{1,1}.trialStartsEnds(1,1)...
      ProjectData.Session_data{1,1}.ACInfo{1,1}.trialStartsEnds(1,2)-ProjectData.Session_data{1,1}.ACInfo{1,1}.trialStartsEnds(1,1)],...
    [0 max(ProjectData.Trial_data{1,1}.TF_in_trial{1,1}(1:end-1))],'--r','LineWidth',0.75);
hold on
plot([TimeLine(min(find(ProjectData.Trial_data{1,1}.TF_in_trial{PerturbationTrialExample,1}(1:end-1)==0))) ...
      TimeLine(max(find(ProjectData.Trial_data{1,1}.TF_in_trial{PerturbationTrialExample,1}(1:end-1)==0)))],...
    [max(ProjectData.Trial_data{1,1}.TF_in_trial{PerturbationTrialExample,1}(1:end-1)) max(ProjectData.Trial_data{1,1}.TF_in_trial{PerturbationTrialExample,1}(1:end-1))],'--b','LineWidth',1);

xlim([-0.5 max(TimeLine)-0.5]);
ylim([0 max(ProjectData.Trial_data{1,1}.TF_in_trial{PerturbationTrialExample,1}(1:end-1))+0.4]);
box off
set(gca,'TickDir','out','XTickLabel',[],'YTick',[0.8 3],'YTickLabel',[0.8 3])    
legend({'Contrast', 'TF'},'location','northwest')
set(gca,'FontSize',15);
colorbar
if value2show == 0
    title(['Diff response of pert. trials - no-pert. trials '])
else
    if cond == 0
        title(['Response to perturbation trials during STILL '])
    elseif cond ==1
        title(['Response to perturbation trials during RUN '])
    elseif cond ==2
        title(['Response to perturbation trials RUN-STILL '])
    end
end
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




  





