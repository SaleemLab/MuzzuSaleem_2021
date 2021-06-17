%% Tomaso Muzzu - UCL - 12/04/2019

% statistical test for significance of perturbation response

clear all
DataFolder = ['X:\DATA\PROJECTS'];
FileNames = uigetfile_n_dir(DataFolder,'Select recording file');
load(FileNames{1,1});

DiffDRPer_total = array2table(cell(1,1),'VariableNames',{'FR_trial_4s'});
for i = 1:size(ProjectData,1) % loop through sessions
    trialSide_seconds = 1; 
    %% 1)
    % extract neural responses from ProjectData Table during trial periods
    clear SingletrialSpikes SpikeCount
    for j = 1:size(ProjectData.Trial_data{i,1},1)
        SingletrialSpikes = ProjectData.Trial_data{i,1}.SpikesTrial{j,1};
        edges = linspace(ProjectData.Session_data{i,1}.ACInfo{1,1}.trialStartsEnds(j,1)-trialSide_seconds,...
            ProjectData.Session_data{i,1}.ACInfo{1,1}.trialStartsEnds(j,2)+trialSide_seconds,...
            round(((ProjectData.Session_data{i,1}.ACInfo{1,1}.trialStartsEnds(j,2)+trialSide_seconds)-...
            (ProjectData.Session_data{i,1}.ACInfo{1,1}.trialStartsEnds(j,1)-trialSide_seconds))*60));
        for k = 1:size(SingletrialSpikes,2)
            [SpikeCount(j,k,1:length(edges)-1), Edges] = histcounts(SingletrialSpikes{k},edges);
        end
    end
    
    % extract neural responses from ProjectData Table during inter-trial periods
    % take all spikes of the recording of interest
    clear SessionSpikes
    for cs = 1:size(ProjectData.Units_Info{i,1},1)    
        SessionSpikes{cs} = (ProjectData.Units_Info{i,1}.Spiketimes{cs}(...
                ProjectData.Units_Info{i,1}.Spiketimes{cs}>ProjectData.Session_data{i,1}.MetaData{1,1}.lims(ProjectData.Session_data{i,1}.RecordingOI{1}-1) & ...
                ProjectData.Units_Info{i,1}.Spiketimes{cs}<(ProjectData.Session_data{i,1}.MetaData{1,1}.lims(ProjectData.Session_data{i,1}.RecordingOI{1})+...
                                   ProjectData.Session_data{i,1}.MetaData{1,1}.lims(ProjectData.Session_data{i,1}.RecordingOI{1}-1))) ...
                - ProjectData.Session_data{i,1}.MetaData{1,1}.lims(ProjectData.Session_data{i,1}.RecordingOI{1}-1))/ProjectData.Session_data{i,1}.SpikeInfo{1,1}{6,2};
    end
    % select only spikes in intertrial periods and compute firing rate
    clear SingletrialSpikes SpikeCount_ITI
    for j = 1:size(ProjectData.Trial_data{i,1},1)-1
            edges = linspace(ProjectData.Session_data{i,1}.ACInfo{1,1}.trialStartsEnds(j,2),...
                                ProjectData.Session_data{i,1}.ACInfo{1,1}.trialStartsEnds(j+1,1),...
                        round(((ProjectData.Session_data{i,1}.ACInfo{1,1}.trialStartsEnds(j+1,1))-...
                        (ProjectData.Session_data{i,1}.ACInfo{1,1}.trialStartsEnds(j,2)))*60));
        for k = 1:size(SessionSpikes,2)
            [SpikeCount_ITI{j,k}, Edges] = histcounts(SessionSpikes{k},edges);
        end
    end
    
    %% 2)
    % shuffle the identity of the trials and recompute the above metrics. Save
    % these values in an array
    OrigPertTrialOrder = ProjectData.Trial_data{i,1}.PerturbationON;
    clear DiffDRPer
    % smooth the FR
    BonvisionFR = (ProjectData.Session_data{1,1}.Time{1,1}(end)-ProjectData.Session_data{1,1}.Time{1,1}(1))/length(ProjectData.Session_data{1,1}.Time{1,1});
    GaussFilter_W = 0.6;
    Width = round(GaussFilter_W/BonvisionFR);
    Sigma = Width/3; % standard deviation in number of samples (converted from time in seconds)
    x_g = linspace(-Width/2, Width/2, Width);
    gaussFilter = exp(-x_g.^2/(2*Sigma^2));
    gaussFilter_ = gaussFilter / sum (gaussFilter); % normalize
    c1 = 1; clear PertStart_i SpikeCount_FR SpikeCount_FR_init SpikeCount_FR_pert SpikeCount_FR_postpert
    for ss = 1: size(SpikeCount,1)
        if ProjectData.Trial_data{i,1}.PerturbationON(ss,1)==1
            PertStart_i(c1,1) = min(find(ProjectData.Trial_data{i,1}.TF_in_trial{ss}==0)); % find index start of perturbation
            PertStart_i(c1,2) = max(find(ProjectData.Trial_data{i,1}.TF_in_trial{ss}==0)); % find  duration of perturbation
            PertStart_i(c1,3) = ss; % perturbation trial ID
            PertStart_i(c1,4) = mean(ProjectData.Trial_data{i,1}.Speed_in_trial{PertStart_i(c1,3)}(PertStart_i(c1,1):PertStart_i(c1,2))); % mean speed
            c1 = c1 + 1;
        end
    end
    for ss = 1: size(SpikeCount,1)    
        for rr = 1: size(SpikeCount,2)
            SpikeCount_FR(ss,rr,:) = conv(squeeze(SpikeCount(ss,rr,:)), gaussFilter_, 'same')*(1/BonvisionFR);
            SpikeCount_FR_init(ss,rr,:) = SpikeCount_FR(ss,rr,round(trialSide_seconds/BonvisionFR):round(trialSide_seconds/BonvisionFR)+60*4);
            if ProjectData.Trial_data{i,1}.PerturbationON(ss,1)==1
                SpikeCount_FR_prepert(ss,rr,:) = SpikeCount_FR(ss,rr,PertStart_i(PertStart_i(:,3)==ss,1)-60:PertStart_i(PertStart_i(:,3)==ss,1));
                SpikeCount_FR_pert(ss,rr,:) = SpikeCount_FR(ss,rr,PertStart_i(PertStart_i(:,3)==ss,1):PertStart_i(PertStart_i(:,3)==ss,1)+60);
                SpikeCount_FR_postpert(ss,rr,:) = SpikeCount_FR(ss,rr,PertStart_i(PertStart_i(:,3)==ss,2)+1:PertStart_i(PertStart_i(:,3)==ss,2)+1+60);
            else
                SpikeCount_FR_prepert(ss,rr,:) = SpikeCount_FR(ss,rr,round(mean(PertStart_i(:,1)))-60:round(mean(PertStart_i(:,1))));
                SpikeCount_FR_pert(ss,rr,:) = SpikeCount_FR(ss,rr,round(mean(PertStart_i(:,1))):round(mean(PertStart_i(:,1)))+60);
                SpikeCount_FR_postpert(ss,rr,:) = SpikeCount_FR(ss,rr,round(mean(PertStart_i(:,2)))+1:round(mean(PertStart_i(:,2)))+1+60);
            end
        end
    end
    clear SpikeCount_FR_ITI SpikeCount_FR_ITI_m
    for c1 = 1:size(SpikeCount_ITI,1)
        for c2 = 1:size(SpikeCount_ITI,2)
            SpikeCount_FR_ITI{c1,c2} = conv(SpikeCount_ITI{c1,c2},gaussFilter_, 'same')*(1/BonvisionFR);
            SpikeCount_FR_ITI_m(c1,c2) = mean(SpikeCount_FR_ITI{c1,c2});
        end
    end
    
    clear DiffDRPert DiffDRPert_post FR_vis FR_bl
    tic
    for p_i = 1:10 %1001
        clear FR_pert FR_postpert
        % compute the difference between FR response to perturbation and the same
        % period of the trials without the perturbation. Can compute other metrics
        % too, possibly.
        if p_i == 1
            PerturbationON = logical(OrigPertTrialOrder);
            PerturbationOFF = ~OrigPertTrialOrder;
        else   
            PerturbationON = logical(OrigPertTrialOrder(randperm(length(OrigPertTrialOrder))));
            PerturbationOFF = ~(OrigPertTrialOrder);
        end 
        % save the chunk of FR relative to three in-trial periods:
        % perturbation period
        FR_pert{1} = squeeze(sum(SpikeCount_FR_pert(PerturbationOFF,:,:),1) ...
            / sum(PerturbationOFF));
        FR_pert{2} = squeeze(sum(SpikeCount_FR_pert(PerturbationON,:,:),1) ...
            / sum(PerturbationON));
        FR_pert{3} = squeeze(sum(SpikeCount_FR_prepert(PerturbationOFF,:,:),1) ...
            / sum(PerturbationOFF));
        FR_pert{4} = squeeze(sum(SpikeCount_FR_prepert(PerturbationON,:,:),1) ...
            / sum(PerturbationON));
        % post perturbation period
        FR_postpert{1} = squeeze(sum(SpikeCount_FR_postpert(PerturbationOFF,:,:),1) ...
            / sum(PerturbationOFF));
        FR_postpert{2} = squeeze(sum(SpikeCount_FR_postpert(PerturbationON,:,:),1) ...
            / sum(PerturbationON));
        FR_postpert{3} = squeeze(sum(SpikeCount_FR_pert(PerturbationOFF,:,:),1) ...
            / sum(PerturbationOFF));
        FR_postpert{4} = squeeze(sum(SpikeCount_FR_pert(PerturbationON,:,:),1) ...
            / sum(PerturbationON));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%% 
        % from here it needs to be updated in order to "quantify" the pert
        % responses
        FR_pert_diff = abs(FR_pert{2}-FR_pert{4});
        FR_pert_diff_ = abs(FR_pert{2}-FR_pert{1});
        FR_postpert_diff = abs(FR_postpert{2}-FR_postpert{4});
        FR_postpert_diff_ = abs(FR_postpert{2}-FR_postpert{1});
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % shuffle the firing rate of every cell and re-compute FR over the first 4 seconds
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if p_i>=2
%             % take all the spiketimes and add an arbitrary time btw 20 and 60
%             % seconds. If the spiketimes fall over the end of the recording
%             % wrap them around.
%             clear SessionSpikes_sh
%             for st_c = 1:size(SessionSpikes,2)
%                 SessionSpikes_temp = SessionSpikes{st_c}+(20+rand*60);
%                 if sum(SessionSpikes_temp > ProjectData.Session_data{i,1}.ACInfo{1,1}.trialStartsEnds(end,2))>1
%                     SessionSpikes_temp(SessionSpikes_temp > ProjectData.Session_data{i,1}.ACInfo{1,1}.trialStartsEnds(end,2)) = ...
%                         SessionSpikes_temp(SessionSpikes_temp > ProjectData.Session_data{i,1}.ACInfo{1,1}.trialStartsEnds(end,2)) - ...
%                         ProjectData.Session_data{i,1}.ACInfo{1,1}.trialStartsEnds(end,2);
%                     SessionSpikes_sh{st_c} = SessionSpikes_temp;
%                 else
%                     SessionSpikes_sh{st_c} = SessionSpikes_temp;
%                 end
%             end
%             % compute the responses during the trials and during ITI's
%             clear SpikeCount_sh SpikeCount_ITI_sh
%             for j = 1:size(ProjectData.Trial_data{i,1},1)
%                 edges_trial = linspace(ProjectData.Session_data{i,1}.ACInfo{1,1}.trialStartsEnds(j,1),...
%                     ProjectData.Session_data{i,1}.ACInfo{1,1}.trialStartsEnds(j,2),...
%                     round(((ProjectData.Session_data{i,1}.ACInfo{1,1}.trialStartsEnds(j,2))-...
%                     (ProjectData.Session_data{i,1}.ACInfo{1,1}.trialStartsEnds(j,1)))*60));
%                 if j < size(ProjectData.Trial_data{i,1},1)
%                     edges_ITI = linspace(ProjectData.Session_data{i,1}.ACInfo{1,1}.trialStartsEnds(j,2),...
%                         ProjectData.Session_data{i,1}.ACInfo{1,1}.trialStartsEnds(j+1,1),...
%                         round(((ProjectData.Session_data{i,1}.ACInfo{1,1}.trialStartsEnds(j+1,1))-...
%                         (ProjectData.Session_data{i,1}.ACInfo{1,1}.trialStartsEnds(j,2)))*60));
%                 end
%                 for k = 1:size(SessionSpikes_sh,2)
%                     [SpikeCount_sh(j,k,1:length(edges_trial)-1), Edges] = histcounts(SessionSpikes_sh{k},edges_trial);
%                     if j < size(ProjectData.Trial_data{i,1},1)
%                         [SpikeCount_ITI_sh{j,k}, Edges] = histcounts(SessionSpikes_sh{k},edges_ITI);
%                     end
%                 end
%             end
%             clear SpikeCount_FR_sh SpikeCount_FR_init_sh
%             for ss = 1: size(SpikeCount,1)
%                 for rr = 1: size(SpikeCount,2)
%                     SpikeCount_FR_sh(ss,rr,:) = conv(squeeze(SpikeCount_sh(ss,rr,:)), gaussFilter_, 'same')*(1/BonvisionFR);
%                     if ss < size(SpikeCount,1)
%                         SpikeCount_FR_init_sh(ss,rr,:) = SpikeCount_FR_sh(ss,rr,1:60*4);
%                     end
%                 end
%             end
%             clear SpikeCount_FR_ITI_sh SpikeCount_FR_ITI_m_sh
%             for c1 = 1:size(SpikeCount_ITI_sh,1)
%                 for c2 = 1:size(SpikeCount_ITI_sh,2)
%                     SpikeCount_FR_ITI_sh{c1,c2} = conv(SpikeCount_ITI_sh{c1,c2},gaussFilter_, 'same')*(1/BonvisionFR);
%                     SpikeCount_FR_ITI_m_sh(c1,c2) = mean(SpikeCount_FR_ITI_sh{c1,c2});
%                 end
%             end
%             % first 4 seconds
%             FR_vis(p_i,:) = mean(squeeze(sum(SpikeCount_FR_init_sh(:,:,:)) ...
%                 / size(SpikeCount_FR_init_sh,1)),2);
%             % avg baseline FR
%             FR_bl(p_i,:) = mean(SpikeCount_FR_ITI_m_sh,1);
        elseif p_i == 1
            % first 4 seconds
            FR_vis(p_i,:) = mean(squeeze(sum(SpikeCount_FR_init(:,:,:)) ...
                        / size(SpikeCount_FR_init,1)),2);
            % avg baseline FR
            FR_bl(p_i,:) = mean(SpikeCount_FR_ITI_m,1);
        end
        
        % compute the difference for each cell
        DiffDRPert(p_i,:) = mean(FR_pert_diff,2);  % pert_ON_period - pre_pert_ON_periods
        DiffDRPert_(p_i,:) = mean(FR_pert_diff_,2); % pert_ON_period - pert_OFF_period
        DiffDRPert_post(p_i,:) = mean(FR_postpert_diff,2); % post_pert_ON_period - pert_ON_periods
        DiffDRPert_post_(p_i,:) = mean(FR_postpert_diff,2); % post_pert_ON_period - post_pert_OFF_period
        p_i
    end
    toc
    
    DiffDRPer_total.FR_trial_4s(i) = {FR_vis};                    % first 4 seconds
    DiffDRPer_total.FR_baseline(i) = {FR_bl};                     % avg baseline FR
    DiffDRPer_total.FR_pert(i) = {DiffDRPert};                            % perturbation period (before-during pert comparison)
    DiffDRPer_total.FR_postpert(i) = {DiffDRPert_post};                   % post-pert period (before-during pert comparison)
    DiffDRPer_total.FR_pert_(i) = {DiffDRPert_};                            % perturbation period (ON-OFF pert comparison)
    DiffDRPer_total.FR_postpert_(i) = {DiffDRPert_post_};                   % post-pert period (ON-OFF pert comparison)
    
%     DiffDRPer_total{1,i} = mean(FR_vis,2);      % first 4 seconds
%     DiffDRPer_total{2,i} = mean(SpikeCount_FR_ITI_m,2);      % avg baseline FR
%     DiffDRPer_total{3,i} = DiffDRPert;          % perturbation period
%     DiffDRPer_total{4,i} = DiffDRPert_post;     % post-pert period
%     
end

save('p_values.mat', 'DiffDRPer_total','-v7.3')

%% 3)
% compute how many times the responses are sign. higher than chance levels
k = 1; c1 = 1; c2 = 1; c3 = 1;
clear P_cell RespCell_pert_Index RespCell_postpert_Index
for i = 1:size(DiffDRPer_total,1) % loop through sessions
    for j = 1:size(DiffDRPer_total.FR_trial_4s{i,1},2)
        %P_cell(1,k) = 1-(sum(abs(DiffDRPer_total.FR_trial_4s{i,1}(1,j))>abs(DiffDRPer_total{2,i}(2:end,j)))/(length(DiffDRPer_total{2,i})-1)); % vis.stim. response
        
        P_cell(2,k) = 1-(sum(abs(DiffDRPer_total.FR_pert{i,1}(1,j))>abs(DiffDRPer_total.FR_pert{i,1}(2:end,j)))/(length(DiffDRPer_total.FR_pert{i,1})-1)); % perturbation
        if P_cell(2,k)<=0.01
            RespCell_pert_Index(c2,:) = [i, j]; 
            c2 = c2 + 1;
        end
        P_cell(3,k) = 1-(sum(abs(DiffDRPer_total.FR_postpert{i,1}(1,j))>abs(DiffDRPer_total.FR_postpert{i,1}(2:end,j)))/(length(DiffDRPer_total.FR_postpert{i,1})-1)); % post perturbation 
        if P_cell(3,k)<=0.01
            RespCell_postpert_Index(c3,:) = [i, j];
            c3 = c3 + 1;
        end
        k = k + 1;
        
    end
end

figure
% vis stim responsive
plot(find(P_cell(1,:)<=0.01),P_cell(1,find(P_cell(1,:)<=0.01)),'or')
hold on
plot(P_cell(1,:),'.k')


figure
% perturbation period
plot(find(P_cell(2,:)<=0.01),P_cell(2,find(P_cell(2,:)<=0.01)),'or')
hold on
plot(P_cell(2,:),'.k')

figure
% post-perturbation period
plot(find(P_cell(3,:)<=0.01),P_cell(3,find(P_cell(3,:)<=0.01)),'or')
hold on
plot(P_cell(3,:),'.k')


%% 4)
% make a graphical summary of the above points i.e. make a nice figure to
% show how you did it.

%% plot best example traces

%% plot example of the values of the diff. response with the shuffled
% responses
Rand_cell = abs(round(randn*length(RespCell_pert_Index)));
Diff_Responses = DiffDRPer_total.FR_pert{RespCell_pert_Index(Rand_cell,1),1}(1:end,RespCell_pert_Index(Rand_cell,2));
edges = 0:0.01:1;
[bincount, edges] = histcounts(Diff_Responses,edges);

figure
h = histogram(Diff_Responses);
xlabel('FR difference (Hz)'); ylabel('count');
hold on
stem(Diff_Responses(1),max(h.BinCounts)/2,'k');
box offfor normalised firing rate
set(gca,'TickDir','out')
title('Diff. perturbation response of shuffled trials');

%% plot Venn diagram or simple histogram with the total nr of cells
CellCount_vis = length(find(P_cell(1,:)<=0.01));
CellCount_pert = length(find(P_cell(2,:)<=0.01));
CellCount_pertpost = length(find(P_cell(3,:)<=0.01));
CellCount_pertANDpp = length(find(P_cell(2,:)<=0.01 & P_cell(3,:)<=0.01));




