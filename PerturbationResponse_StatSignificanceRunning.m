%% Tomaso Muzzu - UCL - 09/05/2019

% statistical test for significance of perturbation response runnning NOT running

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
    SpikeCount = SpikeCount.*60; % multiplied by 60 to get Hz
    
    
    % extract neural responses from ProjectData Table during inter-trial periods
    % take all spikes of the recording of interest
    SessionStart = ProjectData.Session_data{i,1}.ACInfo{1,1}.trialStartsEnds(1,1)-trialSide_seconds;
    SessionEnd = ProjectData.Session_data{i,1}.ACInfo{1,1}.trialStartsEnds(end,2)+trialSide_seconds;
    Time_edges = linspace(SessionStart, SessionEnd, (SessionEnd-SessionStart)*60); % 60 is sampling frequency.

    clear SessionSpikes
    for cs = 1:size(ProjectData.Units_Info{i,1},1)    
         SessionSpikes{cs} = ProjectData.Units_Info{i,1}.Spiketimes{cs,1}(...
                ProjectData.Units_Info{i,1}.Spiketimes{cs,1}>=SessionStart & ...
                ProjectData.Units_Info{i,1}.Spiketimes{cs,1}<=SessionEnd);
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
    clear DiffDRPert
    % smooth the FR
    BonvisionFR = (ProjectData.Session_data{1,1}.Time{1,1}(end)-ProjectData.Session_data{1,1}.Time{1,1}(1))/length(ProjectData.Session_data{1,1}.Time{1,1});
    GaussFilter_W = 0.6;
    Width = round(GaussFilter_W/BonvisionFR);
    Sigma = Width/3; % standard deviation in number of samples (converted from time in seconds)
    x_g = linspace(-Width/2, Width/2, Width);
    gaussFilter = exp(-x_g.^2/(2*Sigma^2));
    gaussFilter_ = gaussFilter / sum (gaussFilter); % normalize
    c1 = 1; clear PertStart_i SpikeCount_FR SpikeCount_FR_init SpikeCount_FR_pert SpikeCount_FR_postpert SpikeCount_FR_prepert
    for ss = 1: size(SpikeCount,1)
        if ProjectData.Trial_data{i,1}.PerturbationON(ss,1)==1
            PertStart_i(c1,1) = min(find(ProjectData.Trial_data{i,1}.TF_in_trial{ss}==0)); % find index start of perturbation
            PertStart_i(c1,2) = max(find(ProjectData.Trial_data{i,1}.TF_in_trial{ss}==0)); % find  duration of perturbation
            PertStart_i(c1,3) = ss; % perturbation trial ID
            PertStart_i(c1,4) = mean(ProjectData.Trial_data{i,1}.Speed_in_trial{PertStart_i(c1,3)}(PertStart_i(c1,1):PertStart_i(c1,2))); % mean speed
            c1 = c1 + 1;
        end
    end
    % save mean speeds over the three periods of interest: initial 4
    % seconds, perturbation time window, 1s post perturbation
    OrigPertTrialOrder = ProjectData.Trial_data{i,1}.PerturbationON;
    clear TrialSpeeds
    for ss = 1:size(SpikeCount,1)
        TrialSpeeds(ss,:) = [ mean(ProjectData.Trial_data{i,1}.Speed_in_trial{ss,1}(round(trialSide_seconds/BonvisionFR):round(trialSide_seconds/BonvisionFR)+60*4)),... % initial 4 s
                              mean(ProjectData.Trial_data{i,1}.Speed_in_trial{ss,1}(round(mean(PertStart_i(:,1))):round(mean(PertStart_i(:,2))))),...
                              mean(ProjectData.Trial_data{i,1}.Speed_in_trial{ss,1}(round(mean(PertStart_i(:,2))):round(mean(PertStart_i(:,2))+60)))];
    end
    for ss = 1: size(SpikeCount,1)    
        for rr = 1: size(SpikeCount,2)
            SpikeCount_FR(ss,rr,:) = conv(squeeze(SpikeCount(ss,rr,:)), gaussFilter_, 'same');
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
    clear DiffDRPert_1 DiffDRPert_2 DiffDRPert_3 DiffDRPert_4 
    clear DiffDRPert_1a DiffDRPert_2a DiffDRPert_3a DiffDRPert_4a DiffDRPert_area
    clear DiffDRPert_post_1 DiffDRPert_post_2 DiffDRPert_post_3 DiffDRPert_post_4 
    clear DiffDRPert_post_1a DiffDRPert_post_2a DiffDRPert_post_3a DiffDRPert_post_4a DiffDRPert_post_area

    tic
    for p_i = 1:1001
        clear FR_pert FR_postpert
        % compute the difference between FR response to perturbation and the same
        % period of the trials without the perturbation. Can compute other metrics
        % too, possibly.
        if p_i == 1
            PerturbationON = logical(OrigPertTrialOrder);
            PerturbationOFF = ~OrigPertTrialOrder;
            % running trials - with pert and w/o pert
            PerturbationON_run = find(logical(OrigPertTrialOrder) & TrialSpeeds(:,2)>3);
            PerturbationOFF_run = find(~OrigPertTrialOrder & TrialSpeeds(:,2)>3);
            % still trials - with pert and w/o pert
            PerturbationON_still = find(logical(OrigPertTrialOrder) & TrialSpeeds(:,2)<3);
            PerturbationOFF_still = find(~OrigPertTrialOrder & TrialSpeeds(:,2)<3);
        else
            PerturbationON = logical(OrigPertTrialOrder(randperm(length(OrigPertTrialOrder))));
            PerturbationOFF = ~(OrigPertTrialOrder);
            % running trials - with pert and w/o pert - SHUFFLED
            trials_temp = [PerturbationON_run; PerturbationOFF_run];
            [PerturbationON_run, idx] = datasample(trials_temp, length(PerturbationON_run) ,'Replace',false);
            trials_temp(idx)=[]; PerturbationOFF_run = trials_temp;
            % still trials - with pert and w/o pert - SHUFFLED
            trials_temp = [PerturbationON_still; PerturbationOFF_still];
            [PerturbationON_still, idx] = datasample(trials_temp, length(PerturbationON_still) ,'Replace',false);
            trials_temp(idx)=[]; PerturbationOFF_still = trials_temp;
        end 
                
        % save the chunk of FR relative to three in-trial periods:
        % perturbation period
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % trials w/o perturbation
        FR_pert{1,1} = squeeze(sum(SpikeCount_FR_pert(PerturbationOFF,:,:),1) ...
            / sum(PerturbationOFF));
        FR_pert{1,2} = squeeze(std(SpikeCount_FR_pert(PerturbationOFF,:,:),1)) ...
            / sqrt(sum(PerturbationOFF));
        FR_pert{1,3} = squeeze(sum(SpikeCount_FR_prepert(PerturbationOFF,:,:),1) ...
            / sum(PerturbationOFF)); 
        FR_pert{1,4} = squeeze(std(SpikeCount_FR_prepert(PerturbationOFF,:,:),1)) ...
            / sqrt(sum(PerturbationOFF));
        % all trials with perturbation
        FR_pert{2,1} = squeeze(sum(SpikeCount_FR_pert(PerturbationON,:,:),1) ...
            / sum(PerturbationON));
        FR_pert{2,2} = squeeze(std(SpikeCount_FR_pert(PerturbationON,:,:),1)) ...
            / sqrt(sum(PerturbationON));
        FR_pert{2,3} = squeeze(sum(SpikeCount_FR_prepert(PerturbationON,:,:),1) ...
            / sum(PerturbationON)); 
        FR_pert{2,4} = squeeze(std(SpikeCount_FR_prepert(PerturbationON,:,:),1)) ...
            / sqrt(sum(PerturbationON));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % running only
        FR_pert{3,1} = squeeze(sum(SpikeCount_FR_pert(PerturbationOFF_run,:,:),1) ...
            / length(PerturbationOFF_run));
        FR_pert{3,2} = squeeze(std(SpikeCount_FR_pert(PerturbationOFF_run,:,:),1)) ...
            / sqrt(length(PerturbationOFF_run));
        FR_pert{3,3} = squeeze(sum(SpikeCount_FR_prepert(PerturbationOFF_run,:,:),1) ...
            / length(PerturbationOFF_run));
        FR_pert{3,4} = squeeze(std(SpikeCount_FR_prepert(PerturbationOFF_run,:,:),1)) ...
            / sqrt(length(PerturbationOFF_run));
        % trials with perturbation - running only
        FR_pert{4,1} = squeeze(sum(SpikeCount_FR_pert(PerturbationON_run,:,:),1) ...
            / length(PerturbationON_run));
        FR_pert{4,2} = squeeze(std(SpikeCount_FR_pert(PerturbationON_run,:,:),1)) ...
            / sqrt(length(PerturbationON_run));
        FR_pert{4,3} = squeeze(sum(SpikeCount_FR_prepert(PerturbationON_run,:,:),1) ...
            / length(PerturbationON_run));
        FR_pert{4,4} = squeeze(std(SpikeCount_FR_prepert(PerturbationON_run,:,:),1)) ...
            / sqrt(length(PerturbationON_run));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % still only
        FR_pert{5,1} = squeeze(sum(SpikeCount_FR_pert(PerturbationOFF_still,:,:),1) ...
            / length(PerturbationOFF_still));
        FR_pert{5,2} = squeeze(std(SpikeCount_FR_pert(PerturbationOFF_still,:,:),1)) ...
            / sqrt(length(PerturbationOFF_still));
        FR_pert{5,3} = squeeze(sum(SpikeCount_FR_prepert(PerturbationOFF_still,:,:),1) ...
            / length(PerturbationOFF_still));
        FR_pert{5,4} = squeeze(std(SpikeCount_FR_prepert(PerturbationOFF_still,:,:),1)) ...
            / sqrt(length(PerturbationOFF_still));
        if ~isempty(PerturbationON_still)
            FR_pert{6,1} = squeeze(sum(SpikeCount_FR_pert(PerturbationON_still,:,:),1) ...
                / length(PerturbationON_still));
            FR_pert{6,2} = squeeze(std(SpikeCount_FR_pert(PerturbationON_still,:,:),1)) ...
                / sqrt(length(PerturbationON_still));
            FR_pert{6,3} = squeeze(sum(SpikeCount_FR_prepert(PerturbationON_still,:,:),1) ...
                / length(PerturbationON_still));
            FR_pert{6,4} = squeeze(std(SpikeCount_FR_prepert(PerturbationON_still,:,:),1)) ...
                / sqrt(length(PerturbationON_still));
        else
            FR_pert{6,1} = 999*ones(size(FR_pert{3,1},1),size(FR_pert{3,1},2));
            FR_pert{6,2} = zeros(size(FR_pert{4,1},1),size(FR_pert{4,1},2));
            FR_pert{6,3} = 999*ones(size(FR_pert{3,1},1),size(FR_pert{3,1},2));
            FR_pert{6,4} = zeros(size(FR_pert{4,1},1),size(FR_pert{4,1},2));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        % post perturbation period
        FR_postpert{1,1} = squeeze(sum(SpikeCount_FR_postpert(PerturbationOFF,:,:),1) ...
            / sum(PerturbationOFF));
        FR_postpert{1,2} = squeeze(std(SpikeCount_FR_postpert(PerturbationOFF,:,:),1)) ...
            / sqrt(sum(PerturbationOFF));
        % all trials with perturbation
        FR_postpert{2,1} = squeeze(sum(SpikeCount_FR_postpert(PerturbationON,:,:),1) ...
            / sum(PerturbationON));
        FR_postpert{2,2} = squeeze(std(SpikeCount_FR_postpert(PerturbationON,:,:),1)) ...
            / sqrt(sum(PerturbationON));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % running only
        FR_postpert{3,1} = squeeze(sum(SpikeCount_FR_postpert(PerturbationOFF_run,:,:),1) ...
            / length(PerturbationOFF_run));
        FR_postpert{3,2} = squeeze(std(SpikeCount_FR_postpert(PerturbationOFF_run,:,:),1)) ...
            / sqrt(length(PerturbationOFF_run));
        FR_postpert{4,1} = squeeze(sum(SpikeCount_FR_postpert(PerturbationON_run,:,:),1) ...
            / length(PerturbationON_run));
        FR_postpert{4,2} = squeeze(std(SpikeCount_FR_postpert(PerturbationON_run,:,:),1)) ...
            / sqrt(length(PerturbationON_run));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % still only
        FR_postpert{5,1} = squeeze(sum(SpikeCount_FR_postpert(PerturbationOFF_still,:,:),1) ...
            / length(PerturbationOFF_still));
        FR_postpert{5,2} = squeeze(std(SpikeCount_FR_postpert(PerturbationOFF_still,:,:),1)) ...
            / sqrt(length(PerturbationOFF_still));
        if ~isempty(PerturbationON_still)
            FR_postpert{6,1} = squeeze(sum(SpikeCount_FR_postpert(PerturbationON_still,:,:),1) ...
                / length(PerturbationON_still));
            FR_postpert{6,2} = squeeze(std(SpikeCount_FR_postpert(PerturbationON_still,:,:),1)) ...
                / sqrt(length(PerturbationON_still));
        else
            FR_postpert{6,1} = 999*ones(size(FR_postpert{3,1},1),size(FR_postpert{3,1},2));
            FR_postpert{6,2} = zeros(size(FR_postpert{4,1},1),size(FR_postpert{4,1},2));
            FR_postpert{6,3} = 999*ones(size(FR_postpert{3,1},1),size(FR_postpert{3,1},2));
            FR_postpert{6,4} = zeros(size(FR_postpert{4,1},1),size(FR_postpert{4,1},2));
        end
        
        clear FR_pert_diff_1 FR_pert_diff_1a FR_pert_diff_2 FR_pert_diff_2a
        clear FR_pert_diff_3 FR_pert_diff_3a FR_pert_diff_4 FR_pert_diff_4a FR_pert_area
        clear FR_postpert_diff_1 FR_postpert_diff_2 FR_postpert_diff_3 FR_postpert_diff_4
        clear FR_postpert_diff_1a FR_postpert_diff_2a FR_postpert_diff_3a FR_postpert_diff_4a FR_postpert_area
        FR_pert_diff_1 = abs(FR_pert{2,1}-FR_pert{1,1}); % diff btw pert and no-pert trials
        FR_pert_diff_1a = abs(FR_pert{2,1}-FR_pert{2,3}); % diff btw pert and pre-pert trials
        
        FR_pert_diff_2 = abs(FR_pert{4,1}-FR_pert{3,1}); % diff btw run-pert and no-pert trials
        FR_pert_diff_2a = abs(FR_pert{4,1}-FR_pert{4,3}); % diff btw run-pert and run pre-pert trials
        
        FR_pert_diff_3 = abs(FR_pert{6,1}-FR_pert{5,1}); % diff btw still-pert and no-pert trials
        FR_pert_diff_3a = abs(FR_pert{6,1}-FR_pert{6,3}); % diff btw still-pert and still pre-pert trials
        
        FR_pert_diff_4 = abs(FR_pert{6,1}-FR_pert{4,1}); % diff btw run-pert and still-pert trials
        FR_pert_diff_4a = abs(FR_pert_diff_2a-FR_pert_diff_3a); % diff btw run-pert and still-pert trials
        
        FR_pert_area = abs(trapz(FR_pert{2},2)-trapz(FR_pert{1},2));
        
        FR_postpert_diff_1 = abs(FR_postpert{2,1}-FR_postpert{1,1}); % diff btw pert and no-pert trials
        FR_postpert_diff_1a = abs(FR_postpert{2,1}-FR_pert{2,1}); % diff btw pert and pre-pert trials
        
        FR_postpert_diff_2 = abs(FR_postpert{4,1}-FR_postpert{3,1}); % diff btw run-pert and no-pert trials
        FR_postpert_diff_2a = abs(FR_postpert{4,1}-FR_pert{4,1}); % diff btw run-pert and run pre-pert trials
        
        FR_postpert_diff_3 = abs(FR_postpert{6,1}-FR_postpert{5,1}); % diff btw still-pert and no-pert trials
        FR_postpert_diff_3a = abs(FR_postpert{6,1}-FR_pert{6,1}); % diff btw still-pert and still pre-pert trials
        
        FR_postpert_diff_4 = abs(FR_postpert{6,1}-FR_postpert{4,1}); % diff btw run-pert and still-pert trials
        FR_postpert_diff_4a = abs(FR_postpert_diff_3a-FR_postpert_diff_2a); % diff btw run-pert and still-pert trials
        
        FR_postpert_area = abs(trapz(FR_postpert{2},2)-trapz(FR_postpert{1},2));
               
        % compute the difference for each cell
        DiffDRPert_1(p_i,:) = mean(FR_pert_diff_1,2);
        DiffDRPert_1a(p_i,:) = mean(FR_pert_diff_1a,2);
        DiffDRPert_2(p_i,:) = mean(FR_pert_diff_2,2);
        DiffDRPert_2a(p_i,:) = mean(FR_pert_diff_2a,2);
        DiffDRPert_3(p_i,:) = mean(FR_pert_diff_3,2);
        DiffDRPert_3a(p_i,:) = mean(FR_pert_diff_3a,2);
        DiffDRPert_4(p_i,:) = mean(FR_pert_diff_4,2);
        DiffDRPert_4a(p_i,:) = mean(FR_pert_diff_4a,2);
        DiffDRPert_area(p_i,:) = mean(FR_pert_area,2);
        
        DiffDRPert_post_1(p_i,:) = mean(FR_postpert_diff_1,2);
        DiffDRPert_post_1a(p_i,:) = mean(FR_postpert_diff_1a,2);
        DiffDRPert_post_2(p_i,:) = mean(FR_postpert_diff_2,2);
        DiffDRPert_post_2a(p_i,:) = mean(FR_postpert_diff_2a,2);
        DiffDRPert_post_3(p_i,:) = mean(FR_postpert_diff_3,2);
        DiffDRPert_post_3a(p_i,:) = mean(FR_postpert_diff_3a,2);
        DiffDRPert_post_4(p_i,:) = mean(FR_postpert_diff_4,2);
        DiffDRPert_post_4a(p_i,:) = mean(FR_postpert_diff_4a,2);
        DiffDRPert_post_area(p_i,:) = mean(FR_postpert_area,2);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % shuffle the firing rate of every cell and re-compute FR over the first 4 seconds
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         if p_i>=2
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
%         else
            if p_i == 1
            % first 4 seconds
            FR_vis(p_i,:) = mean(squeeze(sum(SpikeCount_FR_init(:,:,:)) ...
                        / size(SpikeCount_FR_init,1)),2);
            % avg baseline FR
            FR_bl(p_i,:) = mean(SpikeCount_FR_ITI_m,1);
        end
        
        p_i
    end
    toc
    
    DiffDRPer_total.FR_trial_4s(i) = {FR_vis};                    % first 4 seconds
    DiffDRPer_total.FR_baseline(i) = {FR_bl};                     % avg baseline FR
    DiffDRPer_total.FR_pert1(i) = {[DiffDRPert_1 DiffDRPert_1a]};                            % perturbation period
    DiffDRPer_total.FR_pert2(i) = {[DiffDRPert_2 DiffDRPert_2a]};                            % pert run - no pert
    DiffDRPer_total.FR_pert3(i) = {[DiffDRPert_3 DiffDRPert_3a]};                            % pert still -  no pert
    DiffDRPer_total.FR_pert4(i) = {[DiffDRPert_4 DiffDRPert_4a]};                            % pert run - pert still
    DiffDRPer_total.FR_pertArea(i) = {DiffDRPert_area};                            % pert area
    DiffDRPer_total.FR_postpert1(i) = {[DiffDRPert_post_1 DiffDRPert_post_1a]};                   % postpert period
    DiffDRPer_total.FR_postpert2(i) = {[DiffDRPert_post_2 DiffDRPert_post_2a]};                   % postpert run - no pert
    DiffDRPer_total.FR_postpert3(i) = {[DiffDRPert_post_3 DiffDRPert_post_3a]};                   % postpert still -  no pert
    DiffDRPer_total.FR_postpert4(i) = {[DiffDRPert_post_4 DiffDRPert_post_4a]};                   % postpert run - postpert still
    DiffDRPer_total.FR_postpertArea(i) = {DiffDRPert_post_area};                   % postpert area
    DiffDRPer_total.TrialTypes(i) = {[sum(PerturbationON), sum(PerturbationOFF), ...
                                      length(PerturbationON_run), length(PerturbationOFF_run),...
                                      length(PerturbationON_still), length(PerturbationOFF_still)]};                   % n for the conditions
    i
    pause(0.5)
%     DiffDRPer_total{1,i} = mean(FR_vis,2);      % first 4 seconds
%     DiffDRPer_total{2,i} = mean(SpikeCount_FR_ITI_m,2);      % avg baseline FR
%     DiffDRPer_total{3,i} = DiffDRPert;          % perturbation period
%     DiffDRPer_total{4,i} = DiffDRPert_post;     % post-pert period
%     
end

save('p_values_Sept.mat', 'DiffDRPer_total','-v7.3')

%% 3)
% compute how many times the responses are sign. higher than chance levels
k = 1; c1 = 1; c2 = 1; c3 = 1; c4 = 1; c5 = 1; c6 = 1; c7 = 1; c8 = 1;
clear P_cell RespCell_pert_Index RespCell_postpert_Index RespCell_VisStim_Index RespCell_pert_Index RespCell_postpert_Index 
clear RespCell_pert1_Index RespCell_pert2_Index RespCell_pert3_Index
clear RespCell_postpert1_Index RespCell_postpert2_Index
for i = 1:size(DiffDRPer_total,1) % loop through sessions
    for j = 1:size(DiffDRPer_total.FR_trial_4s{i,1},2)
        FR_diff = (abs(DiffDRPer_total_vis.FR_trial_4s{i,1}(1,j))-abs(DiffDRPer_total_vis.FR_baseline{i,1}(1,j)));
        FR_diff_all = (abs(DiffDRPer_total_vis.FR_trial_4s{i,1}(2:end,j))-abs(DiffDRPer_total_vis.FR_baseline{i,1}(2:end,j)));
        P_cell(1,k) = 1-(sum(FR_diff>FR_diff_all)/(length(DiffDRPer_total_vis.FR_trial_4s{i,1})-1)); % vis.stim. response
        if P_cell(1,k)<=0.01
            RespCell_VisStim_Index(c1,:) = [i, j]; 
            c1 = c1 + 1;
        end
        
        
        % all pert trials
        P_cell(2,k) = 1-(sum(abs(DiffDRPer_total.FR_pert{i,1}(1,j))>abs(DiffDRPer_total.FR_pert{i,1}(2:end,j)))/(length(DiffDRPer_total.FR_pert{i,1})-1)); % perturbation
        if P_cell(2,k)<=0.01
            RespCell_pert_Index(c2,:) = [i, j]; 
            c2 = c2 + 1;
        end
        % 
        P_cell(3,k) = 1-(sum(abs(DiffDRPer_total.FR_postpert{i,1}(1,j))>abs(DiffDRPer_total.FR_postpert{i,1}(2:end,j)))/(length(DiffDRPer_total.FR_postpert{i,1})-1)); % post perturbation
        if P_cell(3,k)<=0.01
            RespCell_postpert_Index(c3,:) = [i, j];
            c3 = c3 + 1;
        end
        
        %% perturbation period
        % moving only
        P_cell(4,k) = 1-(sum(abs(DiffDRPer_total.FR_pert1{i,1}(1,j))>abs(DiffDRPer_total.FR_pert1{i,1}(2:end,j)))/(length(DiffDRPer_total.FR_pert1{i,1})-1)); % post perturbation 
        if P_cell(4,k)<=0.01
            RespCell_pert1_Index(c4,:) = [i, j];
            c4 = c4 + 1;
        end
        
        % still only
        P_cell(5,k) = 1-(sum(abs(DiffDRPer_total.FR_pert2{i,1}(1,j))>abs(DiffDRPer_total.FR_pert2{i,1}(2:end,j)))/(length(DiffDRPer_total.FR_pert2{i,1})-1)); % post perturbation 
        if P_cell(5,k)<=0.01
            RespCell_pert2_Index(c5,:) = [i, j];
            c5 = c5 + 1;
        end
        
        % difference between running and still pert trials
        P_cell(8,k) = 1-(sum(abs(DiffDRPer_total.FR_pert3{i,1}(1,j))>abs(DiffDRPer_total.FR_pert3{i,1}(2:end,j)))/(length(DiffDRPer_total.FR_pert3{i,1})-1)); % post perturbation 
        if P_cell(8,k)<=0.01
            RespCell_pert3_Index(c5,:) = [i, j];
            c8 = c8 + 1;
        end
        
        %% post perturbation
        % moving only
        P_cell(6,k) = 1-(sum(abs(DiffDRPer_total.FR_postpert1{i,1}(1,j))>abs(DiffDRPer_total.FR_postpert1{i,1}(2:end,j)))/(length(DiffDRPer_total.FR_postpert1{i,1})-1)); % post perturbation 
        if P_cell(6,k)<=0.01
            RespCell_postpert1_Index(c6,:) = [i, j];
            c6 = c6 + 1;
        end
        
        % still only
        P_cell(7,k) = 1-(sum(abs(DiffDRPer_total.FR_postpert2{i,1}(1,j))>abs(DiffDRPer_total.FR_postpert2{i,1}(2:end,j)))/(length(DiffDRPer_total.FR_postpert2{i,1})-1)); % post perturbation 
        if P_cell(7,k)<=0.01
            RespCell_postpert2_Index(c7,:) = [i, j];
            c7 = c7 + 1;
        end
        UnitsList(k,:) = [i, j];
        k = k + 1;
    end
end

figure
% vis stim responsive
plot(find(P_cell(1,:)<=0.01),P_cell(1,find(P_cell(1,:)<=0.01)),'or')
hold on
plot(P_cell(1,:),'.k')

figure
plot(find(P_cell(2,:)<=0.01),P_cell(2,find(P_cell(2,:)<=0.01)),'or')
hold on
plot(P_cell(2,:),'.k')
title('pert. all trials')

figure
plot(find(P_cell(3,:)<=0.01),P_cell(3,find(P_cell(3,:)<=0.01)),'or')
hold on
plot(P_cell(3,:),'.k')
title('post pert all trials')

figure
plot(find(P_cell(4,:)<=0.01),P_cell(4,find(P_cell(4,:)<=0.01)),'or')
hold on
plot(P_cell(4,:),'.k')
title('running only')


figure
plot(find(P_cell(5,:)<=0.01),P_cell(5,find(P_cell(5,:)<=0.01)),'or')
hold on
plot(P_cell(5,:),'.k')
title('still only')

figure
plot(find(P_cell(8,:)<=0.01),P_cell(8,find(P_cell(8,:)<=0.01)),'or')
hold on
plot(P_cell(8,:),'.k')
title('run-still only')

figure
histogram(P_cell(8,:),100)

%% 4)
% make a graphical summary of the above points i.e. make a nice figure to
% show how you did it.


%% Venn diagram of responsive cells: stimulus, perturbation, postperturbation
idx_vstim = find(P_cell(1,:)<=0.01);
idx_pert = find(P_cell(2,:)<=0.01);
idx_postpert = find(P_cell(3,:)<=0.01);
idx_pert_run = find(P_cell(4,:)<=0.01);
idx_pert_still = find(P_cell(5,:)<=0.01);
idx_postpert_run = find(P_cell(6,:)<=0.01);
idx_postpert_still = find(P_cell(7,:)<=0.01);

I_1 = intersect(idx_vstim,idx_pert);
I_2 = intersect(idx_vstim,idx_postpert);
I_3 = intersect(idx_pert,idx_postpert);
I_4 = intersect(I_1,idx_postpert);

I_neg = setdiff(idx_vstim,idx_pert); % it is in vis but not pert
I_1_neg = intersect(I_neg,idx_postpert); % vis ~pert ~post-pert

I_neg = setdiff(idx_pert,idx_vstim); % it is in pert but not vis stim
I_2_neg = setdiff(I_neg,idx_postpert); % ~vis pert ~post-pertt

I_neg = setdiff(idx_postpert,idx_pert); % it is in post-pert but not pert
I_3_neg = intersect(I_neg,idx_vstim); % ~vis ~pert post-pertt


%Plot a simple 3-circle venn diagram with custom patch properties
UnitGroups = [length(idx_vstim) length(idx_pert) length(idx_postpert)]; 
I = [length(I_1) length(I_2) length(I_3) length(I_4)]; % [i12 i13 i23 i123]
figure, axis equal, axis off
set(0, 'DefaultFigureRenderer', 'opengl')
venn(UnitGroups,I,'EdgeColor','black','ErrMinMode', 'ChowRodgers')

Recordings_visStim = length(unique(RespCell_VisStim_Index(:,1)));
Recordings_visStim = length(unique(RespCell_pert_Index(:,1)));
Recordings_visStim = length(unique(RespCell_postpert_Index(:,1)));

Recordings_pertrun = length(unique(RespCell_pert1_Index(:,1)));
Recordings_pertstill = length(unique(RespCell_pert2_Index(:,1)));
Recordings_postpertrun = length(unique(RespCell_postpert1_Index(:,1)));
Recordings_postpertstill = length(unique(RespCell_postpert2_Index(:,1)));

% TO DO: 
% - Examples of cell responsive to various stimuli
% visual stimulus only
Units_idx = UnitsList(I_1_neg,:);
value2show = 1; smoothing = 9;
ShowExampleCells(ProjectData,DiffDRPer_total,value2show,Units_idx,smoothing)
% pert only
Units_idx = UnitsList(I_2_neg,:);
value2show = 1; smoothing = 9;
ShowExampleCells(ProjectData,DiffDRPer_total,value2show,Units_idx,smoothing)
% post pert only
Units_idx = UnitsList(I_3_neg,:);
value2show = 1; smoothing = 9;
ShowExampleCells(ProjectData,DiffDRPer_total,value2show,Units_idx,smoothing)
% All Cells together
Units_idx = UnitsList(I_4,:);
value2show = 1; smoothing = 9;
ShowExampleCells(ProjectData,DiffDRPer_total,value2show,Units_idx,smoothing)

% - heat map of all the cells together
value2show = 1; % 1 for normalised firing rate
AllCellResponse(ProjectData,DiffDRPer_total,value2show)

% - heat map of all the cells together differentiating btw run and still 
cond = 0; % 1=run, 0=still, 2=run-diff
value2show = 1;
AllCellResponse_runstill(ProjectData,DiffDRPer_total,value2show,cond)
% this function is plotting the data normalised.

% - heat map of all the cells together distinguishing the directions
value2show = 1; % 0 for pert-nonpert trial,  1 for normalised firing rate of perturbation only
Directions = unique(ProjectData.Trial_data{1,1}.Direction);
for d = 2:length(Directions)
    AllCellsResponse_Direction(ProjectData,DiffDRPer_total,value2show,Directions(d))
end

% example responses during running and stationary
Units_idx = UnitsList(I_4,:);
value2show = 1; smoothing = 9;
ShowExampleCells_RunStill(ProjectData,DiffDRPer_total,value2show,[17*ones(1,64);1:64]' ,smoothing);

% units from recordings with at least 10 trials run/still in pert
% conditions: 
[1*ones(1,42);1:42]
[5*ones(1,41);1:41]
[12*ones(1,22);1:22]
[15*ones(1,30);1:30]
[17*ones(1,64);1:64]

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
box off
set(gca,'TickDir','out')
title('Diff. perturbation response of shuffled trials');

%% plot cell responses to perturbation distinguishing run & stationary
value2show = 1; % 1 for normalised firing rate
cond = 1; % 0 for still, 1 for run, 2 for run-still
AllCellResponse_runstill(ProjectData,DiffDRPer_total,value2show,cond)

%% plot cell responses to perturbation distinguishing grating directions
Units_idx = UnitsList(I_1_neg,:);
value2show = 1; smoothing = 9;
ShowExampleCells_GratingDir(ProjectData,DiffDRPer_total,value2show,Units_idx(1,:) ,smoothing);

%% compute modulation indexes for every 
value2show = 1; % 1 for normalised firing rate
ModulationIndex_AllCells(ProjectData,DiffDRPer_total,value2show,cond)


