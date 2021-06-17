%% Tomaso Muzzu - UCL - 16/05/2019

% function to plot neural responses of units to perturbation trials
% depending on grating direction

function ModulationIndex_AllCellsRunning(ProjectData,DiffDRPer_total,value2show,cond)

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
clear PertRespMean_on PertRespMean_off VisRespMean RespSummay_total
RespSummay_total = array2table(cell(1,7),...
                                'VariableNames',{'PertRespMean_run', ...
                                                    'PertRespMean_still', ...
                                                    'PertRespMean_off_run', ...
                                                    'PertRespMean_off_still',...
                                                    'VisRespMean_run',...
                                                    'VisRespMean_still',...
                                                    'Trial_types'});
for rr = 1:size(ProjectData,1) % scan through sessions
    clear NormFR NormFR_sem NormFR_off NormFR_sem_off
    for i = 1:size(ProjectData.Units_Info{rr,1},1) % scan through units
        clear SingletrialSpikes SpikeCount_temp SpikeCount_smooth edges PertRespMeanCell
        trialcount = 1;
        for j = 1:size(ProjectData.Trial_data{rr,1},1) % scane through trials
            SingletrialSpikes = ProjectData.Trial_data{rr,1}.SpikesTrial{j,1}{1,i};
            TrialStart = ProjectData.Session_data{rr,1}.ACInfo{1,1}.trialStartsEnds(j,1);
            TrialEnd = ProjectData.Session_data{rr,1}.ACInfo{1,1}.trialStartsEnds(j,2);
            edges = linspace(TrialStart-trialSide_seconds,TrialEnd+trialSide_seconds,((TrialEnd+trialSide_seconds)-(TrialStart-trialSide_seconds))*60);
            [SpikeCount_temp, Edges] = histcounts(SingletrialSpikes,edges);
            SpikeCount_smooth(j,1:length(SpikeCount_temp)) = conv(SpikeCount_temp, gaussFilter_, 'same');
            if ProjectData.Trial_data{rr,1}.PerturbationON(j,1)
                Pert_idx_start(trialcount) = min(find(ProjectData.Trial_data{rr,1}.TF_in_trial{j,:}==0));
                Pert_idx_end(trialcount) = max(find(ProjectData.Trial_data{rr,1}.TF_in_trial{j,:}==0));
                PertRespMeanCell(trialcount,:,1) = SpikeCount_smooth(j,Pert_idx_start(trialcount):Pert_idx_start(trialcount)+60); % response during pert. for ranking all response
                PertRespMeanCell(trialcount,:,2) = SpikeCount_smooth(j,Pert_idx_end(trialcount)+1:Pert_idx_end(trialcount)+1+60); % response post perturbation
                PertRespMeanCell(trialcount,:,3) = SpikeCount_smooth(j,Pert_idx_start(trialcount)-60:Pert_idx_start(trialcount)); % response pre perturbation
                trialcount = trialcount + 1;
            end
        end
        
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
                TrialSpeeds(ss,:) = [ mean(ProjectData.Trial_data{rr,1}.Speed_in_trial{ss,1}(round(trialSide_seconds*BonvisionFR)-30:round(trialSide_seconds*BonvisionFR)+60*4)),... % initial 4 s
                                      mean(ProjectData.Trial_data{rr,1}.Speed_in_trial{ss,1}(round(mean(PertStart_i(:,1)))-30:round(mean(PertStart_i(:,2))))),... % pert period
                                      mean(ProjectData.Trial_data{rr,1}.Speed_in_trial{ss,1}(round(mean(PertStart_i(:,2)))-30:round(mean(PertStart_i(:,2))+60)))]; % post-pert period
            end
        end
        
        %% perturbation trials means RUNNING
        for sm = 1:size(TrialSpeeds,2)
            PertON_idx = find(ProjectData.Trial_data{rr,1}.PerturbationON==1);
            PertONRun_idx = find(TrialSpeeds(PertON_idx,sm)>3); % perturbation period
            RespSummay_total.PertRespMean_run{cellcount}(sm,1) = mean(max(PertRespMeanCell(PertONRun_idx,:,1),[],2)); % mean max response during perturbation
            RespSummay_total.PertRespMean_run{cellcount}(sm,2) = mean(max(PertRespMeanCell(PertONRun_idx,:,2),[],2)); % mean max response post perturbation
            RespSummay_total.PertRespMean_run{cellcount}(sm,3) = mean(mean(PertRespMeanCell(PertONRun_idx,:,3),2)); % mean response pre perturbation
            RespSummay_total.PertRespMean_run{cellcount}(sm,4) = mean(mean(PertRespMeanCell(PertONRun_idx,:,1),2)); % mean response during perturbation
            RespSummay_total.PertRespMean_run{cellcount}(sm,5) = mean(mean(PertRespMeanCell(PertONRun_idx,:,2),2)); % mean response post perturbation
            % perturbation trials means STILL
            PertONStill_idx = find(TrialSpeeds(PertON_idx,sm)<3); % perturbation period
            RespSummay_total.PertRespMean_still{cellcount}(sm,1) = mean(max(PertRespMeanCell(PertONStill_idx,:,1),[],2)); % mean max response during perturbation
            RespSummay_total.PertRespMean_still{cellcount}(sm,2) = mean(max(PertRespMeanCell(PertONStill_idx,:,2),[],2)); % mean max response post perturbation
            RespSummay_total.PertRespMean_still{cellcount}(sm,3) = mean(mean(PertRespMeanCell(PertONStill_idx,:,3),2)); % mean response pre perturbation
            RespSummay_total.PertRespMean_still{cellcount}(sm,4) = mean(mean(PertRespMeanCell(PertONStill_idx,:,1),2)); % mean response during perturbation
            RespSummay_total.PertRespMean_still{cellcount}(sm,5) = mean(mean(PertRespMeanCell(PertONStill_idx,:,2),2)); % mean response post perturbation
            
            % non-perturbation trials RUNNING
            PertOFFRun_idx = find(ProjectData.Trial_data{rr,1}.PerturbationON==0 & TrialSpeeds(:,sm)>3); % perturbation period
            RespSummay_total.PertRespMean_off_run{cellcount}(sm,1) = mean(max(SpikeCount_smooth(PertOFFRun_idx,round(mean(Pert_idx_start)):round(mean(Pert_idx_end))),[],2));
            RespSummay_total.PertRespMean_off_run{cellcount}(sm,2) = mean(max(SpikeCount_smooth(PertOFFRun_idx,round(mean(Pert_idx_end))+1:round(mean(Pert_idx_end))+1+60),[],2));
            RespSummay_total.PertRespMean_off_run{cellcount}(sm,3) = mean(mean(SpikeCount_smooth(PertOFFRun_idx,round(mean(Pert_idx_start))-60:round(mean(Pert_idx_start))),2));
            RespSummay_total.PertRespMean_off_run{cellcount}(sm,4) = mean(mean(SpikeCount_smooth(PertOFFRun_idx,round(mean(Pert_idx_start)):round(mean(Pert_idx_end))),2));
            RespSummay_total.PertRespMean_off_run{cellcount}(sm,5) = mean(mean(SpikeCount_smooth(PertOFFRun_idx,round(mean(Pert_idx_end))+1:round(mean(Pert_idx_end))+1+60),2));
            % non-perturbation trials STILL
            PertOFFStill_idx = find(ProjectData.Trial_data{rr,1}.PerturbationON==0 & TrialSpeeds(:,sm)<3); % perturbation period
            RespSummay_total.PertRespMean_off_still{cellcount}(sm,1) = mean(max(SpikeCount_smooth(PertOFFStill_idx,round(mean(Pert_idx_start)):round(mean(Pert_idx_end))),[],2));
            RespSummay_total.PertRespMean_off_still{cellcount}(sm,2) = mean(max(SpikeCount_smooth(PertOFFStill_idx,round(mean(Pert_idx_end))+1:round(mean(Pert_idx_end))+1+60),[],2));
            RespSummay_total.PertRespMean_off_still{cellcount}(sm,3) = mean(mean(SpikeCount_smooth(PertOFFStill_idx,round(mean(Pert_idx_start))-60:round(mean(Pert_idx_start))),2));
            RespSummay_total.PertRespMean_off_still{cellcount}(sm,4) = mean(mean(SpikeCount_smooth(PertOFFStill_idx,round(mean(Pert_idx_start)):round(mean(Pert_idx_end))),2));
            RespSummay_total.PertRespMean_off_still{cellcount}(sm,5) = mean(mean(SpikeCount_smooth(PertOFFStill_idx,round(mean(Pert_idx_end))+1:round(mean(Pert_idx_end))+1+60),2));
            % mean max response during perturbation period - still 
            % mean max response post perturbation period - still
            % mean response pre perturbation period - still
            % means response during perturbation period - still
            % means response during posturbation period - still
            
            % vis stim onset RUNNING
            Vis_Run = find(TrialSpeeds(:,sm)>3);
            RespSummay_total.VisRespMean_run{cellcount}(sm,1) = mean(mean(SpikeCount_smooth(Vis_Run,trialSide_samples+1:trialSide_samples+1+60),2)); % mean max response to visual stimulus onset
            RespSummay_total.VisRespMean_run{cellcount}(sm,2) = mean(mean(SpikeCount_smooth(Vis_Run,1:trialSide_samples),2)); % mean response to 1s before visual stimulus onset
            % vis stim onset STILL
            Vis_Still = find(TrialSpeeds(:,sm)<3);
            RespSummay_total.VisRespMean_still{cellcount}(sm,1) = mean(mean(SpikeCount_smooth(Vis_Still,trialSide_samples+1:trialSide_samples+1+60),2)); % mean max response to visual stimulus onset
            RespSummay_total.VisRespMean_still{cellcount}(sm,2) = mean(mean(SpikeCount_smooth(Vis_Still,1:trialSide_samples),2)); % mean response to 1s before visual stimulus onset
            
            RespSummay_total.Trial_types{cellcount}(sm,:) = [length(PertONRun_idx), length(PertONStill_idx), ...
                                                                length(PertOFFRun_idx), length(PertOFFStill_idx), ...
                                                                length(Vis_Run), length(Vis_Still)];
        end     
        
        cellcount = cellcount + 1;        
    end
%     n_PT_cond(rr,1) = sum(ProjectData.Trial_data{rr,1}.PerturbationON==1 & TrialSpeeds(:,2)>3);
%     n_PT_cond(rr,2) = sum(ProjectData.Trial_data{rr,1}.PerturbationON==1 & TrialSpeeds(:,2)<3);
%     n_PT_cond(rr,3) = size(ProjectData.Units_Info{rr,1},1);
end

% perturbation trials running (during perturbation period)
clear PertRespMean_on_run PertRespMean_on_still PostPertRespMean_on_run PostPertRespMean_on_still VisRespMean_run VisRespMean_still
c_er=1; c_es=1; c_et=1;
for h = 1:size(RespSummay_total,1)
   % select only sessions in which there are at least 10 running and 10 still trials
   if RespSummay_total.Trial_types{h}(2,1)>10 && RespSummay_total.Trial_types{h}(2,2)>10
       PertRespMean_on_run(c_er,:) = RespSummay_total.PertRespMean_run{h}(2,:); % 2nd row contains mean FR of pertON trials with speed in pert period >3cms
       PertRespMean_on_still(c_er,:) = RespSummay_total.PertRespMean_still{h}(2,:); % 2nd row contains mean FR of pertON trials with speed in pert period <3cms
       c_er = c_er+1;
       % see lines 78-82 for values' meaning
   end
   if RespSummay_total.Trial_types{h}(3,1)>10 && RespSummay_total.Trial_types{h}(3,2)>10
       PostPertRespMean_on_run(c_es,:) = RespSummay_total.PertRespMean_run{h}(3,:); % 2nd row contains mean FR of pertON trials with speed in pert period >3cms
       PostPertRespMean_on_still(c_es,:) = RespSummay_total.PertRespMean_still{h}(3,:); % 2nd row contains mean FR of pertON trials with speed in pert period <3cms
       c_es = c_es+1;
   end
   if RespSummay_total.Trial_types{h}(1,6)>10 && RespSummay_total.Trial_types{h}(1,6)>10
       VisRespMean_run(c_et,:) = RespSummay_total.VisRespMean_run{h}(1,:); % 1st row contains mean FR of all trials with speed in 1st 4 sec >3cms
       VisRespMean_still(c_et,:) = RespSummay_total.VisRespMean_still{h}(1,:); % 1st row contains mean FR of all trials with speed in 1st 4 sec >3cms
       c_et = c_et+1;
   end
end

clear ModI_p ModI_pp ModI_v
ModI_p(:,1) = (PertRespMean_on_run(:,4)-PertRespMean_on_run(:,3))./(PertRespMean_on_run(:,4)+PertRespMean_on_run(:,3)); % modulation index to perturbation
ModI_p(:,2) = (PertRespMean_on_still(:,5)-PertRespMean_on_still(:,4))./(PertRespMean_on_still(:,5)+PertRespMean_on_still(:,4)); % modulation index to post perturbation
ModI_pp(:,1) = (PostPertRespMean_on_run(:,5)-PostPertRespMean_on_run(:,4))./(PostPertRespMean_on_run(:,5)+PostPertRespMean_on_run(:,4)); % modulation index to perturbation period in non-pert trials
ModI_pp(:,2) = (PostPertRespMean_on_still(:,5)-PostPertRespMean_on_still(:,4))./(PostPertRespMean_on_still(:,5)+PostPertRespMean_on_still(:,4)); % modulation index to post perturbation period in non-pert trials
ModI_v(:,1) = (VisRespMean_run(:,1)-VisRespMean_run(:,2))./(VisRespMean_run(:,1)+VisRespMean_run(:,2)); % modulation index to post perturbation
ModI_v(:,2) = (VisRespMean_still(:,1)-VisRespMean_still(:,2))./(VisRespMean_still(:,1)+VisRespMean_still(:,2)); % modulation index to post perturbation

% compute hist
edges = 0:0.05:1;
Xaxis =0.05/2:0.05:0.05*20-0.05/2;
% Set up fittype and options.
ft = fittype( 'exp1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Normalize = 'on';
clear fitresult occ
for i = 1:size(ModI_p,2)
    [occ(:,i), edges] = histcounts(abs(ModI_p(:,i)),edges); % Pert period during pert trials
    [xData, yData] = prepareCurveData( Xaxis, occ(:,i)/sum(occ(:,i)) );
    % Fit model to data.
    [fitresult{i}, gof] = fit( xData, yData, ft, opts );
end

figure
set(gcf, 'Renderer', 'painters', 'Position',[-1500 500 1400 400]);
TitlesString = {'Perturbation period MI',...
                'Post-perturbation period MI',...
                'pert OFF Pert. MI',...
                'pert OFF post-pert. MI',...
                'Visual stimulus MI'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1,3,1)
% Pert period during pert trials
b1 = bar(Xaxis, occ(:,1)/sum(occ(:,1)),0.8,'FaceColor','r')
b1.FaceAlpha = 0.5;
hold on
f1 = plot(fitresult{1},'r')
f1.LineWidth = 1;
hold on
% Pert period during non-pert trials
b2 = bar(Xaxis, occ(:,2)/sum(occ(:,2)),0.5,'FaceColor','b')
b2.FaceAlpha = 0.5;
hold on
f2 = plot(fitresult{2},'b')
f2.LineWidth = 1;
xlabel('Modulation Index','FontSize',12); ylabel('Fraction of units','FontSize',12);
set(gca,'TickDir','out')
box off
title(TitlesString{1},'FontSize',15); ylim([0 max(max(occ(:,1)/sum(occ(:,1))),max(occ(:,2)/sum(occ(:,2))))+0.01])
ll = legend([b1 b2],{'Pert. ON - running', 'Pert. ON - still'},'FontSize',10)
ll.Color = 'none'; ll.EdgeColor = 'none'; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1,3,2)
clear fitresult occ
for i = 1:size(ModI_pp,2)
    [occ(:,i), edges] = histcounts(abs(ModI_p(:,i)),edges); % Pert period during pert trials
    [xData, yData] = prepareCurveData( Xaxis, occ(:,i)/sum(occ(:,i)) );
    % Fit model to data.
    [fitresult{i}, gof] = fit( xData, yData, ft, opts );
end
% Pert period during pert trials
b1 = bar(Xaxis, occ(:,1)/sum(occ(:,1)),0.8,'FaceColor','r')
b1.FaceAlpha = 0.5;
hold on
f1 = plot(fitresult{1},'r')
f1.LineWidth = 1;
hold on
% Pert period during non-pert trials
b2 = bar(Xaxis, occ(:,2)/sum(occ(:,2)),0.5,'FaceColor','b')
b2.FaceAlpha = 0.5;
hold on
f2 = plot(fitresult{2},'b')
f2.LineWidth = 1;
xlabel('Modulation Index','FontSize',12); ylabel('Fraction of units','FontSize',12);
set(gca,'TickDir','out')
box off
title(TitlesString{2},'FontSize',15); ylim([0 max(max(occ(:,1)/sum(occ(:,1))),max(occ(:,2)/sum(occ(:,2))))+0.01])
ll = legend([b1 b2],{'Post-pert - running', 'Post-pert - still'},'FontSize',10)
ll.Color = 'none'; ll.EdgeColor = 'none'; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1,3,3)
clear fitresult occ
for i = 1:size(ModI_v,2)
    [occ(:,i), edges] = histcounts(abs(ModI_v(:,i)),edges); % Pert period during pert trials
    [xData, yData] = prepareCurveData( Xaxis, occ(:,i)/sum(occ(:,i)) );
    % Fit model to data.
    [fitresult{i}, gof] = fit( xData, yData, ft, opts );
end
b1 = bar(Xaxis, occ(:,1)/sum(occ(:,1)),0.8,'FaceColor','r')
b1.FaceAlpha = 0.5;
hold on
f1 = plot(fitresult{1},'r')
f1.LineWidth = 1;
hold on
b2 = bar(Xaxis, occ(:,2)/sum(occ(:,2)),0.5,'FaceColor','b')
b2.FaceAlpha = 0.5;
hold on
f2 = plot(fitresult{2},'b')
f2.LineWidth = 1;
xlabel('Modulation Index','FontSize',12); ylabel('Fraction of units','FontSize',12);
set(gca,'TickDir','out')
box off
title(TitlesString{5},'FontSize',15);
ll = legend([b1 b2],{'Vis.Stim. Onset - running', 'Vis.Stim. Onset - still'},'FontSize',10)
ll.Color = 'none'; ll.EdgeColor = 'none'; 
ylim([0 max(max(occ(:,1)/sum(occ(:,1))),max(occ(:,2)/sum(occ(:,2))))+0.01]);

end



