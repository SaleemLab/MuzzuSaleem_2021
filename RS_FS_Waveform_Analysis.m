%% Tomaso Muzzu - UCL - 13/05/2021 - Saleem Lab

% regular vs fast spiking neurons
% for loop that scans through all the recordings and save sequentially the
% data
%%
%% Functions to plot Figure 2 - direction tuning
if ~exist('ProjectData','var')
    [ProjectData AM_UnitResponses AM_Param AM_Speed AM_UOI SelectedResponses AM_UnitResponses_smooth AM_EyeTracking] = LoadDataALL;
end
%% direction
if size(AM_UOI,2)==1
    SelectedCells = AM_UOI;
else
    SelectedCells = AM_UOI(:,1) & AM_UOI(:,2);
end
% select trials of interest and control trials as well
Trials_PertON = AM_Param(:,:,3)==1; % find indexes where perturbation is on
Trials_PertOFF= AM_Param(:,:,3)==0; % find indexes where perturbation is off

BonvisionFR = 60; %Hz
trialSide_samples = 60;
trialSide_seconds = 1;

%%
p_value = 95;

%% first 7 animals
thres = 95;
if  size(ProjectData,1)==37
    CTRL_exp = 0;
    Animal_1st_idx = [1 5 7 12 15 24 31];
    if ~exist('PertResp_units','var')
        % select only perturbation responsive units
        load('AUC_shuffled.mat')
        Sh_responses = AUC_shuffled(:,2:end);
        p_pert_th = prctile(Sh_responses(:),thres);
        PertResp_units = (AUC_shuffled(:,1)>p_pert_th);
        % select only pos. modulated perturbation responsive units
        load('DM_pert_shuffled.mat')
        DM = DM_sh(:,1);
        DM_sign_i(:,1) = DM>0;
        DM_sign_i(:,2) = DM<=0;
        % select only pos. modulated perturbation responsive units
        PertRespUnits_pos = PertResp_units & DM_sign_i(:,1);
        PertRespUnits_neg = PertResp_units & DM_sign_i(:,2);
    end
elseif size(ProjectData,1)==10
    CTRL_exp = 1;
    % naive animals
    Animal_1st_idx = [1 4 7];
    % select only perturbation responsive units
    load('AUC_shuffled_CTRL_1.mat')
    Sh_responses = AUC_shuffled(:,2:end);
    p_pert_th = prctile(Sh_responses(:),thres);
    PertResp_units = (AUC_shuffled(:,1)>p_pert_th);
    % select only pos. modulated perturbation responsive units
    load('DM_pert_shuffled_CTRL.mat')
    DM = DM_sh(:,1);
    DM_sign_i(:,1) = DM>0;
    DM_sign_i(:,2) = DM<=0;
    % select only pos. modulated perturbation responsive units
    PertRespUnits_pos = PertResp_units & DM_sign_i(:,1);
    PertRespUnits_neg = PertResp_units & DM_sign_i(:,2);
else
    % select only perturbation responsive units
    load('AUC_shuffled_CTRL_1.mat')
    AUC_shuffled_CTRL = AUC_shuffled;
    AUC_shuffled_pp_CTRL = AUC_shuffled_pp;
    load('AUC_shuffled.mat')
    AUC_shuffled = cat(1,AUC_shuffled,AUC_shuffled_CTRL);
    AUC_shuffled_pp = cat(1,AUC_shuffled_pp,AUC_shuffled_pp_CTRL);
    Sh_responses = AUC_shuffled(:,2:end);
    p_pert_th = prctile(Sh_responses(:),thres);
    PertResp_units = (AUC_shuffled(:,1)>p_pert_th);
    % select only pos. modulated perturbation responsive units
    load('DM_ALL.mat')
    DM_sign_i(:,1) = DM>0;
    DM_sign_i(:,2) = DM<=0;
    % select only pos. modulated perturbation responsive units
    PertRespUnits_pos = PertResp_units & DM_sign_i(:,1);
    PertRespUnits_neg = PertResp_units & DM_sign_i(:,2);
end


% % name of dat file
% ProjectData.Session_data{1,1}.MetaData{1,1}.FileName
% 
% % cluster ID
% ProjectData.Units_Info{1,1}.k_ID

% INPUT
% gwfparams.dataDir = '/path/to/data/';    % KiloSort/Phy output folder
% gwfparams.fileName = 'data.dat';         % .dat file containing the raw 
% gwfparams.dataType = 'int16';            % Data type of .dat file (this should be BP filtered)
% gwfparams.nCh = 32;                      % Number of channels that were streamed to disk in .dat file
% gwfparams.wfWin = [-40 41];              % Number of samples before and after spiketime to include in waveform
% gwfparams.nWf = 2000;                    % Number of waveforms per unit to pull out
% gwfparams.spikeTimes =    [2,3,5,7,8,9]; % Vector of cluster spike times (in samples) same length as .spikeClusters
% gwfparams.spikeClusters = [1,2,1,1,1,2]; % Vector of cluster IDs (Phy nomenclature)   same length as .spikeTimes


%  for rec = 4:size(ProjectData,1)
%     clear SpikeTimes_def SpikeCluster_def
%     Temp_i_r = ProjectData.Session_data{rec,1}.RecordingOI{1};
%     StartTime_samples = ProjectData.Session_data{rec,1}.MetaData{1,1}.lims(Temp_i_r-1);
%     for j = 1:size(ProjectData.Units_Info{rec,1},1)
%         SpikeTimes_adj_temp = ProjectData.Units_Info{rec,1}.Spiketimes{j,1}*30000; 
%         SpikeCluster_temp = single(ProjectData.Units_Info{rec,1}.k_ID{j,1}) * ones(size(SpikeTimes_adj_temp,1),1);
%         if j == 1
%             SpikeTimes_def = SpikeTimes_adj_temp;
%             SpikeCluster_def = SpikeCluster_temp;
%         else
%             SpikeTimes_def = [SpikeTimes_def; SpikeTimes_adj_temp];
%             SpikeCluster_def = [SpikeCluster_def; SpikeCluster_temp]; 
%         end
%     end
%                    
%     pos_sep = strfind(ProjectData.Session_data{rec,1}.MetaData{1,1}.FileName,filesep);
%     if strfind(ProjectData.Session_data{rec,1}.MetaData{1,1}.FileName,'kilosort')
%         gwfparams.dataDir = ProjectData.Session_data{rec,1}.MetaData{1,1}.FileName(1:pos_sep(end));
%     else
%         gwfparams.dataDir = [ProjectData.Session_data{rec,1}.MetaData{1,1}.FileName(1:pos_sep(end)) 'kilosort\'];
%     end
%     gwfparams.fileName = ProjectData.Session_data{rec,1}.MetaData{1,1}.FileName(pos_sep(end)+1:end);
%     gwfparams.dataType = 'int16';
%     gwfparams.nCh = 32;
%     gwfparams.wfWin = [-60 61];              % -60 samples to +60 (30 samples = 1 ms)
%     gwfparams.nWf = 500;                    % Number of waveforms per unit to pull out
%     gwfparams.spikeTimes = SpikeTimes_def; % Vector of cluster spike times (in samples) same length as .spikeClusters
%     gwfparams.spikeClusters = SpikeCluster_def;  % Vector of cluster IDs (Phy nomenclature)   same length as .spikeTimes
%     
%     clear SpikeTimes_def SpikeCluster_def
%     
%     wf = getWaveForms(gwfparams);
%     clear gwfparams
%     rec
%     % save file every rec individually
%     save(['Waveforms_' num2str(rec) '_f.mat'], 'wf' ,'-v7.3')
%     clear wf
% end
% 


% OUTPUT
% wf.unitIDs                               % [nClu,1]            List of cluster IDs; defines order used in all wf.* variables
% wf.spikeTimeKeeps                        % [nClu,nWf]          Which spike times were used for the waveforms
% wf.waveForms                             % [nClu,nWf,nCh,nSWf] Individual waveforms
% wf.waveFormsMean                         % [nClu,nCh,nSWf]     Average of all waveforms (per channel)
%                                          % nClu: number of different clusters in .spikeClusters
%                                          % nSWf: number of samples per waveform

% for i = 1:47
%     load(['Waveforms_' num2str(i) '_f.mat']);
%     if i == 1
%         wf_all = wf;
%         clear wf;
%     else
%         wf_all = [wf_all wf];
%         clear wf;
%     end
%     i
% end
% 
% save('Waveforms_all_f.mat','wf_all','-v7.3')
% 
% wf = load('Waveforms_all_f.mat')
% wf = wf_light;
% save('Waveforms_all_f_means.mat','wf','-v7.3')

%% Aim: define behaviour of FS and RS neurons
% 1) choose the best units according to the selection criteria during the main
% stimulus (remember the correlation with the kluster ID)
% 2) find the best channel to compute the best waveform mean
% 3) compute the duration of the waveform & cluster two groups
% 4) find a correlation between the responsiveness and the RS/FS feature
AM_Param(1,AM_UOI,1)
exp = 40;
wf_trial = wf(exp);
AM_UOI(AM_Param(1,:,1)==exp)

for i = 1 1:size(wf_trial.waveFormsMean,1)
    figure
    for j = 1:size(wf_trial.waveFormsMean,2)
        subplot(2,16,j)
        plot((squeeze(wf_trial.waveFormsMean(i,j,:))))
        ylim([min(min((squeeze(wf_trial.waveFormsMean(i,:,:))')))...
              max(max((squeeze(wf_trial.waveFormsMean(i,:,:))')))]);
        title(num2str(j))
    end
    [v ind] = max(max(abs(squeeze(wf_trial.waveFormsMean(i,:,30:90))')));
    suptitle(['Channel: ' num2str(ind)])
end

figure
subplot(2,1,1)
volt_signal = (squeeze(wf_trial.waveFormsMean(i,ind,:)));
plot(linspace(-2,2,length(volt_signal)),volt_signal);
ylabel('dv/dt');
subplot(2,1,2)
MeanWaveforms_dt = diff(squeeze(wf_trial.waveFormsMean(i,ind,:)));
plot(linspace(-2,2,length(MeanWaveforms_dt)),MeanWaveforms_dt);
ylabel('v');
xlabel('time (ms)');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% find mean waveform
fs = 30000; % samples per second
% [v ind_t(1)] = max(First_der);
% [v ind_t(2)] = min(First_der);
% Spike_width = abs(diff(ind_t))*(1/fs)*1000;

u_i = 1; clear MeanWaveforms Spike_width k_ID MeanWaveforms_dt
for rec = 1:size(wf,2)
    wf_trial = wf(rec);
    for i = 1:size(wf_trial.waveFormsMean,1)
        % find the best channel
        [v ind] = max(max(abs(squeeze(wf_trial.waveFormsMean(i,:,30:90))')));
        MeanWaveforms_temp = (squeeze(wf_trial.waveFormsMean(i,ind,:))');
        MeanWaveforms(u_i,:) = interp1(linspace(-2,2,122),MeanWaveforms_temp,linspace(-2,2,400));
        % find the trhough-to-peak width (duration) 
        [v ind_t(1)] = min(MeanWaveforms(u_i,:));
        [v ind_t(2)] = max(MeanWaveforms(u_i,floor(size(MeanWaveforms,2)/2)+1:end));
        ind_t(2) = ind_t(2)+floor(size(MeanWaveforms,2)/2);
%         % compute first derivative
%         MeanWaveforms_dt(u_i,:) = diff(MeanWaveforms(u_i,:));
%         % find the width (duration) at half-maximum
%         [v ind_t(1)] = min(MeanWaveforms_dt(u_i,:));
%         [v ind_t(2)] = max(MeanWaveforms_dt(u_i,floor(size(MeanWaveforms_dt,2)/2)+1:end));
%         ind_t(2) = ind_t(2)+floor(size(MeanWaveforms_dt,2)/2);
        Spike_width(u_i) = abs(diff(ind_t))*(1/30)*(122/400); % in ms
        k_ID(u_i) = wf_trial.unitIDs(i);
        u_i = u_i + 1;
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% compute mean FR
% initialise gaussian filter parameters for firing rate
BonvisionFR = 60;
GaussFilter_W = 0.3; % seconds
Width = round(GaussFilter_W*BonvisionFR);
Sigma = Width/3; % standard deviation in number of samples (converted from time in seconds)
x_g = linspace(-Width/2, Width/2, Width);
gaussFilter = exp(-x_g.^2/(2*Sigma^2));
gaussFilter_ = gaussFilter / sum (gaussFilter); % normalize
ISI_edges_int = [0.002 0.1];
ISI_edges_low = 0:ISI_edges_int(1):1; 
ISI_edges = [ISI_edges_low 1:ISI_edges_int(2):10];
n_n = 1; clear Response Mean_FR ISI_n
for rec = 1:size(ProjectData,1)
    %%%%%%%%%%%%%%%%%% 3) prepare data
    ES = ProjectData.Session_data{rec,1};
    % prepare FR's as array (# of neurons) X (# of independent variables)
    clear Response
    for neuron = 1:size(ES.SpikeInfo{1,1},2)-1
        Time_edges = ES.Time{1,1};
        [SpikeCount_temp, Edges] = histcounts(ES.SpikeInfo{1,1}{1,neuron+1},Time_edges);
        SpikeCount = SpikeCount_temp*BonvisionFR; % to get Hz
        UnitFR = [0 conv(SpikeCount, gaussFilter_, 'same')];
        Response(neuron,:) = UnitFR;
        
        Mean_FR(n_n) = nanmean( UnitFR(ES.Contrast{1,1}==0) );         % select only periods of gray screen
        
        ISI_n_temp = histcounts(diff(ES.SpikeInfo{1,1}{1,neuron+1}),ISI_edges);
        ISI_n(n_n,:) = ISI_n_temp/length(ISI_n_temp);
        n_n = n_n + 1;
    end
end


%% select the units of interest
% select units that are active during the main stimulus
Spike_Dur = Spike_width(AM_UOI);

PertRespUnits_pos = PertResp_units & DM_sign_i(:,1);
PertRespUnits_neg = PertResp_units & DM_sign_i(:,2);

[c_idx, waveforms, explained] = classify_ns_bs_cells_T(MeanWaveforms(:,:)');
cl_id = c_idx(AM_UOI);

%% look at waveform duration VS mean FR
Mean_FR_ = Mean_FR(AM_UOI);
figure
plot(Mean_FR_(cl_id==1),Spike_Dur(cl_id==1),'b.')
hold on
plot(Mean_FR_(cl_id==3),Spike_Dur(cl_id==3),'r.')

edges = linspace(0,1.4,30); clear Distr;
Distr(1,:) = histcounts(Spike_Dur(cl_id==1),edges);
Distr(2,:) = histcounts(Spike_Dur(cl_id==3),edges);
figure
x_axis = edges(1:end-1)+mean(diff(edges));
stairs(x_axis,Distr(1,:),'k')
hold on
stairs(x_axis,Distr(2,:),'r')
set(gca,'box','off','TickDir','out')
xlabel('Trough to peak (ms)');
ylabel('Units');
legend({'Regular spiking', 'fast spiking'})

% figure
% boxplot(Spike_Dur(cl_id==1 & PertRespUnits_pos),'.')
% hold on
% boxplot(2,Spike_Dur(cl_id==2 & PertRespUnits_pos),'.')

% figure to plot for FS and RS spiking
Data2Scatter = padcat(Spike_Dur(cl_id==1 & PertRespUnits_pos), Spike_Dur(cl_id==2 & PertRespUnits_pos), Spike_Dur(cl_id==3 & PertRespUnits_pos));
u_ns = [length(Spike_Dur(cl_id==1 & PertRespUnits_pos)), length(Spike_Dur(cl_id==2 & PertRespUnits_pos)), length(Spike_Dur(cl_id==3 & PertRespUnits_pos))];
figure
UnivarScatter(Data2Scatter','Label',{['Cluster 1 n=' num2str(u_ns(1))], ...
                                     ['Cluster 2 n=' num2str(u_ns(2))], ...
                                     ['Cluster 3 n=' num2str(u_ns(3))],},'MarkerFaceColor',[ 0.5 0.5 0.5])
xlabel('Perturbation responsive units');
ylabel('Trough to peak duration (ms)');
set(gca,'TickDir','out','box','off')
title('Pos. mod. units')

% negatively modulated
Data2Scatter = padcat(Spike_Dur(cl_id==1 & PertRespUnits_neg), Spike_Dur(cl_id==2 & PertRespUnits_neg), Spike_Dur(cl_id==3 & PertRespUnits_neg));
u_ns = [length(Spike_Dur(cl_id==1 & PertRespUnits_neg)), length(Spike_Dur(cl_id==2 & PertRespUnits_neg)), length(Spike_Dur(cl_id==3 & PertRespUnits_neg))];
figure
UnivarScatter(Data2Scatter','Label',{['Cluster 1 n=' num2str(u_ns(1))], ...
                                     ['Cluster 2 n=' num2str(u_ns(2))], ...
                                     ['Cluster 3 n=' num2str(u_ns(3))],},'MarkerFaceColor',[ 0.5 0.5 0.5])
xlabel('Perturbation responsive units');
ylabel('Trough to peak duration (ms)');
set(gca,'TickDir','out','box','off')
title('Neg. mod. units')


% figure to plot for FS and RS spiking
Data2Scatter = padcat(Spike_Dur(cl_id==1), Spike_Dur(cl_id==2), Spike_Dur(cl_id==3));
u_ns = [length(Spike_Dur(cl_id==1)), length(Spike_Dur(cl_id==2)), length(Spike_Dur(cl_id==3))];
Colors=ColorCoder(3);
figure
UnivarScatter(Data2Scatter','Label',{['Cluster 1 n=' num2str(u_ns(1))], ...
                                     ['Cluster 2 n=' num2str(u_ns(2))], ...
                                     ['Cluster 3 n=' num2str(u_ns(3))]},...
                                     'MarkerFaceColor',Colors,'SEMColor',Colors/1.5,'StdColor',Colors/2);
hold on
xlabel('All units');
ylabel('Spike width at half max (ms)');
set(gca,'TickDir','out','box','off')
title('All units')

% how many units of each group are perturbation responsive (positive or negative)?
Cluster1 = 1;
Cluster2 = 2;
Cluster3 = 3;
y = [length(Spike_Dur(cl_id == Cluster1 & PertRespUnits_pos)), ...
     length(Spike_Dur(cl_id == Cluster1 & PertRespUnits_neg)), ...
     length(Spike_Dur(cl_id == Cluster1 & PertRespUnits_pos==0 & PertRespUnits_neg==0)); ...
     length(Spike_Dur(cl_id == Cluster2 & PertRespUnits_pos)), ...
     length(Spike_Dur(cl_id == Cluster2 & PertRespUnits_neg)), ...
     length(Spike_Dur(cl_id == Cluster2 & PertRespUnits_pos==0 & PertRespUnits_neg==0));...
     length(Spike_Dur(cl_id == Cluster3 & PertRespUnits_pos)), ...
     length(Spike_Dur(cl_id == Cluster3 & PertRespUnits_neg)), ...
     length(Spike_Dur(cl_id == Cluster3 & PertRespUnits_pos==0 & PertRespUnits_neg==0))];
ClustersTotal = sum(y');
y = y./ClustersTotal'*100;

figure
X = categorical({['p. fast spiking \newline n = ' num2str(ClustersTotal(1))] , ...
                 ['p. regular spiking \newline n = ' num2str(ClustersTotal(2))] , ...
                 ['unclassified \newline n = ' num2str(ClustersTotal(3))] });
b = bar(X,y,'stacked')
for i = 1:3
    subplot(1,3,i)
    pie(y(i,:))
end
legend({'Pert. resp. pos', 'Pert. resp. neg', 'Rest of units'});
set(gca,'TickDir','out','box','off')

xtips1 = [b(1).XEndPoints b(2).XEndPoints b(3).XEndPoints];
ytips1 = [b(1).YEndPoints b(2).YEndPoints b(3).YEndPoints];
labels1 = string([b(1).YData b(2).YData b(3).YData]);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','top','Color',[1 1 1])
ylabel('% units')

colorSet = [];
myColors = [1 0 0; 0 0 1; 0.5 0.5 0.5];
for i = 1:3
    %colorSet = [colorSet myColors];
    b(i).FaceColor = myColors(i,:);
    %b(i).CData = myColors;
end



%% plot mean waveforms for the two main groups to show difference
Cluster1 = 1;
Cluster2 = 3;
MeanWaveforms_OI = MeanWaveforms(AM_UOI,:);
MeanWaveforms_OI = ((MeanWaveforms_OI'-min(MeanWaveforms_OI'))./(max(MeanWaveforms_OI')-min(MeanWaveforms_OI')))';
figure
set(0, 'DefaultFigureRenderer', 'painters')
subplot(1,2,1)
plot(linspace(-2,2,size(MeanWaveforms_OI,2)), MeanWaveforms_OI(cl_id==Cluster1,:),'Color',[0.5 0.5 0.5]);
hold on
plot(linspace(-2,2,size(MeanWaveforms_OI,2)), mean(MeanWaveforms_OI(cl_id==Cluster1,:)),'k','LineWidth',3);
legend({'all waveforms','mean waveform'})
xlabel('Time (ms)');
ylabel('A.U.')
set(gca,'TickDir','out')
title(['pRS neurons, mean dur = ' num2str(mean(Spike_Dur(cl_id==Cluster1))) 'ms'])
box off
set(gca,'TickDir','out')
subplot(1,2,2)
plot(linspace(-2,2,size(MeanWaveforms_OI,2)), MeanWaveforms_OI(cl_id==Cluster2,:),'Color',[0.5 0.5 0.5]);
hold on
plot(linspace(-2,2,size(MeanWaveforms_OI,2)), mean(MeanWaveforms_OI(cl_id==Cluster2,:)),'k','LineWidth',3);
legend({'all waveforms','mean waveform'})
xlabel('Time (ms)');
ylabel('A.U.')
set(gca,'TickDir','out')
title(['pFS neurons, mean dur = ' num2str(mean(Spike_Dur(cl_id==Cluster2))) 'ms'])
box off
set(gca,'TickDir','out')

%% other stuff
mean(Spike_Dur(cl_id==1))
median(Spike_Dur(cl_id==1))
std(Spike_Dur(cl_id==1))

mean(Spike_Dur(cl_id==2))
median(Spike_Dur(cl_id==2))
std(Spike_Dur(cl_id==2))

figure
[N, edges_] = histcounts(Spike_Dur,0:0.04:2,'Normalization','pdf') ;
plot(linspace(edges_(1)+mean(diff(edges_)),edges_(end)-mean(diff(edges_)),length(N)),N,'.');


            

sum(PertRespUnits_pos)


figure
for p = 1:25
    subplot(5,5,p)
    rn = randi(size(MeanWaveforms,1));
    plot(linspace(-2,2,400),MeanWaveforms(rn,:))
%     hold on
%     plot(linspace(-2,2,399),MeanWaveforms_dt(rn,:))
    title(['cl=' num2str(c_idx(rn)) ', ' num2str(Spike_width(rn)) 'ms'])
end















