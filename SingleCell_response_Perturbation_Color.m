%% Tomaso Muzzu - UCL - 06/02/2019

% function to plot basic PSTH to visual responses looking only at a
% specific angle for all cells


function SingleCell_response_Perturbation_Color(Array,SpikesTrial,ES,FileName,RecordingOI,Direction, RunningThreshold, MaxSpeed,StimulusParams,value2show)

trialSide_seconds = 1; % take 60 samples before and after trial
trialSide_samples = trialSide_seconds*round(length(ES.Time)/(max(ES.Time)-min(ES.Time))); % take 60 samples before and after trial
    
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
    
    FR(1,:) = sum(SpikeCount(Array.PerturbationON==0,:)) / sum(Array.PerturbationON==0);
    FR(2,:) = sum(SpikeCount(Array.PerturbationON==1,:)) / sum(Array.PerturbationON==1);
    
    if value2show == 0
        NormFR(i,:) = (FR(2,:)-FR(1,:))/max(max(FR));
    else
        NormFR(i,:) = (FR(2,:))/(max(FR(2,:)));
    end
end

% sort the cells by their response to perturbation 
SortingColumns = [mean(NormFR(:,StimulusParams{min(find(Array.PerturbationON==1))}(:,3)==0),2), ... % perturbation period
                  mean(NormFR(:,trialSide_seconds:trialSide_seconds+60),2),... % first second of trial
                  mean(NormFR(:,max(find(StimulusParams{min(find(Array.PerturbationON==1))}(:,3)==0))+1:...
                                max(find(StimulusParams{min(find(Array.PerturbationON==1))}(:,3)==0))+1+60),2)];   % second after perturbation
[B,I] = sortrows(SortingColumns,[-1 -2 -3]);

% make actual figure
TimeLine = linspace(StimulusParams{min(find(Array.PerturbationON==1))}(1,1)-trialSide_seconds,max(StimulusParams{min(find(Array.PerturbationON==1))}(:,1))-trialSide_seconds,length(FR(1,:)));
figure('Renderer', 'painters', 'Position', [10 10 1000 1000])
subplot(6,1,[2:6])
colormap(parula(256))
imagesc(TimeLine,linspace(1,size(NormFR,1),size(NormFR,1)), ...
    NormFR(I,:));
colorbar
xlabel('seconds'); ylabel('units');
if value2show == 0
    title('Diff response of pert. trials - no-pert. trials')
else
    title('Response during pert. trials')
end

% plot also contrast and TF
subplot(6,1,1)
h1 = plot(StimulusParams{min(find(Array.PerturbationON==1))}(:,1)-trialSide_seconds,StimulusParams{min(find(Array.PerturbationON==1))}(:,2),'k','LineWidth',2)
hold on
h2 = plot(StimulusParams{min(find(Array.PerturbationON==1))}(:,1)-trialSide_seconds,StimulusParams{min(find(Array.PerturbationON==1))}(:,3),'b','LineWidth',2)
hold on
plot([0 0],...
    [0 max(StimulusParams{min(find(Array.PerturbationON==1))}(:,3))],'--r','LineWidth',0.75);
hold on
plot([ES.ACInfo.trialStartsEnds(1,2)-ES.ACInfo.trialStartsEnds(1,1) ES.ACInfo.trialStartsEnds(1,2)-ES.ACInfo.trialStartsEnds(1,1)],...
    [0 max(StimulusParams{min(find(Array.PerturbationON==1))}(:,3))],'--r','LineWidth',0.75);
hold on
plot([ES.ACInfo.PertStartsEnds(end,1)-ES.ACInfo.trialStartsEnds(max(find(Array.PerturbationON)),1) ...
    ES.ACInfo.PertStartsEnds(end,2)-ES.ACInfo.trialStartsEnds(max(find(Array.PerturbationON)),1)],...
    [max(StimulusParams{min(find(Array.PerturbationON==1))}(:,3)) max(StimulusParams{min(find(Array.PerturbationON==1))}(:,3))],'--b','LineWidth',1);

xlim([-trialSide_seconds max(StimulusParams{min(find(Array.PerturbationON==1))}(:,1))-trialSide_seconds]);
ylim([0 max(StimulusParams{min(find(Array.PerturbationON==1))}(:,3))+0.1]);
box off
set(gca,'TickDir','out','XTickLabel',[])    
legend({'Constrast', 'TF'},'location','northwest')
colorbar
    

pattern = 'Matlab';
StartLetter = strfind(FileName,pattern);
nrTrialsPertDir = length(find(ES.PertTrial==1));
nrTrialsNPertDir = length(find(ES.PertTrial==0));
DetailCell = {[FileName(StartLetter+length(pattern)+1:StartLetter+length(pattern)+13) ', ' FileName(end-22:end-4)] , ...
    char(strcat(['All Cells, All Directions, n_p=' num2str(nrTrialsPertDir) ', n_n_p=' num2str(nrTrialsNPertDir)]))};
DetailCell;
a = axes;
t1 = title(DetailCell);
a.Visible = 'off'; % set(a,'Visible','off');
t1.Visible = 'on'; % set(t1,'Visible','on');

end  





