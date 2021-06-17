%% Tomaso Muzzu - UCL - 13/09/2019

% plot basic responses of units for visual perturbation experiments

function PlotResponses_AccelStim_2(ProjectData,SelectedResponses,AM_Param,AM_Speed,Units_Sel,Diff)

% AM_param
% conditions = 1 --> nr of recording
% conditions = 2 --> TF at stimulus onset
% conditions = 3 --> TF at simutlus offset
% conditions = 4 --> dTF shown
% conditions = 5 --> index of when Tf starts changing
% conditions = 6 --> index of when Tf finishes changing
AM_UOI = Units_Sel;

% find all the TF's at the start of the trials
StarTF_2dmask = squeeze(AM_Param(:,:,2));
start_TF_values = unique(StarTF_2dmask(~isnan(StarTF_2dmask(:))));
% separate the responses for the different acceleration values
dTF_2Dmask = squeeze(AM_Param(:,:,4));
dTF_values = unique(dTF_2Dmask(~isnan(dTF_2Dmask(:))));
% find trials in which dTF lasts 0.5s and 1s
dTF_dur_2D = (squeeze(AM_Param(:,:,6))-squeeze(AM_Param(:,:,5)))/60; % 60 is the sampling frequency
dTF_long_idx = dTF_dur_2D>0.8;
% find three trials for each condition for which there is dTF>0; dTF=0;
% dTF<0
clear TF_10_i TF_05_i I_10 I_05
for i = 1:length(dTF_values)
    TF_10_i(i) = min(find((AM_Param(:,:,4).*dTF_long_idx)==dTF_values(i) & (AM_Param(:,:,2).*dTF_long_idx)==6));
    TF_05_i(i) = min(find((AM_Param(:,:,4).*(~dTF_long_idx))==dTF_values(i) & (AM_Param(:,:,2).*(~dTF_long_idx))==6));
end
%convert to 2D indexes
ReferSize = [size(AM_Param,1),size(AM_Param,2)];
[I_10(:,1),I_10(:,2)] = ind2sub(ReferSize,TF_10_i); % from linear to 2D indexes
[I_05(:,1),I_05(:,2)] = ind2sub(ReferSize,TF_05_i); % from linear to 2D indexes


Recs = unique(AM_Param(1,:,1)); clear UnitPerRec
for k = 1:length(Recs)
    UnitPerRec(k) = max(find(AM_Param(1,:,1)==Recs(k)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% determine number of units for each condition:
TotalUnits_05s = find(~isnan(squeeze(SelectedResponses{1,1}(I_05(1,1),:,1))));

% plot the responses for each unit distinguishing the dTF values
% for i = 1 : length(TotalUnits_05s)
%    figure('Renderer', 'painters', 'Position', [50 100 1800 400]);
%    subplot(1,4,1) % plot responses for different pos., neg., and zero acceleration values
%    for dtf = 1:length(dTF_values)
%        plot(linspace(-1,4,299), smooth(nanmean(squeeze(SelectedResponses{dtf,1}(:,TotalUnits_05s(i),:)),1),5),'Color',([1 0 0])*(dtf/length(dTF_values)-1/length(dTF_values)) );
%        FR(dtf) = max(nanmean(squeeze(SelectedResponses{dtf,1}(:,TotalUnits_05s(i),:)),1));
%        hold on
%    end
%    plotTrialMarks(max(FR))
%    subplot(1,4,1) 
%    title(['Rec.: ' num2str(AM_Param(1,TotalUnits_05s(i),1)) ...
%           ', K__id: ' num2str(ProjectData.Units_Info{AM_Param(1,TotalUnits_05s(i),1),1}.k_ID{TotalUnits_05s(i)-UnitPerRec(max(AM_Param(1,TotalUnits_05s(i),1)-1,1)),1} )  ...
%           ', all start TFs']);
%    legend({num2str(dTF_values)});
%    %%%%%%%%%%%%%%%%%%%%%%%
%    for stf = 1:length(start_TF_values)
%        subplot(1,4,1+stf) % plot responses for the seven acceleration values.; Starting TF = 1.5
%        for dtf = 1:length(dTF_values)
%            plot(linspace(-1,4,299), smooth(nanmean(squeeze(SelectedResponses{dtf,1}(find(AM_Param(:,TotalUnits_05s(i),2)==start_TF_values(stf)),TotalUnits_05s(i),:)),1),5),...
%                'Color',([1 0 0])*(dtf/length(dTF_values)-1/length(dTF_values)) )
%            FR(dtf) = max(nanmean(squeeze(SelectedResponses{dtf,1}(find(AM_Param(:,TotalUnits_05s(i),2)==start_TF_values(stf)),TotalUnits_05s(i),:)),1));
%            hold on
%        end
%        plotTrialMarks(max(FR))
%        subplot(1,4,1+stf)
%        title(['Rec.: ' num2str(AM_Param(1,TotalUnits_05s(i),1)) ...
%            ', K__id: ' num2str(ProjectData.Units_Info{AM_Param(1,TotalUnits_05s(i),1),1}.k_ID{TotalUnits_05s(i)-UnitPerRec(max(AM_Param(1,TotalUnits_05s(i),1)-1,1)),1} )  ...
%            ', startTF: ' num2str(start_TF_values(stf)) 'Hz']);
%        legend({num2str(dTF_values)});
%    end
%    %%%%%%%%%%%%%%%%%%%%%%%
%    
%    saveas(gcf, ['X:\DATA\PROJECTS\VisPerturbation\Figures\AccelResponse' filesep ...
%                 'AccelResponse_' ...
%                 'Rec#' num2str(AM_Param(1,TotalUnits_05s(i),1)) ...
%                 '_Kid#' num2str(ProjectData.Units_Info{AM_Param(1,TotalUnits_05s(i),1),1}.k_ID{TotalUnits_05s(i)-UnitPerRec(max(AM_Param(1,TotalUnits_05s(i),1)-1,1)),1})...
%                 '.png']);
%     close all
%     
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% determine number of units for each condition:
TotalUnits_1s = find(~isnan(squeeze(SelectedResponses{1,2}(I_10(1,1),:,1))));
    % plot the responses for each unit distinguishing the dTF values
for i = 1 : length(TotalUnits_1s)
   figure('Renderer', 'painters', 'Position', [50 100 1800 400]);
   subplot(1,4,1) % plot responses for different pos., neg., and zero acceleration values
   for dtf = 1:length(dTF_values)
       plot(linspace(-1,4,299), smooth(nanmean(squeeze(SelectedResponses{dtf,2}(:,TotalUnits_1s(i),:)),1),5),'Color',([1 0 0])*(dtf/length(dTF_values)-1/length(dTF_values)) );
       FR(dtf) = max(nanmean(squeeze(SelectedResponses{dtf,2}(:,TotalUnits_1s(i),:)),1));
       hold on
   end
   plotTrialMarks(max(max(FR),1))
   subplot(1,4,1) 
   title(['Rec.: ' num2str(AM_Param(1,TotalUnits_1s(i),1)) ...
          ', K__id: ' num2str(ProjectData.Units_Info{AM_Param(1,TotalUnits_1s(i),1),1}.k_ID{max(abs(TotalUnits_1s(i)-UnitPerRec(max(AM_Param(1,TotalUnits_1s(i),1)-1,1))),1),1} )  ...
          ', all start TFs']);
   legend({num2str(dTF_values)});
   %%%%%%%%%%%%%%%%%%%%%%%
   for stf = 1:length(start_TF_values)
       subplot(1,4,1+stf) % plot responses for the seven acceleration values.; Starting TF = 1.5
       for dtf = 1:length(dTF_values)
           plot(linspace(-1,4,299), smooth(nanmean(squeeze(SelectedResponses{dtf,2}(find(AM_Param(:,TotalUnits_1s(i),2)==start_TF_values(stf)),TotalUnits_1s(i),:)),1),5),...
               'Color',([1 0 0])*(dtf/length(dTF_values)-1/length(dTF_values)) )
           FR(dtf) = max(nanmean(squeeze(SelectedResponses{dtf,2}(find(AM_Param(:,TotalUnits_1s(i),2)==start_TF_values(stf)),TotalUnits_1s(i),:)),1));
           hold on
       end
       plotTrialMarks(max(max(FR),1))
       subplot(1,4,1+stf)
       title(['Rec.: ' num2str(AM_Param(1,TotalUnits_1s(i),1)) ...
           ', K__id: ' num2str(ProjectData.Units_Info{AM_Param(1,TotalUnits_1s(i),1),1}.k_ID{max(abs(TotalUnits_1s(i)-UnitPerRec(max(AM_Param(1,TotalUnits_1s(i),1)-1,1))),1),1} )  ...
           ', startTF: ' num2str(start_TF_values(stf)) 'Hz']);
       legend({num2str(dTF_values)});
   end
   %%%%%%%%%%%%%%%%%%%%%%%
   
   saveas(gcf, ['X:\DATA\PROJECTS\VisPerturbation\Figures\AccelResponse' filesep ...
                'AccelResponse_' ...
                'Rec#' num2str(AM_Param(1,TotalUnits_1s(i),1)) ...
                '_Kid#' num2str(ProjectData.Units_Info{AM_Param(1,TotalUnits_1s(i),1),1}.k_ID{max(abs(TotalUnits_1s(i)-UnitPerRec(max(AM_Param(1,TotalUnits_1s(i),1)-1,1))),1),1})...
                '.png']);
    close all
    
   
end
    


    

end


function plotTrialMarks(FR)
        hold on
        plot([0.5 0.5],[0 max(FR)*1.1],'g--','LineWidth',2); plot([1.5 1.5],[0 max(FR)*1.1],'g--','LineWidth',2); % dTF period
        plot([0 0],[0 max(FR)*1.1],'k:','LineWidth',4); plot([3 3],[0 max(FR)*1.1],'k:','LineWidth',4); % trial period
        set(gca,'TickDir','out','XTickLabel',[],'XTick',[0 0.5 1 3])
        box off
        ylabel('Hz');
        ylim([0 max(FR)*1.1]); xlim([-0.5 3.5]);
        h = text(-0.1,2, 'trial start'); set(h,'Rotation',90);
        h = text(3.1,2, 'trial end'); set(h,'Rotation',90);
        h = text(0.45,2, 'dTF start'); set(h,'Rotation',90);
        h = text(1.05,2, 'dTF stop'); set(h,'Rotation',90);
end