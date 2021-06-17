%% Tomaso Muzzu - UCL - 30/09/2019

% plot basic responses of units for visual perturbation experiments

function PlotResponses_AccelStim_1(ProjectData,SelectedResponses,AM_Param,AM_Speed,Units_Sel,Diff)

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
dTF_dur = [0.5 1];
for fig = 1:size(SelectedResponses,2)
    figure('Renderer', 'painters', 'Position', [10 10 1800 900]);
    
    for sf = 1:length(start_TF_values)
        clear TF_10_i TF_05_i I_10 I_05
        for i = 1:length(dTF_values)
            TF_10_i(i) = min(find((AM_Param(:,:,4).*dTF_long_idx)==dTF_values(i) & (AM_Param(:,:,2).*dTF_long_idx)==start_TF_values(sf)));
            TF_05_i(i) = min(find((AM_Param(:,:,4).*(~dTF_long_idx))==dTF_values(i) & (AM_Param(:,:,2).*(~dTF_long_idx))==start_TF_values(sf)));
        end
        %convert to 2D indexes
        ReferSize = [size(AM_Param,1),size(AM_Param,2)];
        [I_10(:,1),I_10(:,2)] = ind2sub(ReferSize,TF_10_i); % from linear to 2D indexes
        [I_05(:,1),I_05(:,2)] = ind2sub(ReferSize,TF_05_i); % from linear to 2D indexes
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %% plot stimulus dynamics feature, left is 0.5s and right is 1s
        subplot(10,3,[sf sf+3]); % show trials with dTF duration of 0.5s
        title('sessions with dTF duration 0.5s');
        for i = 1:size(I_05,1)
            plot(ProjectData.Trial_data{AM_Param(I_05(i,1),I_05(i,2),1),1}.time_in_trial{I_05(i,1),1}, ...
                ProjectData.Trial_data{AM_Param(I_05(i,1),I_05(i,2),1),1}.TF_in_trial{I_05(i,1),1},'Color',[0.5 0.5 0.5])
            hold on
        end
        plot([0.5 0.5],[0 9.5],'g--','LineWidth',2); plot([0.5+dTF_dur(fig) 0.5+dTF_dur(fig)],[0 9.5],'g--','LineWidth',2); % dTF period
        plot([0 0],[0 9.5],'k:','LineWidth',4); plot([3 3],[0 9.5],'k:','LineWidth',4); % trial period
        set(gca,'TickDir','out','XTickLabel',[],'XTick',[0 0.5 1.5 3])
        box off
        ylabel('TF');
        ylim([-1 9.5]); xlim([-0.5 3.5]);
        % text(0.25,1,'example trial with starting TF=6Hz & dTF=-3:1:3')
        h = text(-0.1,2, 'trial start'); set(h,'Rotation',90);
        h = text(3.1,2, 'trial end'); set(h,'Rotation',90);
        h = text(0.45,2, 'dTF start'); set(h,'Rotation',90);
        h = text(1.05,2, 'dTF stop'); set(h,'Rotation',90);
        subplot(10,3,[sf sf+3]);
        title(['trials with start_TF=' num2str(start_TF_values(sf)) 'Hz & dTF=-3:1:3; dTF dur. ' num2str(dTF_dur(fig)) 's']);

        %% plot mean population responses for negative, zero and positive dTF's
        ValidUnits = ~isnan(squeeze(SelectedResponses{sf,fig}(:,1,1)));
        SelectedUnits = SelectedResponses{sf,fig}(ValidUnits & AM_UOI,:,:);
        Responses_dTF_neg(:,1) = mean(nansum( (SelectedUnits(:,:,[1 2 3]))./nanmax(squeeze(nanmax(SelectedUnits(:,:,:),[],2)),[],2) ,3)/3);
        Responses_dTF_neg(:,2) = std(nansum( (SelectedUnits(:,:,[1 2 3]))./nanmax(squeeze(nanmax(SelectedUnits(:,:,:),[],2)),[],2) ,3)/3)/sqrt(length(ValidUnits));
        Responses_dTF_zer(:,1) = mean(nansum( (SelectedUnits(:,:,[4]    ))./nanmax(squeeze(nanmax(SelectedUnits(:,:,:),[],2)),[],2) ,3)/1);
        Responses_dTF_zer(:,2) = std(nansum( (SelectedUnits(:,:,[4]    ))./nanmax(squeeze(nanmax(SelectedUnits(:,:,:),[],2)),[],2) ,3)/1)/sqrt(length(ValidUnits));
        Responses_dTF_pos(:,1) = mean(nansum( (SelectedUnits(:,:,[5 6 7]))./nanmax(squeeze(nanmax(SelectedUnits(:,:,:),[],2)),[],2) ,3)/3);
        Responses_dTF_pos(:,2) = std(nansum( (SelectedUnits(:,:,[5 6 7]))./nanmax(squeeze(nanmax(SelectedUnits(:,:,:),[],2)),[],2) ,3)/3)/sqrt(length(ValidUnits));
        subplot(10,3,[sf+6 sf+6+6]); % show trials with dTF duration of 0.5s
        h1 = shadedErrorBar( linspace(-1,4,299),...
            Responses_dTF_neg(:,1),...
            Responses_dTF_neg(:,2),'lineprops',{'b-','markerfacecolor','b'});
        hold on
        h2 = shadedErrorBar( linspace(-1,4,299),...
            Responses_dTF_pos(:,1),...
            Responses_dTF_pos(:,2),'lineprops',{'r-','markerfacecolor','r'});
        hold on
        h3 = shadedErrorBar( linspace(-1,4,299),...
            Responses_dTF_zer(:,1),...
            Responses_dTF_zer(:,2),'lineprops',{'k-','markerfacecolor','k'});
        plot([0.5 0.5],[0 1],'g--','LineWidth',2); plot([0.5+dTF_dur(fig) 0.5+dTF_dur(fig)],[0 1],'g--','LineWidth',2) % dTF period
        plot([0 0],[0 1],'k:','LineWidth',4); plot([3 3],[0 1],'k:','LineWidth',4); % trial period
        set(gca,'TickDir','out','XTickLabel',[],'XTick',[0 0.5 1 3])
        box off
        ylabel('mean norm response');
        xlim([-0.5 3.5]); ylim([0 0.4])
        legend({'mean pop. response dTF<0','mean pop. response dTF>0','mean pop. response dTF=0'})
        subplot(10,3,[sf+6 sf+6+6]);
        title(['norm mean pop. responses for start_TF=' num2str(start_TF_values(sf)) 'Hz; dTF dur. ' num2str(dTF_dur(fig)) 's']);


        %% show responses of all units for the three cases together, possibly the difference
        ValidUnits = ~isnan(squeeze(SelectedResponses{sf,fig}(:,1,1)));
        SelectedUnits = SelectedResponses{sf,fig}(ValidUnits & AM_UOI,:,:);
        Responses_dTF = nansum( (SelectedUnits(:,:,[1 2 3 5 6 7]))./nanmax(squeeze(nanmax(SelectedUnits(:,:,:),[],2)),[],2) ,3)/6;
        MeanResponse = mean(Responses_dTF(:,94:123),2);
        [B, I] = sort(MeanResponse,'descend');
        subplot(10,3,sf+15:3:27+sf);
        if Diff == 0
            imagesc(linspace(-1,4,299),1:size(Responses_dTF,1),Responses_dTF(I,:))
        else
            dTF_zero = SelectedUnits(:,:,4)./nanmax(squeeze(nanmax(SelectedUnits(:,:,:),[],2)),[],2);
            imagesc(linspace(-1,4,299),1:size(Responses_dTF,1),Responses_dTF(I,:)-dTF_zero(I,:))
        end
        set(gca,'TickDir','out','XTick',[0 0.5 1 3])
        box off
        xlabel('seconds'); ylabel('units');
        xlim([-0.5 3.5])
        subplot(10,3,sf+15:3:27+sf);
        title(['norm responses, start_TF=' num2str(start_TF_values(sf)) 'Hz; all dTFs; dTF dur. ' num2str(dTF_dur(fig)) 's']);
         
    end

end



end