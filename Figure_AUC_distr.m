%% Tomaso Muzzu - UCL - 27/03/2020     

%% compute AUC of logistic regressor. 
%% New version which considers one direction at a time or directions are normalised as to compensate for direction tuning


%% Functions to plot Figure 2 - direction tuning
if ~exist('ProjectData','var')
    [ProjectData AM_UnitResponses AM_Param AM_Speed AM_UOI SelectedResponses AM_UnitResponses_smooth] = LoadDataALL;
end

% first 7 animals
if  size(ProjectData,1)>10
    CTRL_exp = 0;
    Animal_1st_idx = [1 5 7 12 15 24 31];
    
    if ~exist('PertResp_units','var')
        % select only perturbation responsive units
        load('AUC_shuffled.mat')
        Sh_responses = AUC_shuffled(:,2:end);
        %p_pert_th = prctile(Sh_responses(:),99);
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
    
else
    CTRL_exp = 1;
    % naive animals
    Animal_1st_idx = [1 4 7];
    
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

PertOnsets = AM_Param(:,:,4); % 2D matrix
PertOffsets = AM_Param(:,:,5); % 2D matrix
[v min_i(1)] = max(PertOnsets(:)); % find the latest moment when pert onset happens
[v min_i(2)] = max(PertOffsets(:)); % find the latest moment when pert offset happens

ReferSize = [size(AM_Param,1),size(AM_Param,2)];
[trial_el(1), trial_el(2)] = ind2sub(ReferSize,min_i(1)); % 2D index of minimum pert onset time

Rec = AM_Param(trial_el(1), trial_el(2),1);
% define timeline of example recording
TrialStart = ProjectData.Session_data{Rec,1}.ACInfo{1,1}.trialStartsEnds(trial_el(1),1);
TrialEnd = ProjectData.Session_data{Rec,1}.ACInfo{1,1}.trialStartsEnds(trial_el(1),2)-TrialStart;
TrialStart = 0;
TimeLine = linspace(TrialStart-trialSide_seconds, ...
                    TrialEnd+trialSide_seconds,...
                    size(AM_UnitResponses_smooth,3)); % -1 seconds
[v po_i(1)] = min(abs(TimeLine-min(PertOnsets(:))));
[v po_i(2)] = min(abs(TimeLine-min(PertOffsets(:))));
SelectedCells = AM_UOI;
param_sel = AM_Param(:,SelectedCells,:);
UnitResponses_smooth = AM_UnitResponses_smooth(:,SelectedCells,:);


%% two options: 1) direction by direction and 2) normalise dir. responses
% 1) One direction VS all (kind of)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2) normalise dir. responses.
Directions = unique(param_sel(:,1,2));
Directions = Directions(~isnan(Directions));
AngleNorm = 1;


%% Old system using 6 metrics relative to the pert. responses
clear DiffResponsePercent_Pert DiffResponsePercent_PostPert DiffResponsePercent_Vis AllRespPeriod_Pert
for i = 1:size(UnitResponses_smooth,2)
    clear UnitResponse_ UnitResponse
    UnitResponse_ = squeeze(UnitResponses_smooth(:,i,:))./max(squeeze(UnitResponses_smooth(:,i,:)),[],2);
    UnitResponse_ = UnitResponse_(1:end-1,:);
    % normalise responses for each direction
    if AngleNorm == 1
        Angle_Info = param_sel(1:end-1,i,2);
        % Normalise each angle response by its max
        clear MaxDirectResp UnitResponse
        for m = 1:length(Directions)
            MaxDirectResp(m) = max(mean(UnitResponse_(Angle_Info==Directions(m),:)));
            UnitResponse(Angle_Info==Directions(m),:) = UnitResponse_(Angle_Info==Directions(m),:)/MaxDirectResp(m);
        end
    else
        UnitResponse = UnitResponse_;
    end
    for k = 1:size(UnitResponse,1)
        AllRespPeriod_Pert(k,i,:) = interp1(1:length(UnitResponse(k,po_i(1)-trialSide_samples:po_i(1))), ...
                                               UnitResponse(k,po_i(1)-trialSide_samples:po_i(1)), ...
                                               1:100);
    end
    
%     DiffResponsePercent_Pert(:,i,1) = (nanmean(UnitResponse(:,po_i(1):po_i(2)),2)./nanmean(UnitResponse(:,po_i(1)-trialSide_samples:po_i(1)),2)); % mean(pert)/mean(prepert)
%     DiffResponsePercent_Pert(:,i,2) = (nanmean(UnitResponse(:,po_i(1):po_i(2)),2)-nanmean(UnitResponse(:,po_i(1)-trialSide_samples:po_i(1)),2)); % mean(pert)-mean(prepert)
%     DiffResponsePercent_Pert(:,i,3) = (nansum(UnitResponse(:,po_i(1):po_i(2)),2)); % integral of mean(pert)
%     DiffResponsePercent_Pert(:,i,4) = (nansum(UnitResponse(:,po_i(1)-trialSide_samples:po_i(2)),2)); % integral of mean(pert) and preceding 1s
%     DiffResponsePercent_Pert(:,i,5) = (nanmean(UnitResponse(:,po_i(1):po_i(2)-trialSide_samples/2),2)-nanmean(UnitResponse(:,po_i(1)-trialSide_samples:po_i(1)),2))./...
%                                       (nanmean(UnitResponse(:,po_i(1):po_i(2)-trialSide_samples/2),2)+nanmean(UnitResponse(:,po_i(1)-trialSide_samples:po_i(1)),2)); % mean(pert)-mean(prepert)/mean(pert+max(prepert) a.k.a depth of modulation
%     DiffResponsePercent_Pert(:,i,6) = nanmean(UnitResponse(:,po_i(1):po_i(2)-trialSide_samples/2),2)./...
%                                       (nanmean(UnitResponse(:,po_i(1):po_i(2)-trialSide_samples/2),2)+nanmean(UnitResponse(:,po_i(1)-trialSide_samples:po_i(1)),2)); % mean(pert)/mean(pert)+mean(prepert) a.k.a modulation index
    DiffResponsePercent_Pert(:,i,1) = (nanmedian(UnitResponse(:,po_i(1):po_i(2)),2)./nanmedian(UnitResponse(:,po_i(1)-trialSide_samples:po_i(1)),2)); % mean(pert)/mean(prepert)
    DiffResponsePercent_Pert(:,i,2) = (nanmedian(UnitResponse(:,po_i(1):po_i(2)),2)-nanmedian(UnitResponse(:,po_i(1)-trialSide_samples:po_i(1)),2)); % mean(pert)-mean(prepert)
    DiffResponsePercent_Pert(:,i,3) = (nansum(UnitResponse(:,po_i(1):po_i(2)),2)); % integral of mean(pert)
    DiffResponsePercent_Pert(:,i,4) = (nansum(UnitResponse(:,po_i(1)-trialSide_samples:po_i(2)),2)); % integral of mean(pert) and preceding 1s
    DiffResponsePercent_Pert(:,i,5) = (nanmedian(UnitResponse(:,po_i(1):po_i(2)-trialSide_samples/2),2)-nanmedian(UnitResponse(:,po_i(1)-trialSide_samples:po_i(1)),2))./...
                                      (nanmedian(UnitResponse(:,po_i(1):po_i(2)-trialSide_samples/2),2)+nanmedian(UnitResponse(:,po_i(1)-trialSide_samples:po_i(1)),2)); % mean(pert)-mean(prepert)/mean(pert+max(prepert) a.k.a depth of modulation
    DiffResponsePercent_Pert(:,i,6) = nanmedian(UnitResponse(:,po_i(1):po_i(2)-trialSide_samples/2),2)./...
                                      (nanmedian(UnitResponse(:,po_i(1):po_i(2)-trialSide_samples/2),2)+nanmedian(UnitResponse(:,po_i(1)-trialSide_samples:po_i(1)),2)); % mean(pert)/mean(pert)+mean(prepert) a.k.a modulation index
    
    
    DiffResponsePercent_PostPert(:,i,1) = (nanmean(UnitResponse(:,po_i(2):po_i(2)+trialSide_samples),2)./nanmean(UnitResponse(:,po_i(1):po_i(2)),2)); % mean(postpert)/mean(epert)
    DiffResponsePercent_PostPert(:,i,2) = (nanmean(UnitResponse(:,po_i(2):po_i(2)+trialSide_samples),2)-nanmean(UnitResponse(:,po_i(1):po_i(2)),2)); % mean(postpert)-mean(pert)
    DiffResponsePercent_PostPert(:,i,3) = (nansum(UnitResponse(:,po_i(2):po_i(2)+trialSide_samples),2)); % integral of mean(postpert)
    DiffResponsePercent_PostPert(:,i,4) = (nansum(UnitResponse(:,po_i(1):po_i(2)+trialSide_samples),2)); % integral of mean(postpert) and preceding 1s
    DiffResponsePercent_PostPert(:,i,5) = (nanmean(UnitResponse(:,po_i(2):po_i(2)+trialSide_samples/2),2)-nanmean(UnitResponse(:,po_i(2)-trialSide_samples/2:po_i(2)),2))./...
                                          (nanmean(UnitResponse(:,po_i(2):po_i(2)+trialSide_samples/2),2)+nanmean(UnitResponse(:,po_i(2)-trialSide_samples/2:po_i(2)),2)); % mean(pert)-mean(prepert)/mean(pert)+mean(prepert) a.k.a depth of modulation
    DiffResponsePercent_PostPert(:,i,6) =  nanmean(UnitResponse(:,po_i(2):po_i(2)+trialSide_samples/2),2)./...
                                          (nanmean(UnitResponse(:,po_i(2):po_i(2)+trialSide_samples/2),2)+nanmean(UnitResponse(:,po_i(2)-trialSide_samples/2:po_i(2)),2)); % mean(pert)/mean(pert)+mean(prepert) a.k.a modulation index
%     
%   DiffResponsePercent_Vis(:,i) = (nanmean(UnitResponse(:,trialSide_samples:trialSide_samples+2*BonVision_SR),2)./nanmean(UnitResponse(:,1:trialSide_samples),2)).^2;

end


%% step 2: use these values with the actual response to fit a logistic model
Param = squeeze(AM_Param(1:end-1,SelectedCells,2:3));
% one VS all classifier 

mode = 1;
clear SingleUnitResponse AUC_units
for i = 1:size(DiffResponsePercent_Pert,2) 
    % perturbation period
    SingleUnitResponse = squeeze(DiffResponsePercent_Pert(:,i,:));
    % select trial responses
    Responses = Param(:,i,2);
    % scrub responses
    Responses = logical(Responses(~isnan(Responses)));
    % scrub predictor
    SingleUnitResponse = SingleUnitResponse(1:length(Responses),:);
    SingleUnitResponse(SingleUnitResponse==Inf,:) = 0;
    SingleUnitResponse(isnan(SingleUnitResponse(:))) = 0;
    switch mode
        case 1 % pert ON or OFF
            % apply logistic model
            mdl = fitglm(SingleUnitResponse,Responses,'Distribution','binomial','Link','logit');
            scores = mdl.Fitted.Probability;
            % step 3: compute the ROC and evaluate the AUC for each unit
            [X,Y,T,AUC1] = perfcurve(Responses,scores,1);
            AUC_units(i,1) = AUC1;
        case 2 % pert ON and specify direction or OFF
            Angles = Param(:,i,1); clear Responses_a
            for d_i = 1:length(Directions)
                a = find(Param(:,i,1)==Directions(d_i) & Responses==1);
                b = find(Param(:,i,1)==Directions(d_i) & Responses==0);
                mdl = fitglm(SingleUnitResponse([a; b],:),Responses([a; b]),'Distribution','binomial','Link','logit');
                scores = mdl.Fitted.Probability;
                % step 3: compute the ROC and evaluate the AUC for each unit
                [X,Y,T,AUC1] = perfcurve(Responses([a; b]),scores,1);
                AUC_units(i,d_i) = AUC1;
            end            
    end
    
    
%     % apply logistic model
%     mdl = fitglm(SingleUnitResponse,Responses,'Distribution','binomial','Link','logit');
%     scores = mdl.Fitted.Probability;
%     % step 3: compute the ROC and evaluate the AUC for each unit
%     [X,Y,T,AUC1] = perfcurve(Responses,scores,1);
%     AUC_units(i,1) = AUC1;
%     if AUC1>0.8 && plotted==0
%         hold on
%         plot(X,Y,'Color',[0.5 0.5 0.5],'LineWidth',2)
%         xlabel('False positive rate')
%         ylabel('True positive rate')
%         title('ROC for Classification by Logistic Regression')
%         saved = [AUC1 i];
%         plotted = 1;
%     end
    
    
    %% post perturbation period 
%     SingleUnitResponse = squeeze(DiffResponsePercent_PostPert(:,i,:));
%     % scrub predictor
%     SingleUnitResponse = SingleUnitResponse(1:length(Responses),:);
%     SingleUnitResponse(SingleUnitResponse==Inf,:) = 0;
%     SingleUnitResponse(isnan(SingleUnitResponse(:))) = 0;
%     % apply logistic model
%     mdl = fitglm(SingleUnitResponse,Responses,'Distribution','binomial','Link','logit');
%     scores = mdl.Fitted.Probability;
%     % step 3: compute the ROC and evaluate the AUC for each unit
%     [X,Y,T,AUC2] = perfcurve(Responses,scores,1);
%     AUC_units(i,2) = AUC2;
    
end



%% step 4: find cutoff value through random shuffling of Pert trials id
tic
clear AUC_shuffled_pp AUC_shuffled AUC
for i = 1:size(DiffResponsePercent_Pert,2)
    clear SingleUnitResponse
    SingleUnitResponse = squeeze(DiffResponsePercent_Pert(:,i,:));
    % select trial responses
    Responses = Param(:,i,2);
    % scrub responses
    Responses = logical(Responses(~isnan(Responses)));
    % scrub predictor
    SingleUnitResponse = SingleUnitResponse(1:length(Responses),:);
    SingleUnitResponse(SingleUnitResponse==Inf,:) = 0;
    SingleUnitResponse(isnan(SingleUnitResponse(:))) = 0;
    
    % post perturbation period
    SingleUnitResponse_pp = squeeze(DiffResponsePercent_PostPert(:,i,:));
    % scrub predictor
    SingleUnitResponse_pp = SingleUnitResponse_pp(1:length(Responses),:);
    SingleUnitResponse_pp(SingleUnitResponse_pp==Inf,:) = 0;
    SingleUnitResponse_pp(isnan(SingleUnitResponse_pp(:))) = 0;
    parfor sh_i = 1:1000
               
        if mode == 1
            % shuffle perturbation trals
            Responses_sh = Responses(randperm(length(Responses)));
            % apply logistic model
            mdl = fitglm(SingleUnitResponse,Responses_sh,'Distribution','binomial','Link','logit');
            scores = mdl.Fitted.Probability;
            % step 3: compute the ROC and evaluate the AUC for each unit
            [X,Y,T,AUC1] = perfcurve(Responses,scores,1);
            AUC_shuffled(i,sh_i) = AUC1;
        elseif mode == 2
%             Angles = Param(:,i,1); clear Responses_a
%             for d_i = 1:length(Directions)
%                 a = find(Param(:,i,1)==Directions(d_i) & Responses==1);
%                 b = find(Param(:,i,1)==Directions(d_i) & Responses==0);
%                 Responses_a = [ones(length(a),1) ; zeros(length(b),1)];
%                 mdl = fitglm(SingleUnitResponse([a; b],:),Responses_a(randperm(length(Responses_a))),'Distribution','binomial','Link','logit');
%                 scores = mdl.Fitted.Probability;
%                 % step 3: compute the ROC and evaluate the AUC for each unit
%                 [X,Y,T,AUC1] = perfcurve(Responses_a(randperm(length(Responses_a))),scores,1);
%                 AUC_shuffled(i,d_i,sh_i) = AUC1;
%             end  
        end
        
%         % post perturbation period
%         mdl = fitglm(SingleUnitResponse_pp,Responses_sh,'Distribution','binomial','Link','logit');
%         scores = mdl.Fitted.Probability;
%         % step 3: compute the ROC and evaluate the AUC for each unit
%         [X,Y,T,AUC2] = perfcurve(Responses,scores,1);
%         AUC_shuffled_pp(i,sh_i) = AUC2;
    end
    i
end
toc

%%
AUC_shuffled = cat(2,AUC_units(:,1),AUC_shuffled);
AUC_shuffled_pp = cat(2,AUC_units(:,2),AUC_shuffled_pp);

save('AUC_shuffled_CTRL_1a.mat','AUC_shuffled')

figure
plot(AUC_shuffled(i,:),'.')

% prctile(AUC_shuffled(:),95) = 0.5922  on 18/09/19 p_sh(1) = 0.5922;
% prctile(AUC_shuffled(:),99) = 0.6754  on 18/09/19 p_sh(2) = 0.6754;
p_sh = prctile(reshape(AUC_shuffled(:,2:end),1,size(AUC_shuffled(:,2:end),1)*size(AUC_shuffled(:,2:end),2)),99);
% prctile(AUC_shuffled_pp(:),95) = 0.6057  on 18/09/19 pp_sh(1) = 0.6057
% prctile(AUC_shuffled_pp(:),99) = 0.7022  on 18/09/19 pp_sh(2) = 0.7022
pp_sh = [prctile(AUC_shuffled_pp(:),95) prctile(AUC_shuffled_pp(:),99)];

sum(AUC_units(:,1)>p_sh)

%% get 99% thresholds for all angles of all directions

for c_i = 1:size(AUC_shuffled,1)
    for d_i = 1:length(Directions)
        p_th_a(c_i,d_i) = prctile(squeeze(AUC_shuffled(c_i,d_i,:)),99);
    end
end

PertRespUnits = (AUC_units-p_th_a)>0;
figure
bar(sum(PertRespUnits,2))
xlabel('Units'); ylabel('# of directions pert. resp.')
set(gca,'box','off','TickDir','out')

UI_idx = 44;
figure
plot((Directions),max(GratingDirRespV(UI_idx,:)+GratingDirResp_BL(UI_idx,:),mean(GratingDirResp_BL(UI_idx,:))),'k','LineWidth',3)
hold on
GR2Plot = GratingDirResp(UI_idx,:);
GR2Plot(GR2Plot(:)<0) = 0;
plot(Directions,max(GratingDirResp(UI_idx,:)+GratingDirResp_BL(UI_idx,:),mean(GratingDirResp_BL(UI_idx,:))),'r','LineWidth',3)
%polarplot([deg2rad(Directions); 0],[GratingDirResp(UI_idx,:)+GratingDirResp_BL(UI_idx,:) GratingDirResp(UI_idx,1)+GratingDirResp_BL(UI_idx,1)],'r')
hold on
% plot baseline
plot([0 315], [mean(GratingDirResp_BL(UI_idx,:)) mean(GratingDirResp_BL(UI_idx,:))],'k:','LineWidth',3)
set(gca,'TickDir','out','box','off','XTick',[0:45:315],'FontSize',13)
xlabel('Direction'); ylabel('Spikes/s');

hold on
for i = 1:size(PertRespUnits(UI_idx,:),2)
    if PertRespUnits(UI_idx,i)==1
        plot(Directions(i), max(GratingDirResp(UI_idx,i)+GratingDirResp_BL(UI_idx,i), mean(GratingDirResp_BL(UI_idx,:))),'or','LineWidth',3)
    end
end

ll = legend({'Grating stim. response','Mean pert. response'})
ll.Color = 'none'; ll.EdgeColor = 'none';


figure
histogram(sum(PertRespUnits,2))
set(gca,'TickDir','out','box','off','FontSize',13)
xlabel('# of directions pert. resp'); ylabel('Units');


%% self reference for each unit
for i = 1:length(AUC_units)
    AUC_p(i) = sum(AUC_units(i)>AUC_shuffled(i,:));
end

% show AUC for perturbation period
figure
plot(AUC_units(:,1),'r.')
hold on
plot(find(AUC_p>=size(AUC_shuffled,2)*0.99),AUC_units(find(AUC_p>=size(AUC_shuffled,2)*0.99),1),'ro')
plot([0 sum(SelectedCells)], [p_sh(1) p_sh(1)], 'k')
xlabel('units'); ylabel('AUC')
ll = legend({'AUC perturbation period','sign. AUC value (s.r.)','95% shuffled trial ID', '99% shuffled trial ID'})
ll = legend({'AUC perturbation period','95% shuffled trial ID', '99% shuffled trial ID'})
ll.Color = 'none'; ll.EdgeColor = 'none';
xlim([0 sum(SelectedCells)+1]); 
% show AUC for post perturbation period
% figure
% plot(AUC_units(:,2),'b.')
% hold on
% plot([0 sum(SelectedCells)], [pp_sh(1) pp_sh(1)], 'k')
% plot([0 sum(SelectedCells)], [pp_sh(2) pp_sh(2)], 'k--')
% xlabel('units'); ylabel('AUC')
% legend({'AUC post-pert period','95% shuffle', '99% shuffle'})
% xlim([0 sum(SelectedCells)+1]); 
% % show AUC1 vs AUC2
% figure
% scatterDiagHist(AUC_units(:,1),AUC_units(:,2),[-1:0.02:1],'.k')
% xlabel('AUC pert period')
% ylabel('AUC post pert period')

%%
% show shuffling process
figure
h = histogram(AUC_shuffled(:,2:end),[0:0.02:1],'FaceColor',[0.5 0.5 0.5],'Normalization','probability','EdgeColor','none')
hold on
h1 = histogram(AUC_units(:,1),[0:0.02:1],'FaceColor',[1 0 0],'Normalization','probability','EdgeColor','none')
plot([p_sh(1) p_sh(1)], [0 max(h.BinCounts/length(AUC_shuffled(:)))*1.1],'k') 
%plot([p_sh(2) p_sh(2)], [0 max(h1.BinCounts/length(AUC_units(:,1)))*1.1],'k--','LineWidth',3)
xlabel('AUC of shuffled data'); ylabel('probability')
ll = legend({['AUC of shuffled data n=' num2str(sum(AM_UOI)) 'k'],['AUC of units n=' num2str(sum(AM_UOI))], '99% shuffled data'})
ll.Color = 'none'; ll.EdgeColor = 'none';
box off
set(gca,'FontSize',13,'TickDir','out');

%% plot modulation indexes
param_sel = AM_Param(:,SelectedCells,:);
Param = squeeze(param_sel(:,:,3));
if ~exist('AUC_units','var')
   AUC_units = AUC_shuffled(:,1); 
   p_sh = prctile(reshape(AUC_shuffled(:,2:end),1,size(AUC_shuffled(:,2:end),1)*size(AUC_shuffled(:,2:end),2)),99)
end
% Depth of modulation
r= 1; clear DM DepthMod
for i = 1:size(Param,2)
    % select trial responses
    Responses = Param(:,i);
    % scrub responses
    Responses = logical(Responses(~isnan(Responses)));
    UnitResponse = mean(squeeze(UnitResponses_smooth(Responses==1,i,:)))*60;

    Dir = sign(sum(UnitResponse(po_i(1):po_i(1)+trialSide_samples))-sum(UnitResponse(po_i(1)-trialSide_samples:po_i(1))));
    switch Dir
        case 1
            DM(i) = (max(UnitResponse(po_i(1):po_i(1)+trialSide_samples))-max(UnitResponse(po_i(1)-trialSide_samples:po_i(1))))/...
                (max(UnitResponse(po_i(1):po_i(1)+trialSide_samples))+max(UnitResponse(po_i(1)-trialSide_samples:po_i(1))));
        case -1
            DM(i) = (min(UnitResponse(po_i(1):po_i(1)+trialSide_samples))-min(UnitResponse(po_i(1)-trialSide_samples:po_i(1))))/...
                (min(UnitResponse(po_i(1):po_i(1)+trialSide_samples))+min(UnitResponse(po_i(1)-trialSide_samples:po_i(1))));
    end
    
    if AUC_units(i,1)>p_sh
        DepthMod(r,:) = [DM(i) AUC_units(i) i];
        r = r + 1;
    end
    
%     parfor sh=1:1000
%         % shuffle perturbation trals
%         Responses_sh = Responses(randperm(length(Responses)));
%         UnitResponse = mean(squeeze(UnitResponses_smooth(Responses_sh==1,i,:)))*60;
%         DM_sh(i,sh) = (sum(UnitResponse(po_i(1):po_i(1)+trialSide_samples))-sum(UnitResponse(po_i(1)-trialSide_samples:po_i(1))))/...
%             (sum(UnitResponse(po_i(1):po_i(1)+trialSide_samples))+sum(UnitResponse(po_i(1)-trialSide_samples:po_i(1))));
%     end

end

% verify disttribution of shuffled DM's
figure
histogram(DM_sh,[-1:0.05:1])
DM_shhh = cat(2,DM',DM_sh);
DM_sh = DM_shhh;
save('DM_pert_shuffled_CTRL.mat','DM_sh')

[unit_nr(:,1) unit_nr(:,2)] = sort(DM);

161 158 185

