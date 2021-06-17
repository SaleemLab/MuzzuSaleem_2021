%% Tomaso Muzzu - UCL - 30/04/2020  - save RF properties of units

% list all recordings
%% load data
if ~exist('ProjectData','var')
    [ProjectData AM_UnitResponses AM_Param AM_Speed AM_UOI SelectedResponses AM_UnitResponses_smooth] = LoadDataALL;
end
if filesep == '\' 
    addpath(['X:' filesep 'CODE' filesep 'DEV' filesep 'general' filesep 'SparseNoise']);
    addpath(['X:' filesep 'CODE' filesep 'STABLE' filesep 'OpenEphys_analysis_tools']);
else
    addpath(['/mnt/pfs09' filesep 'CODE' filesep 'DEV' filesep 'general' filesep 'SparseNoise']);
    addpath(['/mnt/pfs09' filesep 'CODE' filesep 'STABLE' filesep 'OpenEphys_analysis_tools']);
end

%%
if isfile('X:\CODE\DEV\HFSP\RF_temp\RF_mapping.mat')
    load('X:\CODE\DEV\HFSP\RF_temp\RF_mapping.mat')
else   
    if size(ProjectData,2) < 6
        % find the recordings that have sparse noise
        for i = 1:length(ProjectData.Mouse_name)
            SubjectDir = ['X:' filesep '' filesep 'DATA' filesep 'SUBJECTS' filesep ProjectData.Mouse_name{i}];
            ProjectData.Recording_date{i}(1:10);
            % go through all recordings
            ProcessedFiles = dir([SubjectDir filesep 'Processed']);
            f_n = 1; f_n_i = 0;
            while f_n<=length(ProcessedFiles)
                if strfind(ProcessedFiles(f_n).name,'SN')
                    d_i = strfind(ProcessedFiles(f_n).name,['ks_'])+1;
                    if strfind(ProcessedFiles(f_n).name(d_i+2:d_i+11),ProjectData.Recording_date{i}(1:10))
                        clear ES
                        load([ProcessedFiles(f_n).folder filesep ProcessedFiles(f_n).name])
                        ProjectData.SparseNoise{i} = ES;
                    end
                end
                f_n = f_n + 1;
            end
        end       
        % MUA
        plotSingleUnits = 0;
        for i = 1:size(ProjectData,1)
            clear m
            if ~isempty(ProjectData.SparseNoise{i})
                m = ReverseCorrelationSTA(ProjectData.SparseNoise{i},plotSingleUnits);
                m = ReverseCorrelationSTC(ProjectData.SparseNoise{i},m,plotSingleUnits);
                m = ForwardCorrelationOFFMap(ProjectData.SparseNoise{i},m,plotSingleUnits);
                m = ForwardCorrelationONMap(ProjectData.SparseNoise{i},m,plotSingleUnits);
                %         m = MultiMeasureRF_MUA(ProjectData.SparseNoise{i},m);
                %         m = MultiMeasureRF_MUA_smooth(ProjectData.SparseNoise{i},m);
                close all
                ProjectData.SparseNoise{i}.m = m;
            end
        end
        for i = 1:size(ProjectData,1)
            if ~isempty(ProjectData.SparseNoise{i})
                ProjectData.SparseNoise{i}.ACInfo  = [];
                ProjectData.SparseNoise{i}.SpikeInfo  = [];
            end
        end
        save(['X:' filesep 'DATA' filesep 'PROJECTS' filesep 'VisPerturbation' filesep 'AllDataHere_SN.mat'],'ProjectData','-v7.3')
    end
end

%% evaluate signigicance / fit gaussian
if ~isfile('X:\CODE\DEV\HFSP\RF_temp\RF_mapping.mat')    
    clear STA_sh STC_sh OFFMap_sh ONMap_sh
    plotSingleUnits = 0; clear RF
    for i = 1:size(ProjectData,1)
        if ~isempty(ProjectData.SparseNoise{i})
            clear ES
            ES = ProjectData.SparseNoise{i};
            %tic
            for j = 1:1000
                clear m
                ES.SN_sequence = ES.SN_sequence(:,:,randperm(size(ES.SN_sequence,3)));
                m = ReverseCorrelationSTA(ES,plotSingleUnits);
                m = ReverseCorrelationSTC(ES,m,plotSingleUnits);
                m = ForwardCorrelationOFFMap(ES,m,plotSingleUnits);
                m = ForwardCorrelationONMap(ES,m,plotSingleUnits);
                close all
                STA_sh{i}(:,:,:,j) = m.STA;
                STC_sh{i}(:,:,:,j) = m.STC;
                OFFMap_sh{i}(:,:,:,j) = m.OFFMap;
                ONMap_sh{i}(:,:,:,j) = m.ONMap;
                j
            end
            %toc
        end
        i
    end
    save('RF_mapping.mat','STA_sh','STC_sh','OFFMap_sh','ONMap_sh','RFMaps','-v7.3')
else
    if ~exist('STA_sh','var')
        load('X:\CODE\DEV\HFSP\RF_temp\RF_mapping.mat')
    end
end

%%
STA = {}; STC = {}; OffMap = {}; ONMap = {}; SN_sequence = {}; stim_dims = {};
clear RFMaps 
RFMaps = table(STA,STC,OffMap,ONMap,SN_sequence,stim_dims);
k=1;
% save in a linear array all units
for i = 1:size(ProjectData,1)
    UnitsSession = size(ProjectData.Units_Info{i},1);
    if ~isempty(ProjectData.SparseNoise{i})
        RFMaps.SN_sequence{k} = ProjectData.SparseNoise{i}.SN_sequence;
        RFMaps.stim_dims{k} = ProjectData.SparseNoise{i}.stim_dims;
        for j = 1:UnitsSession
            RFMaps.STA{k} = ProjectData.SparseNoise{i}.m.STA(:,:,j);
            RFMaps.STC{k} = ProjectData.SparseNoise{i}.m.STC(:,:,j);
            RFMaps.OffMap{k} = ProjectData.SparseNoise{i}.m.OFFMap(:,:,j);
            RFMaps.ONMap{k} = ProjectData.SparseNoise{i}.m.ONMap(:,:,j);
            RFMaps.Mouse{k} = ProjectData.Mouse_name{i};
            RFMaps.Rec_id{k} = ProjectData.Recording_date{i};
            k = k + 1;
        end  
    else
        for j = 1:UnitsSession
            RFMaps.STA{k} = nan;
            RFMaps.STC{k} = nan;
            RFMaps.OffMap{k} = nan;
            RFMaps.ONMap{k} = nan;
            RFMaps.Mouse{k} = ProjectData.Mouse_name{i};
            RFMaps.Rec_id{k} = ProjectData.Recording_date{i};
            k = k + 1;
        end
    end
end
RFMaps = RFMaps(AM_UOI,:);

%% compute chi-squared test -- discard nan valued frames?
% save in a table the shuffled responses of all the unit
STA = {}; STC = {}; OffMap = {}; ONMap = {}; SN_sequence = {}; stim_dims = {};
clear RFMaps_sh 
RFMaps_sh = table(STA,STC,OffMap,ONMap);
k = 1;
for i = 1:size(ProjectData,1)
    UnitsSession = size(ProjectData.Units_Info{i},1);
    if ~isempty(ProjectData.SparseNoise{i})
        for j = 1:UnitsSession
            RFMaps_sh.STA{k} = squeeze(STA_sh{i}(:,:,j,:));
            RFMaps_sh.STC{k} = squeeze(STA_sh{i}(:,:,j,:));
            RFMaps_sh.OffMap{k} = squeeze(STA_sh{i}(:,:,j,:));
            RFMaps_sh.ONMap{k} = squeeze(STA_sh{i}(:,:,j,:));
            k = k + 1;
        end
    else
        for j = 1:UnitsSession
            RFMaps_sh.STA{k} = nan;
            RFMaps_sh.STC{k} = nan;
            RFMaps_sh.OffMap{k} = nan;
            RFMaps_sh.ONMap{k} = nan;
            k = k + 1;
        end
    end
end
RFMaps_sh = RFMaps_sh(AM_UOI,:);
clear Chi_squared Chi_squared_sh
for i = 1:size(RFMaps_sh,1)
    if ~isempty(RFMaps.SN_sequence{i})
        Stimulus = RFMaps.SN_sequence{i};
        Stimulus = Stimulus/max(Stimulus(:)); % 0 is black, 255 is white
        % # presentations of stimulus per location for STA and STC
        M_sta = sum(Stimulus == 0 | Stimulus == 1,3);
        % # presentations for ON stimulus
        M_on = sum(Stimulus == 1,3);
        % # presentations for OFF stimulus
        M_off = sum(Stimulus == 0,3);
    end
    % STA
    E_i = sum(RFMaps.STA{i}(:))/sum(M_sta(:));
    O_i = RFMaps.STA{i}./M_sta; % this is the receptive field. 2D matrix
    Chi_squared(i,1) = sum(sum(((E_i-O_i).^2)/E_i));
    % STC
    E_i = sum(RFMaps.STC{i}(:))/sum(M_sta(:));
    O_i = RFMaps.STC{i}./M_sta; % this is the receptive field. 2D matrix
    Chi_squared(i,2) = sum(sum(((E_i-O_i).^2)/E_i));
    % ON map
    E_i = sum(RFMaps.ONMap{i}(:))/sum(M_on(:));
    O_i = RFMaps.ONMap{i}./M_on; % this is the receptive field. 2D matrix
    Chi_squared(i,3) = sum(sum(((E_i-O_i).^2)/E_i));
    % OFF map
    E_i = sum(RFMaps.OffMap{i}(:))/sum(M_off(:));
    O_i = RFMaps.OffMap{i}./M_off; % this is the receptive field. 2D matrix
    Chi_squared(i,4) = sum(sum(((E_i-O_i).^2)/E_i));
    
    for j = 1:size(RFMaps_sh.STA{i},3)
        % STA
        E_i = sum(sum(RFMaps_sh.STA{i}(:,:,j)))/sum(M_sta(:));
        O_i = RFMaps_sh.STA{i}(:,:,j)./M_sta;
        Chi_squared_sh(i,1,j) = sum(sum(((E_i-O_i).^2)/E_i));
        % STC
        E_i = sum(sum(RFMaps_sh.STC{i}(:,:,j)))/sum(M_sta(:));
        O_i = RFMaps_sh.STC{i}(:,:,j)./M_sta;
        Chi_squared_sh(i,2,j) = sum(sum(((E_i-O_i).^2)/E_i));
        % ON map
        E_i = sum(sum(RFMaps_sh.ONMap{i}(:,:,j)))/sum(M_on(:));
        O_i = RFMaps_sh.ONMap{i}(:,:,j)./M_on;
        Chi_squared_sh(i,3,j) = sum(sum(((E_i-O_i).^2)/E_i));
        % OFF map
        E_i = sum(sum(RFMaps_sh.OffMap{i}(:,:,j)))/sum(M_off(:));
        O_i = RFMaps_sh.OffMap{i}(:,:,j)./M_off;
        Chi_squared_sh(i,4,j) = sum(sum(((E_i-O_i).^2)/E_i));
    end       
end    
clear Chi_squared_test
for i = 1:size(Chi_squared,1)
    Chi_squared_test(i,:,1) = Chi_squared(i,:) > prctile(squeeze(Chi_squared_sh(i,:,:))',99);
    Chi_squared_test(i,:,2) = Chi_squared(i,:) > prctile(squeeze(Chi_squared_sh(i,:,:))',95);
end

% save units with not nan values
ValidUnits = ~isnan(Chi_squared_test);
ValidUnits = ValidUnits(:,1);

%% first 7 animals
if  size(ProjectData,1)>10
    CTRL_exp = 0;
    Animal_1st_idx = [1 5 7 12 15 24 31];   
    if ~exist('PertResp_units','var')
        % select only perturbation responsive units
        load('AUC_shuffled.mat')
        Sh_responses = AUC_shuffled(:,2:end);
        p_pert_th = prctile(Sh_responses(:),99);
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
    % select only perturbation responsive units
    load('AUC_shuffled_CTRL_1.mat')
    Sh_responses = AUC_shuffled(:,2:end);
    p_pert_th = prctile(Sh_responses(:),99);
    PertResp_units = (AUC_shuffled(:,1)>p_pert_th);
    % select only pos. modulated perturbation responsive units
    load('DM_pert_shuffled_CTRL.mat')
    DM = DM_sh(:,1);
    DM_sign_i(:,1) = DM>0;
    DM_sign_i(:,2) = DM<=0;
    % select only pos. modulated perturbation responsive units
    PertRespUnits_pos = PertResp_units & DM_sign_i(:,1);
    PertRespUnits_neg = PertResp_units & DM_sign_i(:,2);
end

% fit 2d gaussian on sign. resp. units and define RF
addpath('X:\CODE\DEV\NOW_ON_GITHUB\GenericReceptiveFieldAnalysis')
clear fitparamsON optionsON fitparamsOFF optionsOFF
for i = 1:size(RFMaps,1) 
    if ~isnan(RFMaps.ONMap{i})
        forfit_f0 = RFMaps.ONMap{i};
        if ~isempty(RFMaps.stim_dims{i})
            forfit_xval = RFMaps.stim_dims{i}(1)/size(forfit_f0,1)/2 : RFMaps.stim_dims{i}(1)/size(forfit_f0,1) : RFMaps.stim_dims{i}(1)-RFMaps.stim_dims{i}(1)/size(forfit_f0,1)/2;
            forfit_xval = repmat(forfit_xval,size(forfit_f0,2),1);
            forfit_yval = (RFMaps.stim_dims{i}(2)/size(forfit_f0,2)/2 : RFMaps.stim_dims{i}(2)/size(forfit_f0,2) : RFMaps.stim_dims{i}(2)-RFMaps.stim_dims{i}(2)/size(forfit_f0,2)/2)-30;
            forfit_yval = repmat(sort(forfit_yval','descend'),1,size(forfit_f0,1));
        end
        % ON map
        [fitparamsON(i,:), modelout, optionsON(i)] = ibn_fit2DGaussLsq(forfit_f0,forfit_xval,forfit_yval);
        % OFF map
        forfit_f0 = RFMaps.OffMap{i};
        [fitparamsOFF(i,:), modelout, optionsOFF(i)] = ibn_fit2DGaussLsq(forfit_f0,forfit_xval,forfit_yval);
    end
end

% print the pert. responsive units

for i = 1:72 % sum(Chi_squared_test(:,3))
    if mod(i,72) == 1
        figure
    end
    if Chi_squared_test(i,3)
        if mod(i,72) ~= 0
            subplot(8,9,mod(i,72))
        else
            subplot(8,9,72)
        end
      imagesc(RFMaps.ONMap{i})
    end
end


for i = 1:72 % sum(Chi_squared_test(:,3))
    if mod(i,72) == 1
        figure
    end
    if Chi_squared_test(i,3)
        if mod(i,72) ~= 0
            subplot(8,9,mod(i,72))
        else
            subplot(8,9,72)
        end
      imagesc(RFMaps.OffMap{i})
    end
end


for i = 1:72 % sum(Chi_squared_test(:,3))
    if mod(i,72) == 1
        figure
    end
    if Chi_squared_test(i,3)
        if mod(i,72) ~= 0
            subplot(8,9,mod(i,72))
        else
            subplot(8,9,72)
        end
      forfit_f0 = RFMaps.ONMap{i};
      if ~isempty(RFMaps.stim_dims{i})
        forfit_xval = RFMaps.stim_dims{i}(1)/size(forfit_f0,1)/2 : RFMaps.stim_dims{i}(1)/size(forfit_f0,1) : RFMaps.stim_dims{i}(1)-RFMaps.stim_dims{i}(1)/size(forfit_f0,1)/2;
      end
      plot(forfit_xval,sum(RFMaps.ONMap{i})./sum(RFMaps.OffMap{i}))
      hold on
      plot(forfit_xval,zeros(length(forfit_xval),1),'k')
    end
end

% compute the ratio for all units. save it and plot it as a heatmap
clear SumsRatio
k = 1
for i = 1:sum(Chi_squared_test(:,3))
    if Chi_squared_test(i,3) & Chi_squared_test(i,4)
        TempSum = sum(RFMaps.ONMap{i})./sum(RFMaps.OffMap{i});
%         TempSum = TempSum/max(TempSum);
        SumsRatio(k,1:length(TempSum)) = TempSum;
        k = k + 1;
    end
end

figure
imagesc(log10(sortrows(SumsRatio,[-1 -2 -3])))
% zlim([0.5 2])
ylabel('Units')
xlabel('Centre out Azimuth')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear SumsRatio
k = 1
RF_x = fitparamsON(Chi_squared_test(:,3),1);
for i = 1:sum(Chi_squared_test(:,3))
    if Chi_squared_test(i,3) & Chi_squared_test(i,4) & RF_x(i)>=40
        TempSum = sum(RFMaps.ONMap{i})./sum(RFMaps.OffMap{i});
%         TempSum = TempSum/max(TempSum);
        SumsRatio(k,1:length(TempSum)) = TempSum;
        k = k + 1;
    end
end

figure
imagesc(log10(sortrows(SumsRatio,[-1 -2 -3])))
zlim([-1 1])
ylabel('Units')
xlabel('Centre out Azimuth')
RedWhiteBlue;

SumsRatio40 = SumsRatio;
SumsRatio30 = SumsRatio;
SumsRatio30(SumsRatio30==Inf)=nan;
SumsRatio30(SumsRatio30==-Inf)=nan;
SumsRatio40(SumsRatio40==Inf)=nan;
SumsRatio40(SumsRatio40==-Inf)=nan;

nanmean((SumsRatio30(:)))
nanmean(SumsRatio40(:))

%% focus on ON and OFF maps
figure
semilogy(fitparamsOFF(Chi_squared_test(:,4)==1,1), ...
    fitparamsOFF(Chi_squared_test(:,4)==1,3),'ok')
hold on
semilogy(fitparamsON(Chi_squared_test(:,3)==1,1), ...
    fitparamsON(Chi_squared_test(:,3)==1,3),'or')
legend({'Off responses', 'On responses'})
xlabel('Azimuth (deg)');
ylabel('Gain of fit');


edges = 0:10:120;
[N_off, edges] = histcounts(fitparamsOFF(Chi_squared_test(:,4)==1,1),edges)
[N_on, edges] = histcounts(fitparamsON(Chi_squared_test(:,3)==1,1),edges)
figure
bar(mean(diff(edges))/2:mean(diff(edges)):max(edges)-mean(diff(edges))/2,[N_off; N_on],'grouped')
legend({'Off responses', 'On responses'})
xlabel('Azimuth')
ylabel('Units')
set(gca,'box','off','FontSize',13,'TickDir','out')


clear SumsRatio
k = 1
for i = 1:sum(Chi_squared_test(:,3))
    if Chi_squared_test(i,3) & Chi_squared_test(i,4)
        TempSum = sum(optionsON(i).interpmodelout)./sum(optionsOFF(i).interpmodelout);
        TempSum = TempSum/max(TempSum);
        SumsRatio(k,1:length(TempSum)) = TempSum;
        k = k + 1;
    end
end

figure
imagesc(sortrows(SumsRatio,[-1 -2 -3]))
ylabel('Neurons')
xlabel('Centre-out Azimuth (a.u.)')




%% Saving new data for Brice
% 1) provide logical values for units with FR>0.1 Hz during sparse noise
% 2) save logical values for chi-squared test Chi_squared [units x (STA,STC,ON,OFF)]
% 3) convert to table
% 4) provide instructions of how you fitted the 2D gaussian
%%%%%%%%%%%%%%%%%%%%%%%
k=1;
for i = 1:size(ProjectData,1)
    UnitsSession = size(ProjectData.Units_Info{i},1);
    if ~isempty(ProjectData.SparseNoise{i})
        RFMaps.SN_onsets{k} = ProjectData.SparseNoise{i,1}.SN_onsets;
        k = k + UnitsSession;
    else
        k = k + UnitsSession;
    end
end
clear I_dx
for i = 1:size(RFMaps,1)
    I_dx(i) = sum(isnan(RFMaps.STA{i,1}(:)));
end
RFMaps_ = RFMaps{I_dx==0,:};
Chi_squared_test_ = Chi_squared_test(I_dx==0,:,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1)
clear SpikeTimes
for i = 1:size(ProjectData,1)
    if ~isempty(ProjectData.SparseNoise{i})
        if i == 1
            SpikeTimes = ProjectData.SparseNoise{i}.StimSpiketimes;
        else
            SpikeTimes_temp = ProjectData.SparseNoise{i}.StimSpiketimes;
            SpikeTimes = [SpikeTimes SpikeTimes_temp];
        end
    else
        SpikeTimes = [SpikeTimes cell(1,size(ProjectData.Units_Info{i},1))];
    end
end
SpikeTimes_ = SpikeTimes(I_dx==0)';
BonvisionFR = 60;
Time_edges = linspace(0,300,300*BonvisionFR);
GaussFilter_W = 0.3; % seconds
Width = round(GaussFilter_W*BonvisionFR);
Sigma = Width/3; % standard deviation in number of samples (converted from time in seconds)
x_g = linspace(-Width/2, Width/2, Width);
gaussFilter = exp(-x_g.^2/(2*Sigma^2));
gaussFilter_ = gaussFilter / sum (gaussFilter); % normalize
clear UnitFR
for i = 1:length(SpikeTimes_)
    [SpikeCount_temp, Edges] = histcounts(SpikeTimes_{i,1},Time_edges);       
    SpikeCount = SpikeCount_temp*BonvisionFR; % to get Hz
    UnitFR(i,:) = conv(SpikeCount, gaussFilter_, 'same');       
end
I_dx_FR = mean(UnitFR,2)>0.1;
RFMaps_(:,end+1) = SpikeTimes_;
RFMaps_(:,end+1) = num2cell(I_dx_FR);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2)
RFMaps_(:,end+1) = num2cell(Chi_squared_test_,[2 3]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3)
RFMaps__ = [RFMaps_(:,7) RFMaps_(:,8) RFMaps_(:,6) RFMaps_(:,5) RFMaps_(:,9) RFMaps_(:,10) RFMaps_(:,11) RFMaps_(:,1:4)];
SparseNoiseData1 = cell2table(RFMaps__, 'VariableNames', ...
    {'Mouse','Rec_date','SN_dims', 'SN_frames', 'SN_onsets', 'SpikeTimes', 'FR01','STA','STC','ON','OFF'});    


Table = SparseNoiseData(:,1:end-1);
for i = 1:size(Table,1)
    Table.Chi_2_01_05{i,1}= squeeze(SparseNoiseData.Chi_2_01_05(i,:,:))';
end
    
SparseNoiseData = SparseNoiseData1;    
save('SparseNoiseData_Lin_973.mat','SparseNoiseData');

%% verify changes

SparseNoiseData_new = SparseNoiseData;
plot_x = 10; plot_y = 10;
plot_total = plot_x*plot_y;
units = [1:349];
for i = units % sum(Chi_squared_test(:,3))
    if mod(i,plot_total) == 1
        figure
    end
        if mod(i,plot_total) ~= 0
            subplot(plot_x,plot_y,mod(i,plot_total))
        else
            subplot(plot_x,plot_y,plot_total)
        end
      imagesc(SparseNoiseData_new.ON{i})
    
end
suptitle('ON maps - NEW')
for i = units % sum(Chi_squared_test(:,3))
    if mod(i,plot_total) == 1
        figure
    end
        if mod(i,plot_total) ~= 0
            subplot(plot_x,plot_y,mod(i,plot_total))
        else
            subplot(plot_x,plot_y,plot_total)
        end
      imagesc(SparseNoiseData_new.OFF{i})
    
end
suptitle('OFF maps - NEW')


for i = 1:72 % sum(Chi_squared_test(:,3))
    if mod(i,72) == 1
        figure
    end
    
        if mod(i,72) ~= 0
            subplot(8,9,mod(i,72))
        else
            subplot(8,9,72)
        end
      imagesc(SparseNoiseData.ON{i})
    
end
suptitle('ON maps - OLD')
for i = 1:72 % sum(Chi_squared_test(:,3))
    if mod(i,72) == 1
        figure
    end
    
        if mod(i,72) ~= 0
            subplot(8,9,mod(i,72))
        else
            subplot(8,9,72)
        end
      imagesc(SparseNoiseData.OFF{i})
    
end
suptitle('OFF maps - OLD')


