%% Tomaso Muzzu - UCL - 12/12/2018

%% Main script to convert CSV data to MAT
% 1) Get data from csv files saved from Bonsai

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% 1) Get data from csv files saved from Bonsai
clear all
CurrentFolder = cd;
DataFolder = [CurrentFolder(1:end-(length(CurrentFolder)-3)) 'Archive - saleemlab' filesep 'Data' ];
FileName = uigetfile_n_dir(DataFolder,'Select CSV file');
FileName
for ff = 1:length(FileName)
    VisStimLog = csvread(FileName{1,ff},1,0);
    fid = fopen(FileName{1,ff});
    Columns = textscan(fid,'%s',1);
    fclose(fid);
    VisStimLog_Header = strsplit(Columns{1,1}{1,1},',');
    for i = 1:length(VisStimLog_Header)
        ES.(VisStimLog_Header{i}) = VisStimLog(:,i);
    end
    
    % 3) Find Trial starts and stops
    TwoDiffContrast = diff(diff(ES.Contrast));
    TrialExs = find((TwoDiffContrast/max(TwoDiffContrast))>0.1);
    % find the wrong trial onset/offsets
    j = 1; k = 1; clear TrialExs_i
    for i = 1:length(TrialExs)
        % trial onsets
        if TrialExs(i)+400<length(ES.Contrast)
            if sum(ES.Contrast(TrialExs(i)+1:TrialExs(i)+400)>0)>=399
                if i == 1
                    TrialExs_i(j,1) = TrialExs(i)+1;
                    j = j+1;
                else
                    if TrialExs(i)-TrialExs_i(j-1,1)>500
                        TrialExs_i(j,1) = TrialExs(i)+1;
                        j = j+1;
                    end
                end
            end
        end
        % trial offsets
        if i>=2
            if sum(ES.Contrast(TrialExs(i)-400:TrialExs(i)+1)>0)>=399
                if k==1
                    TrialExs_i(k,2) = TrialExs(i)+1;
                    k = k+1;
                else
                    if TrialExs(i)-TrialExs_i(k-1,2)>500
                        TrialExs_i(k,2) = TrialExs(i)+1;
                        k = k+1;
                    end
                end
            end
        end
    end
%     figure
%     plot(ES.Time,ES.Contrast/max(ES.Contrast));
%     hold on
%     plot(ES.Time,ES.TF/max(ES.TF))
%     hold on
%     plot(ES.Time(3:end),TwoDiffContrast/max(TwoDiffContrast))
%     hold on
%     plot(ES.Time(TrialExs(1:2:end)+1),ones(length(TrialExs(1:2:end)),1)*0.4,'or')
%     hold on
%     plot(ES.Time(TrialExs_i(:,1)),ones(length(TrialExs_i(:,1)),1)*0.3,'or') % trial onsets
%     hold on
%     plot(ES.Time(TrialExs_i(:,2)),ones(length(TrialExs_i(:,2)),1)*0.3,'ob') % trial offsets
    
%     DiffTF = diff(ES.TF);
%     PertExs_i = [find(DiffTF<-1), find(DiffTF>1)];
    
    % [MinTrialDur MinTrialDur_i] = min(TrialExs_i(:,2)-TrialExs_i(:,1));
    % if MinTrialDur<mean((TrialExs_i(:,2)-TrialExs_i(:,1)))*0.5
    %     fprintf('\nPlease check the trial on/offsets are correctly detected!\n');
    %     % overalapping of one start with an end of trial
    %     TrialExs = [TrialExs(1:MinTrialDur_i*2-1); TrialExs(MinTrialDur_i*2+1:end)];
    %     TrialExs_i = [TrialExs(1:2:end-1)+1, TrialExs(2:2:end)+1];
    % end
    %
    % [MinPertDur MinPertDur_i] = min(PertExs_i(:,2)-PertExs_i(:,1));
    % if MinPertDur<mean((PertExs_i(:,2)-PertExs_i(:,1)))*0.5
    %     fprintf('\nPlease check the perturbation on/offsets are correctly detected!\n');
    % end
    
    % plot data to verify that on/offsets of trials and perturbation are
    % correctly detected:
    % downsamle timestamps and photodiode signals
%     VerifySignals(ACInfo,ES,TrialExs_i,PertExs_i);
    
    ES.trialStartsEnds = TrialExs_i;
    
    if ff == 1
        % Save all these in a file (this will be added with the spikes later)
        SavingDir = uigetfile_n_dir(DataFolder,'Choose save directory'); 
        save([SavingDir{1,1} filesep 'Cont_TF_grating_' FileName{1,ff}(end-20:end-4) '.mat'], 'ES', '-v7.3');
    else
        save([SavingDir{1,1} filesep 'Cont_TF_grating_' FileName{1,ff}(end-20:end-4) '.mat'], 'ES', '-v7.3');
    end
    clear ES
end





