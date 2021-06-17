%% Tomaso Muzzu - UCL - 12/12/2018

%% version for Bonvision stimulus

%% photodiode signal: save original signal downsampled at 1kHz, timestamps of stimulus onsets and offsets
function [ACInfo] = OE_Analyse_PD_Signal(ACInfo,directory)
    
%     if (nargout<3) & exist([DIRname filesep 'AC_Info.mat'],'file')
%         load([DIRname filesep 'AC_Info.mat'],'stimON','stimOFF');
%         % else load the analog channels now, find when the timestamps of the VS
%         % onsets and offsets
%     else
%         display('Analogue channels have not been converted yet or some error reading the saved file. Please select folder containing relevant OE recording for recovering ACs info.');
        
        % data(:,1) = sync pulse signal
        % data(:,2) = photodiode signal
        % data(:,3) = signal A from rotary encoder
        % data(:,4) = signal B from rotary encoder
        % data(:,5) = Bonvision state signal
       
        
        % resample and normalise signal
        if size(ACInfo.Data,2)>5
            normPD = ACInfo.Data(:,4)/max(ACInfo.Data(:,4));
        else
            normPD = ACInfo.Data(:,2)/max(ACInfo.Data(:,2));
        end
        [temp_pks, temp_locs] = findpeaks(normPD);
        pks = temp_pks(temp_pks>0.5);
        locs = temp_locs(temp_pks>0.5);
        clear temp_pks temp_locs
        % verify position of peaks
%         figure
%         plot(ACInfo.Data(:,2)/max(ACInfo.Data(:,2)));
%         hold on
%         plot(locs,pks,'ro');
      
        
        % find first peaks signalling the onset of the white square
        PksTimeDifference = diff(locs);
        % minimum distance between trial onsets is 7.3+1 seconds and max is
        % 7.3+2.5 seconds. At 30kHz this is 8.3*30000.
        % Perturbation lasts 1 or 1.17 seconds and starts at 4 + 0.5-0.67
        % seconds since the beginning. Worst case is: starts at 4.67 and
        % ends at 5.84 which is 1.4 seconds from the end of the trial.
        temp_VSOnsets = find(PksTimeDifference>0.1*ACInfo.SamplingRateOE); 
        VSOnsetsIndecesOFF = [locs(1); locs(temp_VSOnsets+1)]; % save indexes of stimulus onsets
        % find last peaks signalling the offset of the white square
        VSOnsetsIndecesON = [locs(temp_VSOnsets-1); locs(end)]; % save indexes of stimulus offsets
%         if ~isempty(find(diff(VSOnsetsIndecesOFF(2:end))>11*ACInfo.SamplingRateOE)) || ~isempty(find(diff(VSOnsetsIndecesON(2:end-1))>11*ACInfo.SamplingRateOE))
%             fprintf('Error during photodiode signal extraction. Please check it manually \n');
%             return
%         end
        % find perturbation onsets/offsets
%         PerOnsetsIndicesON = find(diff(VSOnsetsIndecesOFF)<4.1*ACInfo.SamplingRateOE==1);
%         PerOnsetsIndicesOFF = find(diff(VSOnsetsIndecesON)<7.3*ACInfo.SamplingRateOE==1);
        % find perturbation onsets/offsets based on the difference between
        % stim onset/offset timing              
        % check the inter-trial interval to find false events
        DiffITI = abs(VSOnsetsIndecesON(2:end-1)-VSOnsetsIndecesOFF(2:end-1));
        % find events when sync square changed colour by mistake
            %Trials_FalsePert = find(DiffITI<0.9*ACInfo.SamplingRateOE);
        % find trials when contrast stayed at zero
        Trials_ContrastZero = find(DiffITI>9*ACInfo.SamplingRateOE);
        trialON_zero = ACInfo.Timestamps(VSOnsetsIndecesON(Trials_ContrastZero+1));
        trialOFF_zero = ACInfo.Timestamps(VSOnsetsIndecesOFF(Trials_ContrastZero+1));
            %VSOnsetsIndecesON(Trials_FalsePert+1)=[];
            %VSOnsetsIndecesOFF(Trials_FalsePert+1)= [];
        stimON = ACInfo.Timestamps(VSOnsetsIndecesON(1:end-1));
        stimOFF = ACInfo.Timestamps(VSOnsetsIndecesOFF(2:end));
        % compute duration of stimuli
        DiffStim = abs(VSOnsetsIndecesON(1:end-1)-VSOnsetsIndecesOFF(2:end)); % abs(stimON-stimOFF); 
        Trials_idx = find(DiffStim>7*ACInfo.SamplingRateOE);
        Pert_idx = find(DiffStim<7*ACInfo.SamplingRateOE & DiffStim>4*ACInfo.SamplingRateOE);
        PertON = stimOFF(Pert_idx);
        PertOFF = stimON(Pert_idx+1);
        trialON = sort([stimON(Pert_idx); stimON(Trials_idx)]);
        trialOFF = sort([stimOFF(Pert_idx+1); stimOFF(Trials_idx)]);
        
%         PerOnsetsIndicesON = find(diff(VSOnsetsIndecesOFF)<4.1*ACInfo.SamplingRateOE==1);
%         PerOnsetsIndicesOFF = find(diff(VSOnsetsIndecesON)<7.3*ACInfo.SamplingRateOE==1);
%         % save start and end times of perturbations
%         PertON = ACInfo.Timestamps(VSOnsetsIndecesOFF(PerOnsetsIndicesON)+1); % time
%         PertOFF = ACInfo.Timestamps(VSOnsetsIndecesON(PerOnsetsIndicesOFF(1:2:end)+1)); % time
%         PertON = round(PertON,2); PertOFF = round(PertOFF,2);
%         % save start and end times of trials
%         stimON = round(stimON,2); stimOFF = round(stimOFF,2);
%         clear trialON; trialON = 0; k = 1;
%         for i = 1:length(stimON)
%             if sum(stimON(i)==PertOFF)==0
%                 trialON(k) = stimON(i);
%                 k = k+1;
%             end
%         end
%         clear trialOFF; trialOFF = 0; k = 1;
%         for i = 1:length(stimOFF)
%              if sum(stimOFF(i)==PertON)==0
%                 trialOFF(k) = stimOFF(i);
%                 k = k+1;
%              end
%         end
        % verify the events are detected correctly:
        fprintf('Please verify that all trial and perturbation onsets/offsets are detected correctly\n');
        plot_TS_temp = resample(ACInfo.Timestamps,1,30);
        plot_PD_temp = resample(normPD,1,30);
        figure
        plot(plot_TS_temp,plot_PD_temp)
        hold on
        plot(trialON,ones(1,length(trialON))*0.35,'r.','MarkerSize',18)
        hold on
        plot(trialOFF,ones(1,length(trialOFF))*0.4,'k.','MarkerSize',18)
        hold on
        plot(PertOFF,ones(1,length(PertOFF))*0.4,'rd')
        hold on
        plot(PertON,ones(1,length(PertON))*0.35,'kd')
        hold on
        plot(trialON_zero,ones(1,length(trialON_zero))*0.25,'*g')
        hold on
        plot(trialOFF_zero,ones(1,length(trialOFF_zero))*0.25,'*g')
%         figure
%         plot(trialOFF - trialON,'.k')           
        
        ACInfo.trialStartsEnds = [trialON(2:end), trialOFF(2:end)];
        ACInfo.PertStartsEnds = [PertON, PertOFF];
        ACInfo.trialZeroContrast = [trialON_zero trialOFF_zero];
        % save Onsets and Offsets into relevant folder
        %save([DIRname filesep 'AC_Info.mat'],'-struct', 'ACInfo');
%     end
end    




