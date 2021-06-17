function VerifySignals_Accel(ACInfo,ES,TrialExs_i)

% resample and normalise signal
if size(ACInfo.Data,2)>5
    normPD = ACInfo.Data(:,4)/max(ACInfo.Data(:,4));
else
    normPD = ACInfo.Data(:,2)/max(ACInfo.Data(:,2));
end
plot_TS_temp = resample(ACInfo.Timestamps,1,30);
plot_PD_temp = resample(normPD,1,30);

figure
plot(plot_TS_temp(plot_TS_temp>ACInfo.trialStartsEnds(1,1) & plot_TS_temp<ACInfo.trialStartsEnds(end,2))-ACInfo.trialStartsEnds(1,1),...
     plot_PD_temp(plot_TS_temp>ACInfo.trialStartsEnds(1,1) & plot_TS_temp<ACInfo.trialStartsEnds(end,2)))
hold on
plot(ACInfo.trialStartsEnds(:,1)-ACInfo.trialStartsEnds(1,1),ones(length(ACInfo.trialStartsEnds(:,1)),1)*0.5,'sr')
hold on
plot(ACInfo.trialStartsEnds(:,2)-ACInfo.trialStartsEnds(1,1),ones(length(ACInfo.trialStartsEnds(:,1)),1)*0.5,'sk')
hold on
plot(ES.Time,ES.TF/max(ES.TF))
hold on
plot(TrialExs_i,ones(length(TrialExs_i),1)*0.4,'or')

legend({'OE - photodiode', ...
        'OE - trial start', ...
        'OE - trial end', ...
        'BV - norm. TF', ...
        'BV - trial start', ...
        });

TrialsDiff(:,1) = ES.Time(TrialExs_i(:,2))-ES.Time(TrialExs_i(:,1));
TrialsDiff(:,2) = ACInfo.trialStartsEnds(:,2)-ACInfo.trialStartsEnds(:,1)

PertDiff(:,1) = ES.Time(PertExs_i(:,2))-ES.Time(PertExs_i(:,1));
PertDiff(:,2) = ACInfo.PertStartsEnds(:,2)-ACInfo.PertStartsEnds(:,1);

figure
plot(1:length(TrialsDiff(:,2)),TrialsDiff(:,2)-TrialsDiff(:,1),'k.')
xlabel('trial'); ylabel('trial duration diff (BV-OE) - [s]');
title(['avg trial duration diff = ' num2str(mean(TrialsDiff(:,2)-TrialsDiff(:,1)))])

figure
plot(1:length(PertDiff(:,2)),PertDiff(:,2)-PertDiff(:,1),'k.')
xlabel('trial'); ylabel('pertubation duration diff (BV-OE) - [s]');
title(['avg pert duration diff = ' num2str(mean(PertDiff(:,2)-PertDiff(:,1)))])


% time difference between the trial starts
TrialsDiff(:,3) = (ES.Time(TrialExs_i(:,1))-ES.Time(TrialExs_i(1,1)))-(ACInfo.trialStartsEnds(:,1)-ACInfo.trialStartsEnds(1,1));
% time difference between the trial ends
TrialsDiff(:,4) = (ES.Time(TrialExs_i(:,2))-ES.Time(TrialExs_i(1,1)))-(ACInfo.trialStartsEnds(:,2)-ACInfo.trialStartsEnds(1,1));
% time difference between the pertubation starts
PertDiff(:,3) = (ES.Time(PertExs_i(:,1))-ES.Time(TrialExs_i(1,1)))-(ACInfo.PertStartsEnds(:,1)-ACInfo.trialStartsEnds(1,1));
% time difference between the perturbation end
PertDiff(:,4) = (ES.Time(PertExs_i(:,2))-ES.Time(TrialExs_i(1,1)))-(ACInfo.PertStartsEnds(:,2)-ACInfo.trialStartsEnds(1,1));

figure
plot(1:length(TrialsDiff(:,2)),TrialsDiff(:,3),'ro')
hold on
plot(1:length(TrialsDiff(:,2)),TrialsDiff(:,4),'ko')
xlabel('trial'); ylabel('events diff occurrence (BV-OE) - [s]');
legend({'Trial starts diff','Trial ends diff'})
title(['avg trial diff starts= ' num2str(mean(TrialsDiff(:,3))) ', avg trial diff ends= ' num2str(mean(TrialsDiff(:,4)))])

figure
plot(1:length(PertDiff(:,2)),PertDiff(:,3),'go')
hold on
plot(1:length(PertDiff(:,2)),PertDiff(:,4),'yo')
xlabel('trial'); ylabel('events diff occurrence (BV-OE) - [s]');
legend({'Perturbation starts diff','Perturbation ends diff'})
title(['avg pert. diff starts= ' num2str(mean(PertDiff(:,3))) ', avg pert. diff ends= ' num2str(mean(PertDiff(:,4)))])

end

