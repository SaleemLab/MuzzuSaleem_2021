function VerifySignals(ACInfo,ES,TrialExs_i,PertExs_i)

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
plot(ACInfo.trialStartsEnds(:,1)-ACInfo.trialStartsEnds(1,1),ones(length(ACInfo.trialStartsEnds(:,1)),1)*0.3,'sr')
hold on
plot(ACInfo.trialStartsEnds(:,2)-ACInfo.trialStartsEnds(1,1),ones(length(ACInfo.trialStartsEnds(:,2)),1)*0.3,'sk')
hold on
plot(ACInfo.PertStartsEnds(:,1)-ACInfo.trialStartsEnds(1,1),ones(length(ACInfo.PertStartsEnds(:,1)),1)*0.4,'sg')
hold on
plot(ACInfo.PertStartsEnds(:,2)-ACInfo.trialStartsEnds(1,1),ones(length(ACInfo.PertStartsEnds(:,2)),1)*0.4,'sy')
hold on
plot(ES.Time-ES.Time(TrialExs_i(1,1)),ES.Contrast)
hold on
plot(ES.Time-ES.Time(TrialExs_i(1,1)),ES.TF/max(ES.TF))
hold on
plot(ES.Time(TrialExs_i(:,1))-ES.Time(TrialExs_i(1,1)),ones(length(TrialExs_i(:,1)),1)*0.25,'or')
hold on
plot(ES.Time(TrialExs_i(:,2))-ES.Time(TrialExs_i(1,1)),ones(length(TrialExs_i(:,1)),1)*0.25,'ok')
hold on
plot(ES.Time(PertExs_i(:,1))-ES.Time(TrialExs_i(1,1)),ones(length(PertExs_i(:,1)),1)*0.25,'og')
hold on
plot(ES.Time(PertExs_i(:,2))-ES.Time(TrialExs_i(1,1)),ones(length(PertExs_i(:,1)),1)*0.25,'oy')

legend({'OE - photodiode', ...
        'OE - trial start', ...
        'OE - trial end', ...
        'OE - pert start', ...
        'OE - pert end', ...
        'BV - constrast', ...
        'BV - TF', ...
        'BV - trial start', ...
        'BV - trial end', ...
        'BV - pert start', ...
        'BV - pert end', ...
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

