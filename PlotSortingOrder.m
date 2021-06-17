%% Tomaso Muzzu - UCL - 13/09/2019

% order responses to plot

function [B,I, PertTrial_Ex] = PlotSortingOrder(ProjectData,SelectedResponses,AM_Param,AM_Speed,Units_Sel,Diff)

trialSide_samples = 60;
trialSide_seconds = 1;

Param_ = AM_Param(:,Units_Sel,:);
PertOnsets = Param_(:,:,4); % 2D matrix
PertOffsets = Param_(:,:,5); % 2D matrix
[v min_i(1)] = max(PertOnsets(:));
[v min_i(2)] = max(PertOffsets(:));

ReferSize = [size(Param_,1),size(Param_,2)];
[trial_el(1), trial_el(2)] = ind2sub(ReferSize,min_i(1)); % 2D index of minimum pert onset time

Rec = Param_(trial_el(1), trial_el(2),1);
% define timeline of example recording
TrialStart = ProjectData.Session_data{Rec,1}.ACInfo{1,1}.trialStartsEnds(trial_el(1),1);
TrialEnd = ProjectData.Session_data{Rec,1}.ACInfo{1,1}.trialStartsEnds(trial_el(1),2)-TrialStart;
TrialStart = 0;
TimeLine = linspace(TrialStart, ...
                    TrialEnd,...
                    size(SelectedResponses{1},2)); % -1 seconds
[v po_i(1)] = min(abs(TimeLine-min(PertOnsets(:))));
[v po_i(2)] = min(abs(TimeLine-min(PertOffsets(:))));

PertResponse = SelectedResponses{1};
NoPertResponse = SelectedResponses{2};
if Diff ==0
    SortingColumns = [mean(PertResponse(:,po_i(1):po_i(2)-30),2), ... % perturbation period
                      mean(PertResponse(:,trialSide_samples:trialSide_samples+trialSide_samples*2),2),... % first 2 seconds of trial
                      mean(PertResponse(:,po_i(2):po_i(2)+trialSide_samples),2)];   % second after perturbation
else
    SortingColumns = [abs(mean(PertResponse(:,po_i(1):po_i(2)-30),2)-mean(NoPertResponse(:,po_i(1):po_i(2)-30),2)), ... % perturbation period
                      abs(mean(PertResponse(:,trialSide_samples:trialSide_samples+trialSide_samples*2),2)-mean(NoPertResponse(:,trialSide_samples:trialSide_samples+trialSide_samples*2),2)) ,... % first 2 seconds of trial
                      abs(mean(PertResponse(:,po_i(2):po_i(2)+trialSide_samples),2)-mean(NoPertResponse(:,po_i(2):po_i(2)+trialSide_samples),2))];   % second after perturbation
end
[B,I] = sortrows(SortingColumns,[-1 -2 -3]);

% save coordinates of the example trial with perturbation
PertTrial_Ex(1) = Rec; % recording nr
PertTrial_Ex(2) = trial_el(1); % trial nr
PertTrial_Ex(3) = po_i(1); % trial nr
PertTrial_Ex(4) = po_i(2); % trial nr



end