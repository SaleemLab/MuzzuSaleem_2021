%% example of formatting

% 3D matrix
AM_UnitResponses_smooth
% AIM: compute responses for spcific angles during PertON and OFF

PertON_Angle = AM_Param(:,:,2)==0 & AM_Param(:,:,3)==1;
PertOFF_Angle= AM_Param(:,:,2)==0 & AM_Param(:,:,3)==0;

t1 = repmat(PertON_Angle, 1, 1, size(AM_UnitResponses_smooth,3)); % repeating the selection across time
t2 = repmat(PertOFF_Angle, 1, 1, size(AM_UnitResponses_smooth,3)); % repeating the selection across time

TrialsResp_OI = t1.*AM_UnitResponses_smooth;
TrialsRespControl_OI = t2.*AM_UnitResponses_smooth;

TrialsResp_OI_2D = squeeze(nansum(TrialsResp_OI,1))./nansum(PertON_Angle,1);
TrialsRespControl_OI_2D = squeeze(nanmean(TrialsRespControl_OI,1))./nansum(PertOFF_Angle,1);


figure
plot(nansum(PertON_Angle,1),'.')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3D matrix
AM_Speed
% AIM: compute mean speed at different time points during trials 
% find only perturbation 
Trials_PertON = AM_Param(:,:,3)==1; % find indexes where perturbation is on
Trials_PertOFF= AM_Param(:,:,3)==0; % find indexes where perturbation is off

Trials_PertOFF_lim = logical(zeros(size(Trials_PertOFF,1),size(Trials_PertOFF,2)));

% generate 3D matrix for including time responses
t1 = double(repmat(Trials_PertON, 1, 1, size(AM_Speed,3))); % repeating the selection across time
t2 =  double(repmat(Trials_PertOFF, 1, 1, size(AM_Speed,3))); % repeating the selection across time

% nan all the zero values
t1(t1(:)==0) = nan; t2(t2(:)==0) = nan;

% select only the responses during selected trials (pert ON and OFF separately)
Trials_Speed_OI = t1.*AM_Speed;
TrialsControl_Speed_OI = t2.*AM_Speed;

% These are indexes the pert onset and offest times 
PertON_Timings = AM_Param(:,:,4);
PertOFF_Timings = AM_Param(:,:,5)

% select only the chunks from the triaals I need
test = AM_Speed(:,:,ind(1):ind(2));

% or I can: flatten the speed traces onto a 2D matrix while keeping info of the
% recording number and of the unit (not necessary though)
Trials_Speed_OI_2D = reshape(Trials_Speed_OI, size(Trials_Speed_OI,1)*size(Trials_Speed_OI,2),size(Trials_Speed_OI,3));
TrialsControl_Speed_OI_2D = reshape(TrialsControl_Speed_OI, size(TrialsControl_Speed_OI,1)*size(TrialsControl_Speed_OI,2),size(TrialsControl_Speed_OI,3));

% compute the mean speed at different point in times for these trials in
% perturbation trials
AM_SpeedResponseControl_t(:,1) = mean(TrialsControl_Speed_OI_2D(:,trialSide_samples:trialSide_samples+2*BonVision_SR),2);
AM_SpeedResponseControl_t(:,2) = mean(TrialsControl_Speed_OI_2D(:,ind(1)-trialSide_samples/2:ind2(2)),2);
AM_SpeedResponseControl_t(:,3) = mean(TrialsControl_Speed_OI_2D(:,ind2(1)-trialSide_samples/2:ind2(2)+trialSide_samples),2);
AM_SpeedResponseControl_t(:,4) = mean(TrialsControl_Speed_OI_2D(:,trialSide_samples:end-trialSide_samples),2);

