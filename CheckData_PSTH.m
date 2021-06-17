% manual PSTH  to verify units responses


BonvisionFR = 60;


GaussFilter_W = 0.2

Width = round(GaussFilter_W*BonvisionFR);
trialSide_seconds = 1; % take 60 samples before and after trial
trialSide_samples = 60; % take 60 samples before and after trial
Sigma = Width/3; % standard deviation in number of samples (converted from time in seconds)
x_g = linspace(-Width/2, Width/2, Width);
gaussFilter = exp(-x_g.^2/(2*Sigma^2));
gaussFilter_ = gaussFilter / sum (gaussFilter); % normalize

 TimeStart = ProjectData.Session_data{r,1}.MetaData{1,1}.lims(4-1)/ProjectData.Session_data{r,1}.SpikeInfo{1,1}{end,2};
 TimeEnd = (ProjectData.Session_data{r,1}.MetaData{1,1}.lims(4-1)+ProjectData.Session_data{r,1}.MetaData{1,1}.lims(4))/ProjectData.Session_data{r,1}.SpikeInfo{1,1}{end,2};

lims = [11966464,63891456,9791488,60755968,11285504,9989120]./30000 ;
 
 
r = ceil(rand(1)*10);
StartTime = ProjectData.Session_data{r,1}.Time{1,1}(1);
EndTime = ProjectData.Session_data{r,1}.Time{1,1}(end);
for k = 1:size(ProjectData.Units_Info{r,1}.Spiketimes,1)
    edges = linspace(StartTime,EndTime,(EndTime-StartTime)*BonvisionFR);
    [FiringRate edges] = histcounts(ProjectData.Units_Info{r,1}.Spiketimes{k,1},edges);
    FiringRateSM(k,:) = conv(FiringRate, gaussFilter_, 'same');
end

% find responses during trial onsets
for i = 1 : size(ProjectData.Session_data{r,1}.trialStartsEnds{1,1},1)
   
    ProjectData.Session_data{r,1}.trialStartsEnds{1,1}(i,1)-trialSide_seconds
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BonvisionFR = 60;
GaussFilter_W = 0.2;

Width = round(GaussFilter_W*BonvisionFR);
trialSide_seconds = 1; % take 60 samples before and after trial
trialSide_samples = 60; % take 60 samples before and after trial
Sigma = Width/3; % standard deviation in number of samples (converted from time in seconds)
x_g = linspace(-Width/2, Width/2, Width);
gaussFilter = exp(-x_g.^2/(2*Sigma^2));
gaussFilter_ = gaussFilter / sum (gaussFilter); % normalize


for i = 1:size(SpikesTrial,1)
    clear FiringRate FiringRateSM
    for j = 1:size(SpikesTrial,2)
        StartTime = ES.ACInfo.trialStartsEnds(j,1)-trialSide_seconds;
        EndTime = ES.ACInfo.trialStartsEnds(j,2)+trialSide_seconds;
        edges = linspace(StartTime,EndTime,299);
        [FiringRate edges] = histcounts(SpikesTrial{i,j},edges);
        FiringRateSM(j,:) = conv(FiringRate, gaussFilter_, 'same');
    end
    UnitMeanResponse(i,:) = mean(FiringRateSM);
end


% plot average responses as PSTH with raster plot
figure
% subplot(4,1,1)
plot(mean(UnitMeanResponse))

















