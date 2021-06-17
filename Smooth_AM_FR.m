%% Tomaso Muzzu - UCL - 12/09/2019

% smooth spike count in trials

function AM_UnitResponses_smooth = Smooth_AM_FR(ProjectData,AM_UnitResponses,varargin);

BonvisionFR = 60;

if length(varargin)==1
    GaussFilter_W = varargin{1}; % seconds
    Width = round(GaussFilter_W*BonvisionFR);
    Sigma = Width/3; % standard deviation in number of samples (converted from time in seconds)
else 
    GaussFilter_W = varargin{1}; % seconds
    Width = round(GaussFilter_W*BonvisionFR);
    Sigma = varargin{2}; % standard deviation in number of samples (converted from time in seconds)
end

trialSide_seconds = 1; % take 60 samples before and after trial
trialSide_samples = 60; % take 60 samples before and after trial


Sigma = Width/3; % standard deviation in number of samples (converted from time in seconds)
x_g = linspace(-Width/2, Width/2, Width);
gaussFilter = exp(-x_g.^2/(2*Sigma^2));
gaussFilter_ = gaussFilter / sum (gaussFilter); % normalize

AM_UnitResponses_smooth = nan(size(AM_UnitResponses,1),size(AM_UnitResponses,2),size(AM_UnitResponses,3));
for i = 1:size(AM_UnitResponses,2)
    for j = 1:sum(~isnan(AM_UnitResponses(:,i,1)),1)
        AM_UnitResponses_smooth(j,i,:) = conv(squeeze(AM_UnitResponses(j,i,:)), gaussFilter_, 'same');
    end
end


end