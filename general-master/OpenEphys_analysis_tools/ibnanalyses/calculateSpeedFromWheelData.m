% Function to calculate the wheel speed from the optical sensor output
% History
%   30/6/2019 SGS Wrote it based on 'getSpeedSignal' - TM, 10/4/2018

function WdataSpeed = calculateSpeedFromWheelData(Wdata1,Wdata2, options);
if ~exist('options', 'var') || isempty(options) || ~isfield(options,'ResolEncoder')
    % For converting wheel into units into cm
    options.ResolEncoder = 1024; % strobes/revolution
    options.WheelRadius = 10; % cm
    options.sigma = 50; % standard deviation in number of samples, in this case is 50ms
    options.Width = options.sigma*3; % convert size from seconds to number of samples
end
if ~exist('Wdata2', 'var') || isempty(Wdata2)
    options.UseOneSignal = 1;
end
% Some precomputations for smoothing speed signal to computer speed
WheelCircum = options.WheelRadius*2*pi;
x_g = linspace(-options.Width/2, options.Width/2, options.Width);
gaussFilter = exp(-x_g.^2/(2*options.sigma^2));
gaussFilter_ = gaussFilter / sum (gaussFilter); % normalize

% First resample optical encoder data down to 1 kHz (ie. decimate by a
% factor of 30) - TM used resample but for some reason that is producing
% NaNs for SGS
Wdata1_1kHz = decimate(decimate(Wdata1,10),3);
% The data is recorded as an analogue signal, but is a digital signal, so
% make it back to all or nothing (0 or 1)
% There may be variable scale factors on the data so find the mean and then
% adjust so digital
mean_1 = mean(Wdata1_1kHz);
Wdata1_1kHz = 0.5+0.5*sign(Wdata1_1kHz-mean_1);
% Count numbers of rising edges
Wdata1_1kHz_diff = diff(Wdata1_1kHz);
% Find rising edges
plocs_1 = find(Wdata1_1kHz_diff > 0.5);
% Find falling edges
nlocs_1 = find(Wdata1_1kHz_diff < -0.5);
% Now we bring the data together to calculate speed
RFEvents = zeros(size(Wdata1_1kHz));
temp_Position = RFEvents;
RFEvents(plocs_1) = 1; % if A is rising set value to 1
RFEvents(nlocs_1) = -1; % if A is falling set value to -1


if ~options.UseOneSignal
    Wdata2_1kHz = decimate(decimate(Wdata2,10),3);
    mean_2 = mean(Wdata2_1kHz);
    Wdata2_1kHz = 0.5+0.5*sign(Wdata2_1kHz-mean_2);
    Wdata2_1kHz_diff = diff(Wdata2_1kHz);
    plocs_2 = find(Wdata2_1kHz_diff > 0.5);
    nlocs_2 = find(Wdata2_1kHz_diff < -0.5);
    RFEvents(plocs_2) = 2; % if B is rising set value to 2
    RFEvents(nlocs_2) = -2; % if B is falling set value to -2
end

error ('SGS: I dont understand the following')
% % Now 
%     for RF_i = 1:length(temp_RFEvents) 
%         switch temp_RFEvents(RF_i)
%             case 1
%                 if temp_ACInfo.EncoderSignal1kHz(RF_i,2)==0
%                     temp_Position(RF_i)=1;
%                 elseif temp_ACInfo.EncoderSignal1kHz(RF_i,2)==5
%                     temp_Position(RF_i)=-1;
%                 end
%             case -1
%                 if temp_ACInfo.EncoderSignal1kHz(RF_i,2)==5
%                     temp_Position(RF_i)=1;
%                 elseif temp_ACInfo.EncoderSignal1kHz(RF_i,2)==0
%                     temp_Position(RF_i)=-1;
%                 end
%             case 2
%                 if temp_ACInfo.EncoderSignal1kHz(RF_i,1)==5
%                     temp_Position(RF_i)=1;
%                 elseif temp_ACInfo.EncoderSignal1kHz(RF_i,1)==0
%                     temp_Position(RF_i)=-1;
%                 end
%             case -2
%                 if temp_ACInfo.EncoderSignal1kHz(RF_i,1)==0
%                     temp_Position(RF_i)=1;
%                 elseif temp_ACInfo.EncoderSignal1kHz(RF_i,1)==5
%                     temp_Position(RF_i)=-1;
%                 end
%             otherwise
%                 temp_Position(RF_i) = 0;
%         end
%     end    
    
%%%%%
% Convert rotations to speed
temp_PositionSum = -cumsum(temp_Position); % signal A and B are inverted here.
temp_PositionSum = temp_PositionSum*(WheelCircum/ResolEncoder); % [cm]
temp_ACInfo.PositionINcm = temp_PositionSum';

% compute speed
framerate = ((max(temp_ACInfo.timestampsDownsampled)-min(temp_ACInfo.timestampsDownsampled))/length(temp_ACInfo.timestampsDownsampled))^(-1);
temp_speed = -temp_Position.*(WheelCircum/ResolEncoder).*framerate;
% smooth speed signal to computer speed
temp_speedSmoothed = conv(temp_speed, gaussFilter_, 'same');
temp_ACInfo.SpeedSpeedSmoothed = [temp_speed; temp_speedSmoothed]';
temp_ACInfo.SamplingRateSignalsOut = framerate;