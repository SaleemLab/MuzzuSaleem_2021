%% Tomaso Muzzu - UCL - 12/12/2018

%% version for Bonvision stimulus


%% speed signal: save original signal downsampled at 1kHz, position in cm and speed in cm/s at 1kHz
function ACInfo = OE_Compute_Speed_Signal(ACInfo)
    
    % load the mat file containing the speed signal if it already exists
%     if (nargout<3) & exist([DIRname filesep 'AC_Info.mat'],'file')
%         load([DIRname filesep 'AC_Info.mat']); 
%         %%
%         %%
%         %% TO BE CONTINUED ........................
%         %%
%         %%
%         % else load the analog channels now, find when the timestamps of the VS
%         % onsets and offsets
%     else
    data = ACInfo.Data;    
    if size(ACInfo.Data,2)>5
        temp_EncoderSignal(:,1) = data(:,5);
        temp_EncoderSignal(:,2) = data(:,6);
    else
        temp_EncoderSignal(:,1) = data(:,3);
        temp_EncoderSignal(:,2) = data(:,4);
    end
    temp_timestampsDownsampled = linspace(min(ACInfo.Timestamps),max(ACInfo.Timestamps),length(temp_EncoderSignal)); 
    
    % digitalise/normalise the signal
    temp_EncoderSignal(temp_EncoderSignal(:,1)>2.5,1) = 5; % originally it is a digital signal so it is 0-5V
    temp_EncoderSignal(temp_EncoderSignal(:,1)<=2.5,1) = 0; 
    temp_EncoderSignal(temp_EncoderSignal(:,2)>2.5,2) = 5; 
    temp_EncoderSignal(temp_EncoderSignal(:,2)<=2.5,2) = 0; 
    % save it
    temp_ACInfo.EncoderSignal = temp_EncoderSignal;
    % count numbers of rising edges 
    temp_DiffEncoderSignal(:,1) = diff(temp_EncoderSignal(:,1));
    temp_DiffEncoderSignal(:,2) = diff(temp_EncoderSignal(:,2));
    % find rising edges
    [temp_ppks,temp_plocs] = findpeaks( temp_DiffEncoderSignal(:,1));
    ppks_A = temp_ppks(temp_ppks>2.5); % select real edges [Volts]
    plocs_A = temp_plocs(temp_ppks>2.5);
    [temp_ppks,temp_plocs] = findpeaks( temp_DiffEncoderSignal(:,2));
    ppks_B = temp_ppks(temp_ppks>2.5); % select real edges [Volts]
    plocs_B = temp_plocs(temp_ppks>2.5);
    clear temp_ppks temp_plocs
    % find falling edges
    [temp_npks,temp_nlocs] = findpeaks(- temp_DiffEncoderSignal(:,1));
    npks_A = temp_npks(temp_npks>2.5); % select real edges [Volts]
    nlocs_A = temp_nlocs(temp_npks>2.5);
    [temp_npks,temp_nlocs] = findpeaks(- temp_DiffEncoderSignal(:,2));
    npks_B = temp_npks(temp_npks>2.5); % select real edges [Volts]
    nlocs_B = temp_nlocs(temp_npks>2.5);
    clear temp_npks temp_nlocs
    % 
    [v,ind] = min([min(plocs_A) min(nlocs_A) min(plocs_B) min(nlocs_B)]);
    temp_RFEvents = zeros(length(temp_ACInfo.EncoderSignal),1);
    temp_RFEvents(plocs_A) = 1; % if A is rising set value to 1
    temp_RFEvents(nlocs_A) = -1; % if A is falling set value to -1
    temp_RFEvents(plocs_B) = 2; % if B is rising set value to 2
    temp_RFEvents(nlocs_B) = -2; % if B is falling set value to -2
    clear temp_PositionDelta
%     A_set = 0;
%     B_set = 0;
%     for RF_i = 1:length(temp_RFEvents) 
%         switch temp_RFEvents(RF_i)
%             case 1
%                 A_set=1;
%                 if B_set == 0
%                     temp_PositionDelta(RF_i)=1;
% %                 elseif temp_EncoderSignal(RF_i,2)>3
% %                     temp_PositionDelta(RF_i)=-1;
%                 end
%             case -1
%                 A_set=0;
% %                 if temp_EncoderSignal(RF_i,2)>3
% %                     temp_PositionDelta(RF_i)=1;
% %                 elseif temp_EncoderSignal(RF_i,2)<2
% %                     temp_PositionDelta(RF_i)=-1;
% %                 end
%             case 2
% %                 if temp_EncoderSignal(RF_i,1)>3
% %                     temp_PositionDelta(RF_i)=1;
%                 B_set=1;
%                 if A_set == 0
%                     temp_PositionDelta(RF_i)=-1;
%                 end
%             case -2
%                 B_set=0;
% %                 if temp_EncoderSignal(RF_i,1)<2
% %                     temp_PositionDelta(RF_i)=1;
% %                 elseif temp_EncoderSignal(RF_i,1)>3
% %                     temp_PositionDelta(RF_i)=-1;
% %                 end
%             otherwise
%                 temp_PositionDelta(RF_i) = 0;
%         end
%     end
    EncoderPosition = zeros(length(temp_ACInfo.EncoderSignal),1); A_set = 0; B_set = 0;
    for RF_i = 1:length(temp_RFEvents) 
        switch temp_RFEvents(RF_i)
            case 1
                A_set = 1;
                if B_set == 0
                    EncoderPosition(RF_i)=1;
                end
            case -1
                A_set = 0;
            case 2
                B_set = 1;
                if A_set == 0
                    EncoderPosition(RF_i)=-1;
                end
            case -2
                B_set = 0;
        end
    end         
    % verify speed signal
%         figure
%         axis([ACInfo.Timestamps(1) ACInfo.Timestamps(end) -10 10])
%         plot(ACInfo.Timestamps,temp_EncoderSignal(:,1))
%         hold on
%         plot(ACInfo.Timestamps,temp_EncoderSignal(:,1),'k')
%         hold on
%         plot(ACInfo.Timestamps(2:end),temp_DiffEncoderSignal(:,1),'r')
%         hold on
%         plot(ACInfo.Timestamps(plocs_B), ppks_B,'ro')
%         hold on
%         plot(ACInfo.Timestamps(nlocs_B), npks_B,'ko')
%         hold on
%         plot(ACInfo.Timestamps,temp_EncoderSignal(:,2),'c')
%         ylim([-6 6])
%         hold on
%         plot(ACInfo.Timestamps,temp_RFEvents,'g')
%         hold on
%         plot(ACInfo.Timestamps,temp_Position,'y')
    
    %temp_PositionSum = -cumsum(temp_PositionDelta); % signal A and B are inverted here.
    temp_PositionSum = -cumsum(EncoderPosition); % signal A and B are inverted here.
    % convert units into cm
    ResolutionEncoder = 1024; % ticks per revolution
    WheelRadius = 10; % cm
    framerate = ((max(ACInfo.Timestamps)-min(ACInfo.Timestamps))/length(ACInfo.Timestamps))^(-1);
    WheelCircum = WheelRadius*2*pi; % cm
    temp_PositionSum = temp_PositionSum*(WheelCircum/ResolutionEncoder); % [cm]
    temp_ACInfo.Position_cm = temp_PositionSum'; 
    % compute speed
    %temp_speed = -temp_PositionDelta*(WheelCircum/ResolutionEncoder);
    temp_speed = -EncoderPosition*(WheelCircum/ResolutionEncoder);
    % smooth speed signal to computer speed
    sigma = 50; % standard deviation in number of samples, in this case is 50ms
    Width = sigma*3; % convert size from seconds to number of samples
    x_g = linspace(-Width/2, Width/2, Width);
    gaussFilter = exp(-x_g.^2/(2*sigma^2));
    gaussFilter_ = gaussFilter / sum (gaussFilter); % normalize
    temp_speedSmoothed = conv(temp_speed, gaussFilter_, 'same');    
    temp_ACInfo.SpeedSpeedSmoothed = [temp_speed; temp_speedSmoothed]';
    temp_ACInfo.SamplingRateSignalsOut = framerate;
    % verify that position increases in general as the mice move forward mainly 
    
%     figure
%     plot(ACInfo.Timestamps-min(ACInfo.Timestamps),temp_PositionSum)
%         plot(ACInfo.Timestamps-min(ACInfo.Timestamps),-cumsum(EncoderPosition))
%     figure
%     plot(ACInfo.Timestamps-min(ACInfo.Timestamps),temp_speed)
%     hold on
%     plot(ACInfo.Timestamps-min(ACInfo.Timestamps),temp_speedSmoothed*framerate,'k')  
    
    ACInfo.RotEncSpeed = [temp_speed, temp_speedSmoothed];
    ACInfo.RotEncPos = temp_PositionSum;
    clear temp_ACInfo;
end