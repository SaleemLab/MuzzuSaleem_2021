%% Tomaso Muzzu - UCL - 16/01/2019

% function to plot general info about sessions of single mice

function Stats = GeneralSessionStats_plot(Array,ES,StimulusParams,FileNames, RunningThreshold, MaxSpeed)

% Reminder of how Array is structured
% Array{i,1} = i; % order of trials of display
% Array{i,2} = ES.OrienTrial(i); % grating direction
% Array{i,3} = ind_1(i); % order of trials by ranking of choice
% Array{i,4} = SpeedTrial(i,:); % speed of every trial
% Array{i,5} = ES.PertTrial(i); % perturbation
% Array{i,6} = ind_2(cc);
% Array{i,7} = SpeedTrial_PertStartEnds(cc,:);

if ~exist('MaxSpeed','var')
    MaxSpeed = 50;
end


figure
set(gcf,'Renderer', 'painters', 'Position', [-1000 0 1050 500])
p1 = subplot(3,3,[1 2 3])
% mouse name
text(0,0.9,['Mouse: ' FileNames{1,1}(strfind(FileNames{1,1},'Cont_')-14:strfind(FileNames{1,1},'Cont_')-2)]);
% session nr and date
text(0,0.8,['Hab. sessions: ' num2str(size(ES,1)) ', from ' FileNames{1,1}(end-20:end-13) ' to ' FileNames{1,1}(end-20:end-13)]);
% number of trials
text(0,0.7,['Tot nr trials: ' num2str(size(Array,1)) ', trials/session: ' num2str(diff([find(Array.Trial_ID==1);length(Array.Trial_ID)]'))]);
% session duration
for i = 1:size(ES,1) 
    SessionDurations(1,i) = ES(i).Time(end);
end
text(0,0.60,['Session durations: ' num2str([floor(SessionDurations/60)+(floor(mod(SessionDurations,60))/100)]) ' min']);
% nr of trials for each direction
GratingDirections = unique(Array.Direction);
for i = 1:length(GratingDirections)
    GratingDirections(i,2) = sum(Array.Direction==GratingDirections(i)); 
end
text(0,0.5,['Grating directions: ' num2str(GratingDirections(:,1)')]);
text(0,0.40,['Trials per grating dir.: ' num2str(GratingDirections(:,2)')]);
% parameters of the grating: TF, SF, dimension and position of stimulus
text(0,0.3,['max contrast: ' num2str(max(ES(1).Contrast)) ', TF: ' num2str(max(ES(1).TF)) ', SF: ' num2str(max(ES(1).SF))]);
text(0,0.2,['Stimulus span(x,y): ' num2str([max(ES(1).X_Ext) max(ES(1).Y_Ext)])]); 
text(0,0.1,['Stimulus position(x,y): ' num2str([min(ES(1).X_Pos+1) min(ES(1).Y_Pos+1)])]);
box off
set(p1,'Visible','off')

% Gaussian filter
GaussFilter_W = 0.6; % seconds
BonvisionFR = (ES(1).Time(end)-ES(1).Time(1))/length(ES(1).Time); 
Width = round(GaussFilter_W/BonvisionFR);
Sigma = Width/3; % standard deviation in number of samples (converted from time in seconds)
x_g = linspace(-Width/2, Width/2, Width);
gaussFilter = exp(-x_g.^2/(2*Sigma^2));
gaussFilter_ = gaussFilter / sum (gaussFilter); % normalize
for i=1:size(ES,1)
    MouseSpeed = conv(ES(i).MouseSpeed, gaussFilter_, 'same');
    PerCentRun(i) = (sum((MouseSpeed(MouseSpeed>RunningThreshold & MouseSpeed<MaxSpeed)))/length(MouseSpeed));
    if PerCentRun(i) == 0
        MeanRunSpeed(i) = 0;
        SEMRunSpeed(i) = 0;
        MaxRunSpeed(i) = 0;
    else
        MeanRunSpeed(i) = mean(MouseSpeed(MouseSpeed>RunningThreshold & MouseSpeed<MaxSpeed));
        SEMRunSpeed(i) = std(MouseSpeed(MouseSpeed>RunningThreshold & MouseSpeed<MaxSpeed))/sqrt(length(MouseSpeed(MouseSpeed>RunningThreshold & MouseSpeed<MaxSpeed)));
        MaxRunSpeed(i) = max(MouseSpeed(MouseSpeed>RunningThreshold & MouseSpeed<MaxSpeed));
    end
        clear MouseSpeed
end
% time spent running
p2 = subplot(3,3,[4 7])
plot(1:size(ES,1), PerCentRun,'-o')
xlabel('session'); ylabel(['% time running > ' num2str(RunningThreshold) ' cm/s']);
box off; set(gca,'TickDir', 'out');

% other speed stats
p3 = subplot(3,3,[5 8])
hold off
h1 = errorbar(1:size(ES,1), MeanRunSpeed, SEMRunSpeed,'k')
hold on
h2 = plot(1:size(ES,1), MaxRunSpeed,'b-o')
xlabel('session'); ylabel(['% time running > ' num2str(RunningThreshold) ' cm/s']);
box off; set(gca,'TickDir', 'out');
legend({'Mean speed (SEM)', 'Max speed'})


p4 = subplot(3,3,[6 9])
t4 = title('other plot for other stats?')
set(p4,'Visible','off')
set(t4,'Visible','on')

a = axes;
MouseName_init = strfind(FileNames{1,1},'\M');
t1 = title(['Speed profiles during ' num2str(size(Array,1)) ' trials of ' num2str(length(FileNames)) ' sessions, '  ...
              FileNames{1,1}(MouseName_init+1:MouseName_init+1+6) 'TM' FileNames{1,1}(MouseName_init+1+11:MouseName_init+1+12)]);
a.Visible = 'off'; % set(a,'Visible','off');
t1.Visible = 'on'; % set(t1,'Visible','on');

Stats = table({1},{1},{1},...
        'VariableNames',{'MeanRunSpeed';'SEMRunSpeed';'TimeRunning'})
Stats{1,1} = {MeanRunSpeed};
Stats{1,2} = {SEMRunSpeed};
Stats{1,3} = {PerCentRun};

%EOF
end