%% Tomaso Muzzu - UCL - 05/02/2019

% function to plot general stats from 


function MultiMouseStats_Behaviour(StatsM,RunningThreshold)


figure
set(gcf,'Renderer', 'painters', 'Position', [-1000 0 1050 500])

for i = 1:size(StatsM,1)

    % time spent running
    p2 = subplot(1,2,1)
    plot(1:length(StatsM.TimeRunning{i}), StatsM.TimeRunning{i},'-o')
    hold on
    
    % other speed stats
    p3 = subplot(1,2,2)
    h1 = errorbar(1:length(StatsM.MeanRunSpeed{i}), StatsM.MeanRunSpeed{i}, StatsM.SEMRunSpeed{i}) ...
        %'Color', (i/(size(StatsM,1)+1))*([1 1 1]));
    hold on
    
    LegendNames{i} = ['Mouse ' num2str(i)];
end
subplot(1,2,1)
legend(LegendNames,'Location','northwest')
xlabel('session'); ylabel(['% time running > ' num2str(RunningThreshold) ' cm/s']);
box off; set(gca,'TickDir', 'out','XTick',[1:length(StatsM.MeanRunSpeed{1})]);

subplot(1,2,2)
box off; set(gca,'TickDir', 'out');
xlabel('session'); ylabel(['mean speed [cm/s]']);
set(gca,'TickDir', 'out','XTick',[1:length(StatsM.MeanRunSpeed{1})]);

end
