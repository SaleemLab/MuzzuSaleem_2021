function boxplot_scatter_2(varargin)

% nargin(varargin)

FigPlot=figure;
ModIndex_Pos=[];Category=[];
for ii=1:size(varargin,2)
    var=varargin{1,ii};
    ModIndex_Pos=[ModIndex_Pos;var];
    Category=[Category;ii*ones(size(var))];

    scatter(ii*ones(size(var)).*(1+(rand(size(var))-0.5)/10),var,20,'k','filled')
    hold on
end
boxplot(ModIndex_Pos,Category,'Colors','k');

% ModIndex_Pos=[var1;var2];
% Category=[ones(size(var1));2*ones(size(var2))];
% ArenaTimeFig=figure;
% boxplot(ModIndex_Pos,Category);
% hold on
% scatter(ones(size(var1)).*(1+(rand(size(var1))-0.5)/10),var1,'r','filled')
% scatter(ones(size(var1)).*(1+(rand(size(var1))-0.5)/10),var1,'r')
% scatter(2*ones(size(var2)).*(1+(rand(size(var2))-0.5)/10),var2,'b')
end