%% Tomaso Muzzu - UCL - 28/05/2020 - plot ori and dir angles 

%% subordinate to Figure_2_D

Cycles2plot = 2; 
OriMax = 180;
%% ORIENTATION
figure
BinningDir = 0:15:Cycles2plot*OriMax;
p_thres = 0.01; 
%set(gcf,'Position',[100 100 650 600])
subplot(4,4,[9 10 13 14])
title('Ori. tuned units (vis.)')
%plot(TuningPropsV.P_ori,TuningProps.P_ori,'.k')
%load('AUC_shuffled.mat')
%Sh_responses = AUC_shuffled(:,2:end);
p_pert_th = prctile(Sh_responses(:),99);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PertResp_units = (AUC_shuffled(:,1)>p_pert_th);
%PertResp_units = UnitsPertPos;
PertResp_units = PertRespUnits_pos;
PertResp_units_V = PertRespUnits_pos & TuningPropsV.p_HT2<p_thres;
%PertResp_units = PertResp_units_V;
%PertResp_units = 1:189; PertResp_units = logical(PertResp_units');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on
% plot(TuningProps.P_ori(PertResp_units),TuningPropsV.P_ori(PertResp_units),'or','MarkerSize',5)
plot(TuningProps.P_ori(PertResp_units_V),TuningPropsV.P_ori(PertResp_units_V),'or','MarkerSize',5,'markerfacecolor','r')
for Cy = 2:Cycles2plot
    % plot higher angles for vis. stim. responses
    plot(TuningProps.P_ori(PertResp_units_V),TuningPropsV.P_ori(PertResp_units_V)+OriMax*(Cy-1),'or','MarkerSize',5,'markerfacecolor','r')
    % plot higher angles for pert. responses
    plot(TuningProps.P_ori(PertResp_units_V)+OriMax*(Cy-1),TuningPropsV.P_ori(PertResp_units_V),'or','MarkerSize',5,'markerfacecolor','r')
    % bottom left corner
    plot(TuningProps.P_ori(PertResp_units_V)+OriMax*(Cy-1),TuningPropsV.P_ori(PertResp_units_V)+OriMax*(Cy-1),'or','MarkerSize',5,'markerfacecolor','r')
end
ylabel('Visual preferred orientation'); xlabel('Perturbation max resp. orientation'); 
set(gca,'YAxisLocation', 'left','XTick',[0:45:Cycles2plot*OriMax],'XTickLabel',[0:45:Cycles2plot*OriMax],'YTick',[0:45:Cycles2plot*OriMax],'YTickLabel',[0:45:Cycles2plot*OriMax],'TickDir','out'); 
box off;
xlim([0 Cycles2plot*OriMax])
ylim([0 Cycles2plot*OriMax])

subplot(4,4,[1 2 5 6])
DataCat = TuningProps.P_ori(PertResp_units_V);
for Cy = 2:Cycles2plot
    DataCat = [DataCat; TuningProps.P_ori(PertResp_units_V)+OriMax*(Cy-1)]
end
histogram(DataCat, BinningDir,'EdgeColor','k','FaceColor',[229 44 37]/255)
set(gca,'YAxisLocation','left','XTick',[0:45:Cycles2plot*OriMax],'XTickLabel',[],'TickDir','out'); box off;
ylabel('Units')
xlim([0 Cycles2plot*OriMax])

subplot(4,4,[11 12 15 16])
DataCat = TuningPropsV.P_ori(PertResp_units_V);
for Cy = 2:Cycles2plot
    DataCat = [DataCat; TuningPropsV.P_ori(PertResp_units_V)+OriMax*(Cy-1)]
end
histogram(DataCat, BinningDir,'EdgeColor','k','FaceColor',[229 44 37]/255)
set(gca,'XTick',[0:45:Cycles2plot*OriMax],'XTickLabel',[],'TickDir','out','XAxisLocation', 'bottom'); box off;
ylabel('Units');
set(gca, 'XDir','reverse')
xlim([0 Cycles2plot*OriMax])
view(90,90)
hold off

% chance level
Distr1 = TuningPropsV.P_ori(PertResp_units_V); % preferred direction for grating
Distr2 = TuningProps.P_ori(PertResp_units_V);  % preferred direction for perturbation
clear Distr_uni N_sh
for sh = 1:1000
    Distr1_sh = Distr1(randperm(length(Distr1)));
    Distr2_sh = Distr1(randperm(length(Distr2)));
    Distr_uni(sh,:) = rad2deg(angdiff(deg2rad(Distr1_sh*2),deg2rad(Distr2_sh*2)))/2;
    %Distr_uni_(sh,:) = histcounts(Distr1_sh-Distr2_sh,-90:18:90);
    %[N_sh(:,:,sh), Xedges, Yedges] = histcounts2(Distr1_sh,Distr2,0:10:180,0:10:180);
end
subplot(4,4,[3 4 7 8])
hold off
Ori_diff = rad2deg(angdiff(deg2rad(TuningPropsV.P_ori(PertResp_units_V)*2),deg2rad(TuningProps.P_ori(PertResp_units_V)*2)))/2;
histogram(Ori_diff,-90:18:90,'EdgeColor','r','FaceColor',[229 44 37]/255,'Normalization','probability')
hold on
histogram(Distr_uni,-90:18:90,'EdgeColor','k','FaceColor',[0.5 0.5 0.5],'Normalization','probability')
% pd = fitdist(TuningPropsV.P_ori(PertResp_units)-TuningProps.P_ori(PertResp_units),'Normal')
% hold on
% x_values = -180:1.8:180;
% y_values = (1/(pd.sigma*sqrt(2*pi))).*exp((-(x_values-pd.mean).^2)./(2*pd.sigma^2));
% plot(x_values,y_values,'LineWidth',2)
xlim([-90 90])
set(gca,'xTick',-90:45:90,'TickDir','out','View',[45,90])
box off;
xlabel('Orientation'); ylabel('Units')
set(gcf,'Position',[10 10 800 700])

% [p h] = signrank(Ori_diff)


% figure
% histogram(Distr_uni,-90:18:90,'EdgeColor','k','FaceColor',[0.5 0.5 0.5],'Normalization','probability')
% hold on
% histogram(Ori_diff,-90:18:90,'EdgeColor','k','FaceColor',[229 44 37]/255,'Normalization','probability')
% xlim([-180 180])
% set(gca,'xTick',-90:45:90,'TickDir','out','View',[45,90])
% box off;
% xlabel('Orientation'); ylabel('Units')
% set(gcf,'Position',[10 10 800 700])
% title(['Delta ori. angle, n=' num2str(sum(PertResp_units_V))])


[h p] = kstest2(histcounts(Distr_uni,-90:5:90)/length(Distr_uni(:)) , histcounts(Ori_diff,-90:5:90)/length(Distr1-Distr2))

p_1 = cumsum(histcounts(Distr_uni,-180:36:180))/(length(Distr_uni(:)));
p_2 = cumsum(histcounts(Distr1-Distr2,0:36:180))/length(Distr1-Distr2);
figure
plot(p_1,'k')
hold on
plot(p_2,'r')

% make 2D histogram out of the scatter plot
[N, Xedges, Yedges] = histcounts2(Distr1,Distr2,0:10:180,0:10:180);
scal_f = 10; sigma = 3;
I = N/sum(N(:));
clear J
J = imresize(I,scal_f);
J = imgaussfilt(J,sigma);
%J = (J/max(J(:)));
I_sh = mean(N_sh,3)/sum(sum(mean(N_sh,3)));
J_sh = imresize(I_sh,scal_f);
J_sh = imgaussfilt(J_sh,sigma);
%J_sh = (J_sh/max(J_sh(:)));
figure
imagesc(J_sh)
set(gca, 'YDir','normal','TickDir','out','box','off')
ylabel('Grating Pref Ori.');
xlabel('Perturbation Pref Ori.');
title('Shuffled')
title('Observed')
title('Observed-shuffled')
% figure
% histfit(TuningPropsV.P_dir(PertResp_units)-TuningProps.P_dir(PertResp_units),12) %-360:45:360)

% figure
% subplot(1,3,1)
% histogram(TuningPropsV.P_ori(PertResp_units_V)-TuningProps.P_ori(PertResp_units_V),-180:36:180,'EdgeColor','k','FaceColor',[229 44 37]/255,'Normalization','probability')
% xlim([-180 180])
% set(gca,'xTick',-180:45:180,'TickDir','out')
% box off;
% xlabel('Orientation'); ylabel('Units')
% set(gcf,'Position',[10 10 800 700])
% title(['Observed, n=' num2str(sum(PertResp_units_V))])
% ylim([0 0.3])
% subplot(1,3,2)
% histogram(Distr_uni(:),-180:36:180,'EdgeColor','k','FaceColor','k','Normalization','probability')
% xlim([-180 180])
% set(gca,'xTick',-180:45:180,'TickDir','out')
% box off;
% xlabel('Orientation'); ylabel('Units')
% set(gcf,'Position',[10 10 800 700])
% title(['Shuffled, n=' num2str(length(Distr_uni(:)))])
% ylim([0 0.3])
% subplot(1,3,3)
% Bars_p(1,:) = histcounts(TuningPropsV.P_ori(PertResp_units_V)-TuningProps.P_ori(PertResp_units_V),-180:36:180)/length(TuningPropsV.P_ori(PertResp_units_V)-TuningProps.P_ori(PertResp_units_V));
% Bars_p(2,:) = histcounts(Distr_uni(:),-180:36:180)/length(Distr_uni(:));
% bar(-180+18:36:180-18,Bars_p(1,:)-Bars_p(2,:),'EdgeColor','k','FaceColor',[229 44 37]/255,'BarWidth', 1)
% xlim([-180 180])
% set(gca,'xTick',-180:45:180,'TickDir','out')
% box off;
% xlabel('Orientation'); ylabel('Units')
% set(gcf,'Position',[10 10 800 700])
% title(['Observed-Shuffled'])


%% DIRECTION
figure
Binning = -90:30:360;
p_thres = 0.01; PertResp_units = (AUC_shuffled(:,1)>0);
%set(gcf,'Position',[100 100 650 600])
subplot(4,4,[9 10 13 14])
title('Dir. tuned units (vis.)')
%plot(TuningPropsV.P_dir,TuningProps.P_dir,'.k')
%load('AUC_shuffled.mat')
%Sh_responses = AUC_shuffled(:,2:end);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PertResp_units = (AUC_shuffled(:,1)>p_pert_th);
%PertResp_units = UnitsPertPos;
PertResp_units = PertRespUnits_pos;
PertResp_units_V = PertRespUnits_pos & TuningPropsV.p_tt<p_thres;
PertResp_units = PertResp_units_V; 
% PertResp_units = 1:189; PertResp_units = logical(PertResp_units');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on
% plot(TuningProps.P_dir(PertResp_units),TuningPropsV.P_dir(PertResp_units),'^r','MarkerSize',5)
plot(TuningProps.P_dir(PertResp_units_V),TuningPropsV.P_dir(PertResp_units_V),'^r','MarkerSize',5,'markerfacecolor','r')
xlim([-90 360])
ylim([-90 360])
% plot higher angles for vis. stim. responses
%plot(TuningPropsV.P_dir(TuningPropsV.P_dir>270)-360,TuningProps.P_dir(TuningPropsV.P_dir>270),'.k')
% plot(TuningProps.P_dir(PertResp_units & (TuningPropsV.P_dir>270)),TuningPropsV.P_dir(PertResp_units & (TuningPropsV.P_dir>270))-360,'^r','MarkerSize',5)
plot(TuningProps.P_dir(PertResp_units_V & (TuningPropsV.P_dir>270)),TuningPropsV.P_dir(PertResp_units_V & (TuningPropsV.P_dir>270))-360,'^r','MarkerSize',5,'markerfacecolor','r')
% plot higher angles for pert. responses
%plot(TuningPropsV.P_dir(TuningProps.P_dir>270),TuningProps.P_dir(TuningProps.P_dir>270)-360,'.k')
% plot(TuningProps.P_dir(PertResp_units & (TuningProps.P_dir>270))-360,TuningPropsV.P_dir(PertResp_units & (TuningProps.P_dir>270)),'^r','MarkerSize',5)
plot(TuningProps.P_dir(PertResp_units_V & (TuningProps.P_dir>270))-360,TuningPropsV.P_dir(PertResp_units_V & (TuningProps.P_dir>270)),'^r','MarkerSize',5,'markerfacecolor','r')
% bottom left cornere
%plot(TuningPropsV.P_dir(TuningProps.P_dir>270 & TuningPropsV.P_dir>270)-360,TuningProps.P_dir(TuningProps.P_dir>270 & TuningPropsV.P_dir>270)-360,'.k')
% plot(TuningProps.P_dir(PertResp_units & (TuningProps.P_dir>270) & TuningPropsV.P_dir>270)-360,TuningPropsV.P_dir(PertResp_units & (TuningProps.P_dir>270)& TuningPropsV.P_dir>270)-360,'^r','MarkerSize',5)
plot(TuningProps.P_dir(PertResp_units_V & (TuningProps.P_dir>270) & TuningPropsV.P_dir>270)-360,TuningPropsV.P_dir(PertResp_units_V & (TuningProps.P_dir>270)& TuningPropsV.P_dir>270)-360,'^r','MarkerSize',5,'markerfacecolor','r')
% scatterDiagHist(TuningPropsV.P_dir(PertResp_units & TuningPropsV.p_HT2<=p_thres),...
%      TuningProps.P_dir(PertResp_units & TuningPropsV.p_HT2<=p_thres)) 
ylabel('Visual preferred direction'); xlabel('Perturbation max resp. direction'); 
% legend({['Ori. tuned vis. units n=' num2str(length(TuningPropsV.P_dir(TuningPropsV.p_HT2<=p_thres)))], ...
%         ['Pert. responsive units n=' num2str(length(TuningPropsV.P_dir(PertResp_units & TuningPropsV.p_HT2<=p_thres)))]})
set(gca,'YAxisLocation', 'left','XTick',-90:90:360,'XTickLabel',[270 0:90:360],'YTick',-90:90:360,'YTickLabel',[270 0:90:360],'TickDir','out'); box off;

subplot(4,4,[1 2 5 6])
%histogram([TuningProps.P_dir(TuningProps.P_dir>270)-360; TuningProps.P_dir],-90:22.5:360,'EdgeColor','none','FaceColor','k')
hold on
% histogram([TuningProps.P_dir(PertResp_units & TuningProps.P_dir>270)-360; TuningProps.P_dir(PertResp_units)],Binning,'EdgeColor','none','FaceColor',[229 44 37]/255)
histogram([TuningProps.P_dir(PertResp_units_V & TuningProps.P_dir>270)-360; TuningProps.P_dir(PertResp_units_V)],Binning,'EdgeColor','k','FaceColor',[229 44 37]/255)
set(gca,'YAxisLocation', 'left','XTick',0:90:360,'XTickLabel',[],'TickDir','out'); box off;
%set(gca, 'XDir','reverse')
ylabel('Units')
xlim([-90 360])

subplot(4,4,[11 12 15 16])
%histogram([TuningPropsV.P_dir(TuningPropsV.P_dir>270)-360; TuningPropsV.P_dir],-90:22.5:360,'EdgeColor','none','FaceColor','k')
hold on
% histogram([TuningPropsV.P_dir(PertResp_units & TuningPropsV.P_dir>270)-360; TuningPropsV.P_dir(PertResp_units)],Binning,'EdgeColor','none','FaceColor',[229 44 37]/255)
histogram([TuningPropsV.P_dir(PertResp_units_V & TuningPropsV.P_dir>270)-360; TuningPropsV.P_dir(PertResp_units_V)],Binning,'EdgeColor','k','FaceColor',[229 44 37]/255)
xlim([-90 360])
set(gca,'XTick',-90:90:360,'XTickLabel',[],'TickDir','out','XAxisLocation', 'bottom'); box off;
ylabel('Units');
set(gca, 'XDir','reverse')
view(90,90)

% chance level
Distr1 = TuningPropsV.P_dir(PertResp_units_V);
Distr2 = TuningProps.P_dir(PertResp_units_V);
clear Distr_uni
for sh = 1:1000
    Distr1_sh = Distr1(randperm(length(Distr1)));
    Distr2_sh = Distr2(randperm(length(Distr2)));
    Distr_uni(sh,:) = Distr1_sh-Distr2_sh;
end
subplot(4,4,[3 4 7 8])
hold off
Dir_diff = rad2deg(angdiff(deg2rad(TuningPropsV.P_dir(PertResp_units_V)),deg2rad(TuningProps.P_dir(PertResp_units_V))));
histogram(Dir_diff,-180:22.5:180,'EdgeColor','r','FaceColor',[229 44 37]/255,'Normalization','probability')
hold on
histogram(Distr_uni(:),-180:22.5:180,'EdgeColor','k','FaceColor',[0.5 0.5 0.5],'Normalization','probability')
% pd = fitdist(TuningPropsV.P_dir(PertResp_units)-TuningProps.P_dir(PertResp_units),'Normal')
% hold on
% x_values = -360:1.8:360;
% y_values = (1/(pd.sigma*sqrt(2*pi))).*exp((-(x_values-pd.mean).^2)./(2*pd.sigma^2));
% plot(x_values,y_values,'LineWidth',2)
xlim([-180 180])
set(gca,'xTick',-180:45:180,'TickDir','out','View',[45,90])
box off;
xlabel('Direction'); ylabel('Units')
set(gcf,'Position',[10 10 800 700])


% chance level
% figure
% histogram(Distr_uni(:),-360:45:360,'EdgeColor','k','FaceColor',[0.5 0.5 0.5],'Normalization','probability')
% hold on
% histogram(TuningPropsV.P_dir(PertResp_units_V)-TuningProps.P_dir(PertResp_units_V),-360:45:360,'EdgeColor','k','FaceColor',[229 44 37]/255,'Normalization','probability')
% xlim([-360 360])
% set(gca,'xTick',-360:90:360,'TickDir','out','View',[45,90])
% box off;
% xlabel('Direction'); ylabel('Units')
% set(gcf,'Position',[10 10 800 700])
% title(['Delta dir. angle, n=' num2str(sum(PertResp_units_V))])

[p h] = kstest2(histcounts(Distr_uni,-360:45:360)/length(Distr_uni(:)) , histcounts(Distr1-Distr2,-360:45:360)/length(Distr1-Distr2))
p_1 = cumsum(histcounts(Distr_uni,-360:45:360))/(length(Distr_uni(:)));
p_2 = cumsum(histcounts(Distr1-Distr2,-360:45:360))/length(Distr1-Distr2);
figure
plot(p_1,'k')
hold on
plot(p_2,'r')


% (Z, smth_win, M, bins1, bins2, Fcircular1, Fcircular2,FnormXonly)
% PertResp_units = 1:189;
% Z(:,1) = histcounts(TuningProps.P_dir(PertResp_units),0:90:360);
% Z(:,2) = histcounts(TuningPropsV.P_dir(PertResp_units),0:90:360);
% smoothhist2D_corrected2_MM(Z,[2 2],[length(0:90:360) length(0:90:360)],0:90:360,0:90:360,1,1,1)
% 
% figure
% histfit(TuningPropsV.P_dir(PertResp_units)-TuningProps.P_dir(PertResp_units),12) %-360:45:360)
% 
% figure
% histfit(TuningPropsV.P_dir(PertResp_units))