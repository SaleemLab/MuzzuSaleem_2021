
figure
p_thres = 0.01; PertResp_units = (AUC_shuffled(:,1)>0);
set(gcf,'Position',[100 100 650 600])
subplot(4,4,[9 10 13 14])
title('Ori. tuned units (vis.)')
plot(TuningPropsV.P_ori,TuningProps.P_ori,'.k')
%load('AUC_shuffled.mat')
%Sh_responses = AUC_shuffled(:,2:end);
p_pert_th = prctile(Sh_responses(:),99);
%PertResp_units = (AUC_shuffled(:,1)>p_pert_th);
PertResp_units = UnitsPertPos;
hold on
plot(TuningPropsV.P_ori(PertResp_units),TuningProps.P_ori(PertResp_units),'or','MarkerSize',5)
xlim([-45 180])
ylim([-45 180])
% plot higher angles for vis. stim. responses
plot(TuningPropsV.P_ori(TuningPropsV.P_ori>135)-180,TuningProps.P_ori(TuningPropsV.P_ori>135),'.k')
plot(TuningPropsV.P_ori(PertResp_units & (TuningPropsV.P_ori>135))-180,TuningProps.P_ori(PertResp_units & (TuningPropsV.P_ori>135)),'or','MarkerSize',5)
% plot higher angles for pert. responses
plot(TuningPropsV.P_ori(TuningProps.P_ori>135),TuningProps.P_ori(TuningProps.P_ori>135)-180,'.k')
plot(TuningPropsV.P_ori(PertResp_units & (TuningProps.P_ori>135)),TuningProps.P_ori(PertResp_units & (TuningProps.P_ori>135))-180,'or','MarkerSize',5)
% bottom left cornere
plot(TuningPropsV.P_ori(TuningProps.P_ori>135 & TuningPropsV.P_ori>135)-180,TuningProps.P_ori(TuningProps.P_ori>135 & TuningPropsV.P_ori>135)-180,'.k')
plot(TuningPropsV.P_ori(PertResp_units & (TuningProps.P_ori>135)& TuningPropsV.P_ori>135)-180,TuningProps.P_ori(PertResp_units & (TuningProps.P_ori>135) & TuningPropsV.P_ori>135)-180,'or','MarkerSize',5)
% scatterDiagHist(TuningPropsV.P_ori(PertResp_units & TuningPropsV.p_HT2<=p_thres),...
%      TuningProps.P_ori(PertResp_units & TuningPropsV.p_HT2<=p_thres)) 
xlabel('Visual preferred orientation'); ylabel('Perturbation max resp. orientation'); 
% legend({['Ori. tuned vis. units n=' num2str(length(TuningPropsV.P_ori(TuningPropsV.p_HT2<=p_thres)))], ...
%         ['Pert. responsive units n=' num2str(length(TuningPropsV.P_ori(PertResp_units & TuningPropsV.p_HT2<=p_thres)))]})
set(gca,'YAxisLocation', 'left','XTick',-45:45:180,'XTickLabel',[135 0:45:180],'YTick',-45:45:180,'YTickLabel',[135 0:45:180],'TickDir','out'); box off;

subplot(4,4,[1 2 5 6])
histogram([TuningPropsV.P_ori(TuningPropsV.P_ori>135)-180; TuningPropsV.P_ori],-45:11.25:180,'EdgeColor','none','FaceColor','k')
hold on
histogram([TuningPropsV.P_ori(PertResp_units & TuningPropsV.P_ori>135)-180; TuningPropsV.P_ori(PertResp_units)],-45:11.25:180,'EdgeColor','none','FaceColor',[229 44 37]/255)
set(gca,'YAxisLocation','left','XTick',0:45:180,'XTickLabel',[],'TickDir','out'); box off;
%set(gca, 'XDir','reverse')
ylabel('Units')
xlim([-45 180])

subplot(4,4,[11 12 15 16])
histogram([TuningProps.P_ori(TuningProps.P_ori>135)-180; TuningProps.P_ori],-45:11.25:180,'EdgeColor','none','FaceColor','k')
hold on
histogram([TuningProps.P_ori(PertResp_units & TuningProps.P_ori>135)-180; TuningProps.P_ori(PertResp_units)],-45:11.25:180,'EdgeColor','none','FaceColor',[229 44 37]/255)
xlim([-45 180])
set(gca,'XTick',-45:45:180,'XTickLabel',[],'TickDir','out','XAxisLocation', 'bottom'); box off;
ylabel('Units');
set(gca, 'XDir','reverse')
view(90,90)


subplot(4,4,[3 4 7 8])
scatterDiagHist(TuningPropsV.P_ori(PertResp_units),TuningProps.P_ori(PertResp_units))  

figure
histogram(TuningPropsV.P_ori(PertResp_units)-TuningProps.P_ori(PertResp_units),-180:18:180)

figure
histogram(TuningPropsV.P_dir(PertResp_units)-TuningProps.P_dir(PertResp_units),-360:45:360)

%% DIRECTION

figure
p_thres = 0.01; PertResp_units = (AUC_shuffled(:,1)>0);
set(gcf,'Position',[100 100 650 600])
subplot(4,4,[9 10 13 14])
title('Ori. tuned units (vis.)')
plot(TuningPropsV.P_dir,TuningProps.P_dir,'.k')
%load('AUC_shuffled.mat')
%Sh_responses = AUC_shuffled(:,2:end);
p_pert_th = prctile(Sh_responses(:),99);
%PertResp_units = (AUC_shuffled(:,1)>p_pert_th);
PertResp_units = UnitsPertPos;
hold on
plot(TuningPropsV.P_dir(PertResp_units),TuningProps.P_dir(PertResp_units),'^r','MarkerSize',5)
xlim([-90 360])
ylim([-90 360])
% plot higher angles for vis. stim. responses
plot(TuningPropsV.P_dir(TuningPropsV.P_dir>270)-360,TuningProps.P_dir(TuningPropsV.P_dir>270),'.k')
plot(TuningPropsV.P_dir(PertResp_units & (TuningPropsV.P_dir>270))-360,TuningProps.P_dir(PertResp_units & (TuningPropsV.P_dir>270)),'^r','MarkerSize',5)
% plot higher angles for pert. responses
plot(TuningPropsV.P_dir(TuningProps.P_dir>270),TuningProps.P_dir(TuningProps.P_dir>270)-360,'.k')
plot(TuningPropsV.P_dir(PertResp_units & (TuningProps.P_dir>270)),TuningProps.P_dir(PertResp_units & (TuningProps.P_dir>270))-360,'^r','MarkerSize',5)
% bottom left cornere
plot(TuningPropsV.P_dir(TuningProps.P_dir>270 & TuningPropsV.P_dir>270)-360,TuningProps.P_dir(TuningProps.P_dir>270 & TuningPropsV.P_dir>270)-360,'.k')
plot(TuningPropsV.P_dir(PertResp_units & (TuningProps.P_dir>270)& TuningPropsV.P_dir>270)-360,TuningProps.P_dir(PertResp_units & (TuningProps.P_dir>270) & TuningPropsV.P_dir>270)-360,'^r','MarkerSize',5)
% scatterDiagHist(TuningPropsV.P_dir(PertResp_units & TuningPropsV.p_HT2<=p_thres),...
%      TuningProps.P_dir(PertResp_units & TuningPropsV.p_HT2<=p_thres)) 
xlabel('Visual preferred orientation'); ylabel('Perturbation max resp. orientation'); 
% legend({['Ori. tuned vis. units n=' num2str(length(TuningPropsV.P_dir(TuningPropsV.p_HT2<=p_thres)))], ...
%         ['Pert. responsive units n=' num2str(length(TuningPropsV.P_dir(PertResp_units & TuningPropsV.p_HT2<=p_thres)))]})
set(gca,'YAxisLocation', 'left','XTick',-90:90:360,'XTickLabel',[270 0:90:360],'YTick',-90:90:360,'YTickLabel',[270 0:90:360],'TickDir','out'); box off;

subplot(4,4,[1 2 5 6])
histogram([TuningPropsV.P_dir(TuningPropsV.P_dir>270)-360; TuningPropsV.P_dir],-90:22.5:360,'EdgeColor','none','FaceColor','k')
hold on
histogram([TuningPropsV.P_dir(PertResp_units & TuningPropsV.P_dir>270)-360; TuningPropsV.P_dir(PertResp_units)],-90:22.5:360,'EdgeColor','none','FaceColor',[229 44 37]/255)
set(gca,'YAxisLocation', 'left','XTick',0:90:360,'XTickLabel',[],'TickDir','out'); box off;
%set(gca, 'XDir','reverse')
ylabel('Units')
xlim([-90 360])

subplot(4,4,[11 12 15 16])
histogram([TuningProps.P_dir(TuningProps.P_dir>270)-360; TuningProps.P_dir],-90:22.5:360,'EdgeColor','none','FaceColor','k')
hold on
histogram([TuningProps.P_dir(PertResp_units & TuningProps.P_dir>270)-360; TuningProps.P_dir(PertResp_units)],-90:22.5:360,'EdgeColor','none','FaceColor',[229 44 37]/255)
xlim([-90 360])
set(gca,'XTick',-90:90:360,'XTickLabel',[],'TickDir','out','XAxisLocation', 'bottom'); box off;
ylabel('Units');
set(gca, 'XDir','reverse')
view(90,90)


subplot(4,4,[3 4 7 8])
scatterDiagHist(TuningPropsV.P_dir(PertResp_units),TuningProps.P_dir(PertResp_units))  

