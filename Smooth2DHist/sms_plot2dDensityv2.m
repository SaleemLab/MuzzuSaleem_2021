function sms_plot2dDensityv2(animal, iseries, es, obj, samples, criteria, cells, pred_cv, smthSpikeCount, FnormXonly, save_flag)
if nargin<2
    criteria='full';
end

switch criteria
    case 'full'
        cellID=[1:length(es.spikeIDs)];
        
    case 'thres'
        cellID=cells.thres_based;
        
    case 'zero'
        cellID=cells.zero_based;     
end
%% Plot 2D density (new)
cd 'C:\Users\m.morimoto\Documents\GitHub\VisNavGUI';

[X1, ~] = normalise1var(es.traj, obj.numBins);
Ndecwin=round([0.25/(1/60)]);% in samples decoding time window

figure;
conds=[samples.no, samples.swap, samples.skip1, samples.skip2];

for iconds=1:4
    
    subplot(2,2,iconds)
    
%         cd 'C:\Users\m.morimoto\Documents\GitHub\Saleem-Lab-Code\Decoding' ;
        [pred_new, X_new, Posterior, nPosterior, nonNormPosterior] = obj.predictBayesDecoder(X1(conds(:,iconds)),smthSpikeCount(conds(:,iconds),cellID),(1/60)*ones(size(X1(conds(:,iconds))))*Ndecwin , 'mean');% T=60*ones(size(X))
%         if iconds==1 %for no manipulation trials use pred from training
%             Y=pred_cv;
%         else
%             Y=pred_new';
%         end
        X=es.traj(conds(:,iconds));
        Y=pred_new;
        Z=[X(1:length(Y))/2,Y];%4
        smth=[4 4]; % need two lambdas for VisNavGUI version
        nbins=[100 100];%[50,50];
        bins1=1:100;%50;
        bins2=1:100;%50;
        Fcircular1=[];
        Fcircular2=[];
        
%         cd 'C:\Users\m.morimoto\Documents\GitHub\VisNavGUI\scripts';
        [nF, F, c1, c2, con, H] = smoothhist2D_corrected2_MM(Z, smth, nbins, bins1, bins2, Fcircular1, Fcircular2, FnormXonly);
        
        condsF(:,:,iconds)=log(nF);
        imagesc(log(nF), [-max(max(log(nF)))*0.8 max(max(log(nF)))*0.8]);%[0 5]
        axis xy equal tight
        xlabel('actual position')
        ylabel('decoded position')
        
     switch iconds
         case 1
            title([animal '_' iseries ' no manipulation ' criteria],'Interpreter','none')
         case 2
            title('swap L2 and L3','Interpreter','none')
         case 3
            title('skip L2','Interpreter','none')
         case 4
            title('skip L3','Interpreter','none')
     end
        RedWhiteBlue;
        colorbar
        set(gca, 'box','off','TickDir','out','fontsize',14,'color','none'); 
end

%% Plot 2D density (old)
% figure;
% 
% conds=[samples.no, samples.swap, samples.skip1, samples.skip2];
% 
% for iconds=1:4
%     
%     subplot(2,2,iconds)
%     
%         cd 'C:\Users\m.morimoto\Documents\GitHub\Saleem-Lab-Code\Decoding' ;
%         [pred_new, Posterior, nPosterior, nonNormPosterior] = obj.predictBayesDecoder(smthSpikeCount(conds(:,iconds),cellID), 0, 'mean');
%         
%         X=es.traj(conds(:,iconds));
% %         if iconds==1 %for no manipulation trials use pred from training
% %             Y=pred_cv;
% %         else
% %             Y=pred_new';
% %         end
%         Y=pred_new';
%         Z=[X(1:length(Y))/4,Y];
%         smth=[4 4]; % need two lambdas for VisNavGUI version
%         nbins=[50,50];
%         bins1=1:50;
%         bins2=1:50;
%         Fcircular1=[];
%         Fcircular2=[];
%         
%         cd 'C:\Users\m.morimoto\Documents\GitHub\VisNavGUI\scripts';
%         [nF, F, c1, c2, con, H] = smoothhist2D_corrected2_MM(Z, smth, nbins, bins1, bins2, Fcircular1, Fcircular2, FnormXonly);
%         
%         condsF(:,:,iconds)=log(nF);
%         imagesc(log(nF), [-max(max(log(nF)))*0.8 max(max(log(nF)))*0.8]);%[0 5]
%         axis xy equal tight
%         xlabel('actual position')
%         ylabel('decoded position')
%         
%      switch iconds
%          case 1
%             title([animal '_' iseries ' no manipulation ' criteria],'Interpreter','none')
%          case 2
%             title('swap L2 and L3','Interpreter','none')
%          case 3
%             title('skip L2','Interpreter','none')
%          case 4
%             title('skip L3','Interpreter','none')
%      end
%         RedWhiteBlue;
%         colorbar
%         set(gca, 'box','off','TickDir','out','fontsize',14,'color','none'); 
% end

%% Plot difference of logs
figure;
for iconds=1:4
        subplot(2,2,iconds)
        if iconds==1
            f=condsF(:,:,iconds);
            imagesc(f, [-max(max(f))*0.8 max(max(f))*0.8]);
        else
            f=condsF(:,:,iconds)-condsF(:,:,1);
            imagesc(f, [-max(max(f))*0.8 max(max(f))*0.8]);
        end
        
        axis xy equal tight
        xlabel('actual position')
        ylabel('decoded position')
        
        switch iconds
         case 1
            title([animal '_' iseries ' no manipulation ' criteria],'Interpreter','none')
         case 2
            title({'swap L2 and L3','difference'}, 'Interpreter','none')
         case 3
            title({'skip L2','difference'},'Interpreter','none')
         case 4
            title({'skip L3','difference'},'Interpreter','none')
        end
        RedWhiteBlue;
        colorbar
        set(gca, 'box','off','TickDir','out','fontsize',14,'color','none'); 
end

%% Plot the normalization matrix - no manipulation
figure;
%   cd 'C:\Users\m.morimoto\Documents\GitHub\Saleem-Lab-Code\Decoding' ;
%     [pred, Posterior, nPosterior, nonNormPosterior] = obj.predictBayesDecoder(smthSpikeCount(samples.no,cellID), 0, 'mean');
    [pred, X_new, Posterior, nPosterior, nonNormPosterior] = obj.predictBayesDecoder(X1(samples.no),smthSpikeCount(samples.no,cellID), 1/60 , 'mean');% T=60*ones(size(X))

    X=es.traj(samples.no);
    Y=pred;
    Z=[X(1:length(Y))/2,Y];%4
    smth=[4 4]; % need two lambdas for VisNavGUI version
    nbins=[100,100];%[50,50]
    bins1=1:100;%1:50
    bins2=1:100;%1:50
    Fcircular1=[];
    Fcircular2=[];

    [nF, F, c1, c2, con, H] = smoothhist2D_corrected2_MM(Z, smth, nbins, bins1, bins2, Fcircular1, Fcircular2, FnormXonly);

    animal=animal{1};
    iseries=iseries{1};
    
subplot(2,2,1)
    imagesc(H, [min(min(H)) max(max(H))*0.3] ); colorbar
    axis xy equal tight
    xlabel('actual position')
    ylabel('decoded position')
    title({[animal '_' iseries ' no manipulation ' criteria],'probability raw',},'Interpreter','none')
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none'); 
    
subplot(2,2,2)
    imagesc(F); colorbar
    axis xy equal tight
    xlabel('actual position')
    ylabel('decoded position')
    title('probability smoothed')
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none'); 
    
subplot(2,2,3)
    imagesc(con); colorbar
    axis xy equal tight
    xlabel('actual position')
    ylabel('decoded position')
    title('chance=sqrt(2dX*2dY)')
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none'); 

subplot(2,2,4)
    imagesc(log(nF), [-max(max(log(nF)))*0.8 max(max(log(nF)))*0.8]); colorbar
    axis xy equal tight
    xlabel('actual position')
    ylabel('decoded position')
    title('probability smoothed/chance')
    RedWhiteBlue;
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none'); 


%% old
%         cd 'X:\Archive - saleemlab\Code\general'; % old version
%         [F, c1, c2, con, H] = smoothhist2D_corrected(Z, smth, nbins, bins);
%         cd 'C:\Users\m.morimoto\Documents\GitHub\VisNavGUI\scripts';
%         [F, c1, c2, con, H] = smoothhist2D_corrected(Z, smth, nbins, bins1, bins2, Fcircular1, Fcircular2, FnormXonly);
%     imagesc(F, [min(min(F)) max(max(F))*0.8]);%[0 5]

end