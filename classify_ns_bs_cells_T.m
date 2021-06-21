function [c_idx, waveforms, explained] = classify_ns_bs_cells_T(waveforms)


% Normalize waveforms
waveforms = (waveforms-min(waveforms))./(max(waveforms)-min(waveforms))-1;

%Align waveforms to their trough
[a, ia] = min(waveforms);
diff_col = ia - round(mean(ia));

for i=1:size(waveforms,2)
    if diff_col(i)>0
        waveforms(:,i) = [waveforms((diff_col(i)+1):end,i); zeros(1,diff_col(i))'];
    elseif diff_col(i)<0
        waveforms(:,i) = [zeros(1,abs(diff_col(i)))'; waveforms(1:(end-abs(diff_col(i))),i)];
    else
        waveforms(:,i) = waveforms(:,i);
    end
end

%Decide between two or three clusters
%pca
[coeff, score, latent, tsquared, explained] = pca(waveforms);
%kmeans clusterring
c_idx2 = kmeans(waveforms',2);
figure; gplotmatrix(coeff(:,1:5),[],c_idx2);
c_idx3 = kmeans(waveforms',3);
figure; gplotmatrix(coeff(:,1:5),[],c_idx3);
% c_idx4 = kmeans(coeff(:,1:3),2);
% figure; gplotmatrix(coeff(:,1:5),[],c_idx4);

c_idx = c_idx3;

% clust = zeros(size(coeff(:,1:5),1),6);
% for i=1:6
%     clust(:,i) = kmeans(coeff(:,1:5),i,'emptyaction','singleton','replicate',5);
% end
% 
% eva = evalclusters(waveforms',clust,'CalinskiHarabasz');
% 
% eva = evalclusters(coeff(:,1:3),'kmeans','CalinskiHarabasz','KList',[1:6])
% 
% c_idx = clust(:,eva.OptimalK);

% accept = input('Enter 2 for two clusters, 3 for three clusters: ');
% 
% if accept == 2
%     c_idx = c_idx2;
% elseif accept == 3
%     c_idx = c_idx3;
% end
