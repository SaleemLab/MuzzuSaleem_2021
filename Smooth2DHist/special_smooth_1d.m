function [output] = special_smooth_1d(input, win, bins, nGrid)

output = zeros(size(input));

x = (-nGrid:nGrid)/nGrid;
Smoother = exp(-x.^2/win^2/2);
Smoother = Smoother./sum(Smoother);

if win==1;
    output(:) = nanmean(input(:));
else
%     if gpuDeviceCount>0
%         Ga = gpuArray(single(input));
%         Gb = gpuArray(single(Smoother));
%         Gc = conv(Ga, Gb, 'same');
%         output = gather(Gc);
%     else
        output = conv(input, Smoother, 'same');
%     end
end