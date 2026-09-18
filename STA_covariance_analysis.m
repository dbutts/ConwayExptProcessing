close all;
nLags = size(STA.DKL, 5);
nUnits = size(STA.DKL,1);
nx = size(STA.DKL,3);
ny = size(STA.DKL,2);
nchrom_chans = size(STA.DKL,4);

unit = 19;

thisSTA = squeeze(STA.DKL(unit,:,:,:,:));%

%thisSTA = circshift(thisSTA, [30 30]);


C_spatial = cov(transpose(reshape(thisSTA, [], nLags)));

figure;
imagesc(C_spatial); colormap gray
title('Spatial covariance matrix', 'FontSize', 16)


C_temporal = cov(reshape(permute(thisSTA, [1 2 4 3]), [], nchrom_chans * nLags)); colormap gray

figure;
imagesc(C_temporal); colormap gray
title('Temporal covariance matrix', 'FontSize', 16)


C_temporal_diag = diag(C_temporal);

colors = {'k', 'r', 'b'};
figure; hold on
for i = 1:nchrom_chans
    plot(C_temporal_diag((nLags*(i-1) + 1): nLags*i ), colors{i}, 'linewidth',2) 
end

% filtering based on C_spatial

D = C_spatial .*( abs(C_spatial - mean(C_spatial(:))) > 3 * std(C_spatial(:)));

D_sum = sum(abs(D),1);

D_sum = D_sum./max(D_sum(:));

mask = reshape(D_sum, ny,nx,3);
mask = repmat(mask, [1 1 1 nLags]);
thisSTA_masked =thisSTA .* mask;


figure; hold on
for i = 1:nchrom_chans
    subplot(1,nchrom_chans,i)
    imagesc(squeeze(thisSTA(:,:,i,5))); colormap gray; axis square

end

figure; hold on
for i = 1:nchrom_chans
    subplot(1,nchrom_chans,i)
    imagesc(squeeze(thisSTA_masked(:,:,i,5))); colormap gray; axis square

end

%  blocks = reshape(C_spatial, 60, 180, 60, 180);
% blocks = permute(blocks, [1 3 2 4]);
% blocks = reshape(blocks, 60, 60, []);


for i = 1:nchrom_chans
    idx = (nx*ny)*(i-1) +1 : (nx*ny)*i;

     blocks = reshape(C_spatial(idx,idx), 60, 60, 60, 60);
blocks = permute(blocks, [1 3 2 4]);
blocks = reshape(blocks, 60, 60, []);


    
    M{i} = squeeze(mean(abs(blocks),3));
figure, imagesc(M{i});
end

%% fft

unit = 1;
lag  =5;
F = fftshift(fftshift(fft2(permute(STA.DKL, [2 3 1 4 5])), 1), 2);
F = permute(F, [3 1 2 4 5]);

V = squeeze(var(abs(F(unit,:,:,1,:)),[], 5));
H = V > prctile(V(:), 99);
%H = V./max(V(:));
%H = V > 0.25.*max(V(:));
I_filtered_FT = fftshift(fftshift(H).*squeeze((fftshift(F(unit,:,:,1,lag)))));

I_filtered = real(ifft2(fftshift(I_filtered_FT)));

%I_filtered = circshift(I_filtered, [30 30]);
figure, imagesc(I_filtered); colormap gray

figure, imagesc(squeeze(STA.DKL(unit,:,:,1,lag))); colormap gray