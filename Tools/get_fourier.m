
lag = 5;
cell = 3;
% pull the STA data from the current (cell, lag, pixel)
cur_STA = stas(cell,lag,:);
% Reshape data to be in form (pixels, pixels, (Lum, L-M, S)) (60,60,3)
cloud = reshape(cur_STA, [60, 60, 3]);
% Find LUM data
STA_LUM = cloud(:,:,1);
% Find L-M data
STA_LM = cloud(:,:,2);
% Find S data
STA_S = cloud(:,:,3);
% Find L and M data
cur_StA2=reshape(cur_STA(1,:),60*60,3);
% Find L data
STA_L = cur_StA2*[-0.040, 0.5220,0]'; % from 01_2022 calib
STA_L=reshape(STA_L,60,60);
% Find M data
StA_M = cur_StA2*[0.0351, -0.4625, 0]'; % from 01_2022 calib
StA_M=reshape(StA_M,60,60);
% plot each 
figure
subplot(2, 3, 1)
imagesc(STA_LUM)
title('Luminance')
subplot(2, 3, 2)
imagesc(STA_LM)
title('L-M')
subplot(2, 3, 3)
imagesc(STA_S)
title('S')
subplot(2, 3, 4)
imagesc(STA_L)
title('L')
subplot(2, 3, 5)
imagesc(StA_M)
title('M')
sgtitle(['STAs Lag: ', num2str(lag), ' SU # ', num2str(cell)])

% 1) find fourier transform 2) shift it so that 0 frequency is in middle 3)
% apply the absolute value so that there are no complex numbers
STA_LUM_F = abs(fftshift(fft2(STA_LUM)));
STA_LM_F = abs(fftshift(fft2(STA_LM)));
STA_S_F = abs(fftshift(fft2(STA_S)));
STA_L_F = abs(fftshift(fft2(STA_L)));
STA_M_F = abs(fftshift(fft2(StA_M)));

% Plot each fourier transform
figure
subplot(2, 3, 1)
imagesc(STA_LUM_F)
title('Luminance')
subplot(2, 3, 2)
imagesc(STA_LM_F)
title('L-M')
subplot(2, 3, 3)
imagesc(STA_S_F)
title('S')
subplot(2, 3, 4)
imagesc(STA_L_F)
title('L')
subplot(2, 3, 5)
imagesc(STA_M_F)
title('M')
sgtitle(['Fourier Transformed STAs: ', num2str(lag), ' SU # ', num2str(cell)])


