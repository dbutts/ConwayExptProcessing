
target_SUs = [5];
lag = 3;

figure
binned_SU = [single(data.Robs'), single(data.RobsMU')];
use_inds = data.valid_data;
use_inds(end-10:end) = []; %cut last few indices to avoid artifacts
stim_shift=permute(data.stim,[4 1 2 3]);
NT=size(data.ETtrace,2);
stim_deltas = zeros(2,NT); 
stim_shift = shift_stim( stim_shift, data.ETtrace, stim_deltas );
stim2 = single(reshape(stim_shift,size(stim_shift,1),3*60*60))./127;
cur_STA1(lag,:) = binned_SU(use_inds+lag,target_SUs)' * stim2(use_inds,:);
curlim=max(abs(cur_STA1(:)')); 
    if curlim==0; curlim=0.1; end % avoids plotting bugs if a bad STA is included in a large set of plots
cur_STA2=reshape(cur_STA1(lag,:),60,180);
%cur_STA2(:,[60 120])=curlim; % adds lines to separate lum, L-M, and S


% plot lum, L-M, S together
subplot(4,2,1)
imagesc(cur_STA2);  clim([-curlim curlim]); pbaspect([3 1 1]); title('Lum, L-M, S');

%plot lum portion
Lum = cur_STA2(:,1:60);
subplot(6,2,2); imagesc(Lum); clim([-curlim curlim]); pbaspect([1 1 1]); title('Lum');
subplot(6,2,3); imhist(Lum); title('Lum histogram')

%plot L-M portion
LM = cur_STA2(:,60:120);
subplot(6,2,4); imagesc(LM); clim([-curlim curlim]); pbaspect([1 1 1]); title('L-M');
subplot(6,2,5); imhist(LM); title('L-M histogram')

% plot S portion
S = cur_STA2(:,120:180);
subplot(6,2,6); imagesc(S); clim([-curlim curlim]); pbaspect([1 1 1]); title('S'); colormap('gray'); % change the color map to be grey
subplot(6,2,7); imhist(S); title('S histogram')

% plot L portion
cur_StA3=reshape(cur_STA1(curlag,:),60*60,3);
cur_StA_L = cur_StA3*[-0.040, 0.5220,0]'; % from 01_2022 calib
cur_StA_L2=reshape(cur_StA_L,60,60);
subplot(6,2,8)
imagesc(cur_StA_L2); clim([-curlim curlim].*.522); pbaspect([1 1 1]); title('L'); % color limit based on L-isolating max
subplot(6,2,9); imhist(cur_StA_L2); title('L histogram')

% Plot M portion
cur_StA_M = cur_StA3*[0.0351, -0.4625, 0]'; % from 01_2022 calib
cur_StA_M2=reshape(cur_StA_M,60,60);
subplot(6,2,10)
imagesc(cur_StA_M2); clim([-curlim curlim].*.4625); pbaspect([1 1 1]); title('M')
subplot(6,2,11); imhist(cur_StA_M2); title('M histogram')

%% Plot lum, L-M, S, L, and M all together
all = [cur_STA2 cur_StA_L2 cur_StA_M2];
subplot(4,2,7); imagesc(all); pbaspect([5 1 1]); title('lum, L-M, S, L, and M')
xlabel('x')

