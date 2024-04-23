%% Generate STAS
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

cur_STA=reshape(cur_STA1(lag,:),60,180);

%% Color Imhist Entropy
cur_STA2 = cur_STA; %redefine variable so we can re-run this cell without increasing size of cur_STA2

% Get L portion 
cur_StA3=reshape(cur_STA1(curlag,:),60*60,3);
cur_StA_L = cur_StA3*[-0.040, 0.5220,0]'; % from 01_2022 calib
cur_StA_L2=reshape(cur_StA_L,60,60);
% Add L to cur_STA2
cur_STA2 = [cur_STA2, cur_StA_L2];

% Get M portion
cur_StA_M = cur_StA3*[0.0351, -0.4625, 0]'; % from 01_2022 calib
cur_StA_M2=reshape(cur_StA_M,60,60);
subplot(6,2,10)
% Add M to cur_STA2
cur_STA2 = [cur_STA2, cur_StA_M2];

% Draw lines to divide each Lum, L-M ect.
curlim=max(abs(cur_STA1(:)')); 
    if curlim==0; curlim=0.1; end % avoids plotting bugs if a bad STA is included in a large set of plots
cur_STA2(:,[60 120 180 240])=curlim; % draws lines to divide Lum, L-M, and S


% rescale all of the plots at the same time so later we have 
STA2_rescale = rescale(cur_STA2); % rescale all elements to the interval [0 1]
STA2_array = reshape(STA2_rescale', 1, []);
STA2Kur = kurtosis(STA2_array);
STA2Skew = skewness(STA2_array);

% plot lum, L-M, S, L, and M together
subplot(6,2,1); imagesc(STA2_rescale);  clim([0 1]); pbaspect([5 1 1]); title('Lum, L-M, S, L, and M')

eAll = entropy(double(cur_STA2));
subplot(6,2,2); imhist(STA2_rescale); title(['All graphs with kurtosis = ' num2str(STA2Kur) ' skewness = ' num2str(STA2Skew)]);

%plot lum portion
Lum = STA2_rescale(:,1:59);
subplot(6,2,3); imagesc(Lum); clim([0 1]); pbaspect([1 1 1]); title('Lum');

Lum_array = reshape(Lum', 1, []); %change Lum into an array
LumSkew = skewness(Lum_array);
LumKur = kurtosis(Lum_array);
Lum_range = max(Lum_array) - min(Lum_array);
subplot(6,2,4); imhist(Lum); title(['Lum with skewness = ' num2str(LumSkew) ' kurtosis = ' num2str(LumKur)])

%plot L-M portion
LM = STA2_rescale(:,61:119);
subplot(6,2,5); imagesc(LM); clim([0 1]); pbaspect([1 1 1]); title('L-M');

LM_array = reshape(LM', 1, []);
LMSkew = skewness(LM_array);
LMKur = kurtosis(LM_array);
LM_range = max(LM_array) - min(LM_array);
subplot(6,2,6); imhist(LM); title(['LM with skewness = ' num2str(LMSkew) ' kurtosis = ' num2str(LMKur)])

% plot S portion
S = STA2_rescale(:,121:179);
subplot(6,2,7); imagesc(S); clim([0 1]); pbaspect([1 1 1]); title('S'); 

S_array = reshape(S', 1, []);
SSkew = skewness(S_array);
SKur = kurtosis(S_array);
S_range = max(S_array) - min(S_array);
subplot(6,2,8); imhist(S); title(['S with skewness = ' num2str(SSkew) ' kurtosis = ' num2str(SKur)])

% plot L portion
L = STA2_rescale(:, 181:239);
subplot(6,2,9); imagesc(L); clim([0 1]); pbaspect([1 1 1]); title('L'); 

L_array = reshape(L', 1, []);
LSkew = skewness(L_array);
LKur = kurtosis(L_array);
L_range = max(L_array) - min(L_array);
subplot(6,2,10); imhist(L); title(['L with skewness = ' num2str(LSkew) ' kurtosis = ' num2str(LKur)])

% plot M portion
M = STA2_rescale(:, 241:299);
subplot(6,2,11); imagesc(M); clim([0 1]); pbaspect([1 1 1]); title('M'); 

M_array = reshape(M', 1, []);
MSkew = skewness(M_array);
MKur = kurtosis(M_array);
M_range = max(M_array) - min(M_array);
subplot(6,2,12); imhist(M); title(['M with skewness = ' num2str(MSkew) ' kurtosis = ' num2str(MKur)])

colormap('gray'); % change the color map to be grey

range_total = Lum_range + LM_range + S_range; %+ L_range + M_range;
drive_Lum = Lum_range / range_total * 100;
drive_LM = LM_range / range_total * 100;
drive_S = S_range / range_total * 100;
% drive_L = L_range / range_total * 100;
% drive_M = M_range / range_total * 100;
sgtitle(['% Min/Max Drives: Lum = ' num2str(drive_Lum) ' LM = ' num2str(drive_LM) ' S = ' num2str(drive_S) '\newline .' ])
    %' L = ' num2str(drive_L) ' M = ' num2str(drive_M) '\newline .'])