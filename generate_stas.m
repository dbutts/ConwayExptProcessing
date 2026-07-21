function STA = generate_stas(Robs, stim, nlags, v_bkg, varargin)

% Robs is expected to be a # units x # frames matrix, while stim is
% expected to be Y pixels x X pixels x 3 x # frames

% STA are # units x Y pixels x X pixels x 3 chromatic dimensions x nlags  

% Need psychtoolbox on path
%addpath(genpath('~/Psychtoolbox-3/Psychtoolbox/'))

% Example usage
% rx   = 0.6280;
% ry   = 0.3310;
% gx   = 0.3059;
% gy   = 0.5826;
% bx   = 0.1639;
% by   = 0.0617;

% R_maxLumCdm2 = 18.6469;
% G_maxLumCdm2 = 75.8449;
% B_maxLumCdm2 = 10.5313;

% xyMat = [rx ry; gx gy; bx by] 
% lumMat = [R_maxLumCdm2;
% G_maxLumCdm2; B_maxLumCdm2]

% nlags = 10
% v_bkg = [0.5; 0.5; 0.5];

% [STA, STA_CC] = generate_stas(data.Robs, data.stim, nlags, v_bkg, 'xyMat', xyMat, 'lumMat', lumMat)

% if nargin < 5
%     saveflag = 0;
% end
% if nargin < 4
%     plotflag = 0;
% end

%nLags = 10;
tempSTA = [];
S = transpose(single(reshape(stim,prod(size(stim,1:3)), [])))./127;

tic;
fprintf('Computing STAs\n');

for lag = 0:nlags-1
    tempSTA(:,:,:,:,lag+1) = (Robs(:,lag+1:end) * S(1:end-lag,:))./sum(Robs(:,lag+1:end),2);
end

STA.DKL = reshape(tempSTA, size(tempSTA,1), size(stim,1), size(stim,2), 3, nlags);


% convert STA from DKL to LMS

if nargin > 5
    T_RGB2LMS = compute_T_RGB2LMS(varargin{:});
    T_DKL2RGB = compute_T_DKL2RGB(varargin{:});
    bkg_LMS = T_RGB2LMS*(v_bkg.*(T_DKL2RGB*[0;0;0])+v_bkg); % v_bkg = [0.5;0.5;0.5]

   % for each STA and each lag, convert from DKL to RGB

    for unit = 1:size(STA.DKL,1)
        for lag = 1:nlags
            thisSTA_DKL = squeeze(STA.DKL(unit,:,:,:,lag));
            thisSTA_DKL = transpose(reshape(thisSTA_DKL, [],3));

            thisSTA_RGB = v_bkg.*(T_DKL2RGB * thisSTA_DKL)+v_bkg;
            thisSTA_RGB = reshape(thisSTA_RGB, 3, size(stim,1),size(stim,2));
            thisSTA_RGB = permute(thisSTA_RGB, [2 3 1]);

            STA.RGB(unit,:,:,:,lag) = thisSTA_RGB;

            thisSTA_LMS = T_RGB2LMS * (v_bkg .* (T_DKL2RGB * thisSTA_DKL) + v_bkg);

            thisSTA_CC = (thisSTA_LMS - bkg_LMS)./bkg_LMS;
            thisSTA_CC = reshape(thisSTA_CC, 3, size(stim,1),size(stim,2));
            thisSTA_CC = permute(thisSTA_CC, [2 3 1]);

            STA.coneContrast(unit,:,:,:,lag) = thisSTA_CC;

            thisSTA_LMS = reshape(thisSTA_LMS, 3, size(stim,1), size(stim,2));
            thisSTA_LMS = permute(thisSTA_LMS, [2 3 1]);

            STA.LMS(unit,:,:,:,lag) = thisSTA_LMS;

        end
    end


end