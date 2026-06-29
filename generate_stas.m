function [STA, STA_CC] = generate_stas(Robs, stim, nlags, plotflag, saveflag, v_bkg, varargin)

% Robs is expected to be a # units x # frames matrix, while stim is
% expected to be Y pixels x X pixels x 3 x # frames

if nargin < 5
    saveflag = 0;
end
if nargin < 4
    plotflag = 0;
end

nLags = 10;
tempSTA = [];
S = transpose(single(reshape(stim,prod(size(stim,1:3)), [])))./127;

tic;
fprintf('Computing STAs\n');

for lag = 0:nLags-1
    tempSTA(:,:,:,:,lag+1) = (Robs(:,lag+1:end) * S(1:end-lag,:))./sum(Robs(:,lag+1:end),2);
end

STA = reshape(tempSTA, size(tempSTA,1), size(stim,1), size(stim,2), 3, nLags);


% convert STA from DKL to LMS
STA_CC =[];
if nargin > 5
    T_RGB2LMS = compute_T_RGB2LMS(varargin{:});
    T_DKL2RGB = compute_T_DKL2RGB(varargin{:});

    bkg_LMS = T_RGB2LMS*(v_bkg.*(T_DKL2RGB*[0;0;0])+v_bkg); % v_bkg = [0.5;0.5;0.5]

   % for each STA and each lag, convert from DKL to RGB

    for unit = 1:size(STA,1)
        for lag = 1:nlags
            thisSTA_DKL = squeeze(STA(unit,:,:,:,lag));
            thisSTA_DKL = transpose(reshape(thisSTA_DKL, [],3));
            thisSTA_LMS = T_RGB2LMS * (v_bkg .* (T_DKL2RGB * thisSTA_DKL) + v_bkg);
            thisSTA_CC = (thisSTA_LMS - bkg_LMS)./bkg_LMS;
            thisSTA_CC = reshape(thisSTA_CC, 3, size(stim,1),size(stim,2));
            thisSTA_CC = permute(thisSTA_CC, [2 3 1]);
            STA_CC(unit,:,:,:,lag) = thisSTA_CC;

        end
    end



end