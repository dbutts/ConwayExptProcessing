function CCs = generate_color_cloud_binary( NT, NX, CLRscale, SPscale, THRESH )
%
% Usage: CCs = generate_color_cloud_binary( NT, NX, CLRscale, SPscale, THRESH )

RGBscale = 0.70;

if nargin < 5
	THRESH = -0.2;
end

%% Generate raw components
SPclouds = mk_spatialcloud( NX, NX, NT, SPscale );

CC_DKL=mk_spatialcloudRGB( NX, NX, NT, CLRscale);
nrm = sqrt(sum(CC_DKL.^2, 4));
nrm(nrm==0) = 1e-6;

CCclr = CC_DKL./repmat(nrm,[1,1,1,3]) * RGBscale;

mskB = SPclouds;
mskB(mskB > THRESH) = 1;
mskB(mskB < THRESH) = 0;

% Apply mask
CCs = CCclr .* repmat(mskB,[1,1,1,3]);