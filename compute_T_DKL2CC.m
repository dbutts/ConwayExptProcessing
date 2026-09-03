function T_DKL2CC = compute_T_DKL2CC(varargin)

% Same arguments as for compute_T_RGB2LMS and compute_T_DKL2RGB

T_RGB2LMS = compute_T_RGB2LMS(varargin{:});
T_DKL2RGB = compute_T_DKL2RGB(varargin{:});

bkg_RGB = [0.5; 0.5; 0.5];
bkg_coneExcitations = T_RGB2LMS * bkg_RGB;

T_DKL2CC = 0.5.*(T_RGB2LMS * T_DKL2RGB)./bkg_coneExcitations;

end
