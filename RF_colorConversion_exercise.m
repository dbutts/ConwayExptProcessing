
% L isolating
L_CC = [0.15; 0;0];
L_LMS = bkg_LMS.*L_CC + bkg_LMS;
L_RGB = T_RGB2LMS \ L_LMS;

L_DKL = T_DKL2RGB \ ((L_RGB - v_bkg)./v_bkg);

% M isloating
M_CC = [0; -0.15;0];
M_LMS = bkg_LMS.*M_CC + bkg_LMS;
M_RGB = T_RGB2LMS \ M_LMS;

M_DKL = T_DKL2RGB \ ((M_RGB - v_bkg)./v_bkg);


% visualization
% ground truth RF drawing input from 1 L cone, 1 M cone

[X,Y] = meshgrid(1:60);

s = 1.5;
L_CC_RF = 0.15 .* exp(-((X-29)/s).^2) .* exp(-((Y-29)/s).^2);
M_CC_RF = -0.15 .* exp(-((X-31)/s).^2) .* exp(-((Y-31)/s).^2);
S_CC_RF = 0.* ones(size(X));

RF_CC = cat(3, L_CC_RF, M_CC_RF, S_CC_RF);

RF_CC = transpose(reshape(RF_CC, [],3)); % now 3 x (NX x NY)

RF_LMS = bkg_LMS.*RF_CC + bkg_LMS;


% transform from LMS to DKL

% first, LMS to RGB

RF_RGB = T_RGB2LMS \ RF_LMS;   

% next  RGB to DKL :D

RF_DKL = T_DKL2RGB \ ((RF_RGB - v_bkg)./v_bkg);


RF_LMS2 = T_RGB2LMS * (v_bkg.*(T_DKL2RGB * RF_DKL)+v_bkg);
RF_CC2 = (RF_LMS2 - bkg_LMS)./bkg_LMS;

RF_DKL = reshape(transpose(RF_DKL), 60,60,3);

RF_CC2 = reshape(transpose(RF_CC2), 60,60,3);