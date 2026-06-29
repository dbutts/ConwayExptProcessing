function T_DKL2RGB = compute_T_DKL2RGB(varargin)

% Required inputs: R_bkgIntensity, G_bkgIntensity, B_bkgIntensity
% Option set 1: wavelengths, red_spd, green_spd, blue_spd
% Option set 2: xyMat = [rx ry; gx gy; bx by], lumMat = [R_maxLumCdm2;
% G_maxLumCdm2; B_maxLumCdm2]
% B_maxLumCdm2

% written by mjg 6/27/26

% Note: spds must be raw spds in units of W/Sr/m2
% function to be run in folder containing required PTB files:
% T_xyzJuddVos.mat, SplineCmf.m, SToWls.m

%% Parse inputs
p = inputParser;
addParameter(p, 'RGBmat', [0.5; 0.5; 0.5]);
addParameter(p, 'wavelengths', []);
addParameter(p, 'red_spd', []);
addParameter(p, 'green_spd', []);
addParameter(p, 'blue_spd', []);
addParameter(p, 'xyMat', []);
addParameter(p, 'lumMat', []);

parse(p, varargin{:});

R_bkgIntensity = p.Results.RGBmat(1);
G_bkgIntensity = p.Results.RGBmat(2);
B_bkgIntensity = p.Results.RGBmat(3);
wavelengths = p.Results.wavelengths;
red_spd = p.Results.red_spd;
green_spd = p.Results.green_spd;
blue_spd = p.Results.blue_spd;

if size(p.Results.xyMat, 1) < size(p.Results.xyMat,2)
    p.Results.xyMat = transpose(p.Results.xyMat);
end

if ~isempty(p.Results.xyMat)
    rx = p.Results.xyMat(1,1);
    ry = p.Results.xyMat(1,2);
    gx = p.Results.xyMat(2,1);
    gy = p.Results.xyMat(2,2);
    bx = p.Results.xyMat(3,1);
    by = p.Results.xyMat(3,2);
end

if ~isempty(p.Results.lumMat)
    R_maxLumCdm2 = p.Results.lumMat(1);
    G_maxLumCdm2 = p.Results.lumMat(2);
    B_maxLumCdm2 = p.Results.lumMat(3);
end


% If computing from raw spds, make sure they are column vectors
if size(wavelengths, 1) < size(wavelengths,2)
    wavelengths = transpose(wavelengths);
end

if size(red_spd,1) < size(red_spd,2)
    red_spd = transpose(red_spd);
end

if size(green_spd,1) < size(green_spd,2)
    green_spd = transpose(green_spd);
end

if size(blue_spd,1) < size(blue_spd,2)
    blue_spd = transpose(blue_spd);
end

load T_xyzJuddVos.mat

%% Do the calcuations

Wr = R_bkgIntensity;
Wg = G_bkgIntensity;
Wb = B_bkgIntensity;

P_LMS = compute_T_RGB2LMS(varargin{:});

Lr = P_LMS(1,1); Lg = P_LMS(1,2); Lb = P_LMS(1,3);
Mr = P_LMS(2,1); Mg = P_LMS(2,2); Mb = P_LMS(2,3);
Sr = P_LMS(3,1); Sg = P_LMS(3,2); Sb = P_LMS(3,3);

dBx_dGx = [Sb Sg; Lb+Mb Lg+Mg] \ [Sr*Wr; (Lr+Mr)*Wr];
dRx = -Wr; dBx = dBx_dGx(1); dGx = dBx_dGx(2);

dGy_dRy = [Lg Lr; Mg Mr] \ [Lb*Wb; Mb*Wb];
dBy = -Wb; dGy = dGy_dRy(1); dRy = dGy_dRy(2);

e_Lum = [1; 1; 1];
e_LM = [-dRx/Wr; -dGx/Wg; -dBx/Wb];
e_S = [-dRy/Wr; -dGy/Wg; -dBy/Wb];

T_DKL2RGB = [e_Lum e_LM e_S];

end

