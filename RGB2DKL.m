function T_DKL2RGB = DKL2RGB(varargin)


p = inputParser;
addParameter(p, 'wavelengths', []);
addParameter(p, 'red_spd', []);
addParameter(p, 'green_spd', []);
addParameter(p, 'blue_spd', []);
addParameter(p, 'rx', []);
addParameter(p, 'ry', []);
addParameter(p, 'gx', []);
addParameter(p, 'gy', []);
addParameter(p, 'bx', []);
addParameter(p, 'by', []);


%addpath(genpath('~/Git/Psychtoolbox')) 
% function to be run in folder containing required PTB files:
% T_xyzJuddVos.mat, SplineCmf.m, SToWls.m

% Method 1, compute T_DKL2RGB from raw spectral power distributions
red_spd = mean(getfield(load('/mnt/isilon/users/greenemj/measurements_MJG/PrimaryMeasurements_20-Nov-2025_MJG/red_spd.mat'), 'red_spd'),2);
green_spd = mean(getfield(load('/mnt/isilon/users/greenemj/measurements_MJG/PrimaryMeasurements_20-Nov-2025_MJG/green_spd.mat'), 'green_spd'),2);
blue_spd = mean(getfield(load('/mnt/isilon/users/greenemj/measurements_MJG/PrimaryMeasurements_20-Nov-2025_MJG/blue_spd.mat'), 'blue_spd'),2);
wavelengths = getfield(load('/mnt/isilon/users/greenemj/measurements_MJG/PrimaryMeasurements_20-Nov-2025_MJG/wavelengths.mat'), 'wls');

% If computing from raw spds, make sure they are column vectors
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
T_XYZ = SplineCmf(SToWls(S_xyzJuddVos), T_xyzJuddVos, wavelengths);
RGB = [red_spd green_spd blue_spd];

T_XYZ2LMS = [0.11514 0.54312 -0.03286
            -0.15514 0.45684 0.03286
             0 0 0.01608];

P_XYZ = T_XYZ * RGB;
P_LMS = T_XYZ2LMS * P_XYZ;

Wr = 0.5; Wg = 0.5; Wb = 0.5;
W = [Wr; Wg; Wb];

% Method 2 compute from xyz values of primaries (in matrix P_xyz) and
% luminances at max power in (vector RGB_maxLum), bkg Lum

P_xyz = [rx gx bx; ry gy by; 1-rx-ry, 1-gx-gy, 1-bx-by];
P_LMS = T_XYZ2LMS * P_xyz;
P_LMS = P_LMS ./ sum(P_LMS(1:2,:), 1);

Wr = RGB_maxLum(1); Wg = RGB_maxLum(2); Wb = RGB_maxLum(3);

% main calculation

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

