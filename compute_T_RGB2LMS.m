function T_RGB2LMS = compute_T_RGB2LMS(varargin)

% Option set 1: wavelengths, red_spd, green_spd, blue_spd
% Option set 2: xyMat = [rx ry; gx gy; bx by], lumMat = [R_maxLumCdm2;
% G_maxLumCdm2; B_maxLumCdm2]
% B_maxLumCdm2
% Output: Matrix containing cone sensitivities to display primaries
% for conversion from RGB to LMS. Rows correspond to L, M, S and
% columns to R, G, B

% written by mjg 6/27/26

% Note: spds must be raw spds in units of W/Sr/m2
% function to be run in folder containing required PTB files:
% T_xyzJuddVos.mat, SplineCmf.m, SToWls.m

%% Parse inputs
p = inputParser;

addParameter(p, 'wavelengths', []);
addParameter(p, 'red_spd', []);
addParameter(p, 'green_spd', []);
addParameter(p, 'blue_spd', []);
addParameter(p, 'xyMat', []);
addParameter(p, 'lumMat', []);

parse(p, varargin{:});

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

T_XYZ2LMS = [0.15514 0.54312 -0.03286;
    -0.15514 0.45684 0.03286;
    0 0 0.01608];

if ~isempty(wavelengths) && ~isempty(red_spd) && ~isempty(green_spd) && ~isempty(blue_spd)
    T_XYZ = SplineCmf(SToWls(S_xyzJuddVos), T_xyzJuddVos, wavelengths);
    RGB = [red_spd green_spd blue_spd];
    P_XYZ = T_XYZ * RGB;
    P_xyz = P_XYZ./sum(P_XYZ,1);

    load T_xyz1931
    T_xyz1931_splined = SplineCmf(SToWls(S_xyz1931), T_xyz1931, wavelengths);

    R_maxLumCdm2=683 * dot(red_spd, T_xyz1931_splined(2,:));
    G_maxLumCdm2=683 * dot(green_spd, T_xyz1931_splined(2,:));
    B_maxLumCdm2=683 * dot(blue_spd, T_xyz1931_splined(2,:));
    %T_RGB2LMS = T_XYZ2LMS * P_XYZ;
else
    try
        P_xyz = [rx gx bx;
            ry gy by;
            1-rx-ry, 1-gx-gy, 1-bx-by];


    catch
        error('Missing/incorrect input arguments')
    end
end

P_MB_num = T_XYZ2LMS * P_xyz;
P_MB_denom = sum(P_MB_num(1:2,:), 1);
P_MB = P_MB_num ./ P_MB_denom;
T_RGB2LMS = P_MB .* [R_maxLumCdm2 G_maxLumCdm2 B_maxLumCdm2];
% Normalize
%T_RGB2LMS = T_RGB2LMS ./ max(T_RGB2LMS(:));

end


