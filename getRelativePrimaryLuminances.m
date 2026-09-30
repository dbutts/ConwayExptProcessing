function lumMat = getRelativePrimaryLuminances(xyMat)

% given xyz chromaticity coordinates of primaries and white point, find the
% relative luminances of the primaries assuming they contribute equal power
% to the white point

rx = xyMat(1,1); ry = xyMat(1,2); 
gx = xyMat(2,1); gy = xyMat(2,2);
bx = xyMat(3,1); by = xyMat(3,2);
wx = xyMat(4,1); wy = xyMat(4,2);

T_XYZ2LMS = [0.15514 0.54312 -0.03286;
            -0.15514 0.45684 0.03286;
             0 0 0.01608];

 PW_xyz = [rx gx bx wx;
          ry gy by wy;
          1-rx-ry, 1-gx-gy, 1-bx-by, 1 - wx - wy];

PW_MB_num = T_XYZ2LMS * PW_xyz;
PW_MB_denom = sum(PW_MB_num(1:2,:), 1);

% columns of P_MB give R, G, B, W MacLeod-Boynton coordinates
PW_MB = PW_MB_num ./ PW_MB_denom;

rR = PW_MB(1,1); rG = PW_MB(1,2); rB = PW_MB(1,3); rW = PW_MB(1,4);
bR = PW_MB(3,1); bG = PW_MB(3,2); bB = PW_MB(3,3); bW = PW_MB(3,4);

slope_BW = (bW - bB)/(rW - rB);
intercept_BW = (rW*bB - rB*bW)/(rW - rB);

slope_GR = (bR - bG)/(rR - rG);
intercept_GR = (rR*bG - rG*bR)/(rR - rG);

% point  a long GR axis that specifies relative amounts of R and G in white
% point
W_GR = [-slope_BW 1; -slope_GR 1] \ [intercept_BW; intercept_GR];

% ratio of green lum to red lum -- distance from R to W_GR divided by
% distance from G to W_GR

GR_lumRatio = sqrt(sum((W_GR - [rR; bR]).^2)) / sqrt(sum((W_GR - [rG; bG]).^2));

% ratio of R+G lum to blue lum

GRB_lumRatio = sqrt(sum(([rB; bB] - [rW; bW]).^2)) / sqrt(sum((W_GR-[rW; bW]).^2));

% let lumR be 1

lumR = 1;
lumG = GR_lumRatio;
lumB = (GR_lumRatio + lumR)./(GRB_lumRatio);

Rrel = lumR/(lumR+lumG+lumB);
Grel = lumG/(lumR+lumG+lumB);
Brel = lumB/(lumR+lumG+lumB);

lumMat = [Rrel; Grel; Brel];