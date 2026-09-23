function H = makeHartleyFromParams(hartleys60_meta)

% recreate Hartleys

% Felix's method

sf = hartleys60_meta(1);
ori = hartleys60_meta(2);
shift =  hartleys60_meta(3);
chrom = hartleys60_meta(4);

use_nPix = 60;
L = 120;
xs = (0:L-1);
M = floor(use_nPix)/2;
kx = sf * use_nPix/L;
xs2 = shift + xs*2*pi*(kx)/M;

if sf == 30 && shift == 0
    xs2 = pi/2 + xs*2*pi*(kx)/M;
elseif sf == 30 && shift == pi
    xs2 = 3*pi/2 + xs*2*pi*(kx)/M;
end

Hs = sin(xs2);
Hs2dcis2 = imrotate(repmat(Hs, L,1), -ori, 'nearest', 'crop');
corevec = floor((1:use_nPix)+(L-use_nPix)/2);
hartley60 = Hs2dcis2(corevec, corevec);

temp = zeros([size(hartley60) 3]);
temp(:,:,chrom) = hartley60;
H = temp;


% What would be sensible:
% 
% [X,Y] = meshgrid(1:60);
% Xp = X.*cos(deg2rad(ori)) + Y.*sin(deg2rad(ori));
% hartley60_sensible = sin(2*pi*sf/60 .* Xp - deg2rad(ori));
% 
end