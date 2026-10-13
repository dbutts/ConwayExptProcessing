function H = makeHartleyFromParams(hartleys_meta, sideLength, scalar, method)
% mjg 10/1/26

% Defaults: use Felix's method, assume 60 pixel hartley
if nargin < 4
    method = 1;
end

if nargin < 3
    scalar = 1;
end

if nargin < 2
    sideLength = 60;
end

% Initialize array of hartleys
num_chrom_channels = 3;
num_hartleys = size(hartleys_meta,2);
H = zeros(sideLength,sideLength,num_chrom_channels, num_hartleys, 'uint8');

% recreate Hartleys

for i = 1:num_hartleys
    % get hartley parameters
    sf = hartleys_meta(1,i);
    ori = hartleys_meta(2,i);
    shift =  hartleys_meta(3,i);
    chrom = hartleys_meta(4,i);

    if method == 1 % Felix's method

        use_nPix = sideLength;
        L = 2*sideLength;
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
        hartley = uint8(127 * Hs2dcis2(corevec, corevec));

    else % the sensible way
        [X,Y] = meshgrid((0:sideLength-1) - sideLength/2);
        Xp = X.*cos(deg2rad(ori)) + Y.*sin(deg2rad(ori));
        hartley = uint8(127*sin(2*pi*sf/60 .* Xp + shift));
    end

    temp = zeros([size(hartley) 3]);
    temp(:,:,chrom) = hartley;
    H(:,:,:,i) = temp;
end

% resize if necessary, I think Psychtoolbox uses bilinear interpolation
if scalar ~= 1
    H = imresize(H, scalar, 'bilinear');
else
    sz = size(H);
    newRows = ceil(scalar * sz(1));
    newCols = ceil(scalar * sz(2));
    [x,y] = meshgrid(((1:newCols) - 0.5) / scalar + 0.5, ...
        ((1:newRows) - 0.5) / scalar + 0.5);
    x = min(max(x, 1), sz(2));
    y = min(max(y, 1), sz(1));
    
    slices = reshape(H, sz(1), sz(2), []);
    out = zeros(newRows, newCols, size(slices,3), 'like', H);
    
    for k = 1:size(slices,3)
        out(:,:,k) = interp2(double(slices(:,:,k)), x, y, 'linear');
    end
    
    H = reshape(out, [newRows, newCols, sz(3:end)]);
end
end