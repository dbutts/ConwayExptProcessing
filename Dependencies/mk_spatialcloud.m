function [ stim ] = mk_spatialcloud( stimw, stimh, num_frames, spatial_scale, seedn )
%mk_spatialcloudRGB creates spatial cloud stimulus with three color
%channels
if nargin==5; rng(seedn); end

stim = zeros(stimw,stimh,num_frames);
%spatial_scale = 5;

% Start with gaussian white noise (on larger square)
L = max([stimw stimh]);
noise_stim = randn(L,L,num_frames);

% 2-D Gaussian mask
xs=(0:(L-1))-L/2;
r2s = xs'.^2*ones(1,L) + ones(L,1)*xs.^2;
rad1 = 2*L/pi/spatial_scale;
mask1 = exp(-r2s/(2*rad1^2));
mask1 = mask1/max(mask1(:));

for k=	1:num_frames
	im1 = squeeze(noise_stim(:,:,k));
	manip1 = fftshift(fft2(im1)); 
	manip2 = mask1.*manip1;
%	im2 = abs(ifft2(ifftshift(manip2)));
	im2 = ifft2(ifftshift(manip2),'symmetric');
	stim(:,:,k) = im2(1:stimw,1:stimh)./max(abs(im2(:)));
end

%stim = (stim/std(stim(:)))/max(abs(stim(:)));
%stim = stim./max(abs(stim(:)));
end

