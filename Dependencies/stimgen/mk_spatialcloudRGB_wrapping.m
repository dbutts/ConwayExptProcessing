function [ stim ] = mk_spatialcloudRGB_wrapping( stimw, stimh, num_frames, spatial_scale, seedn )
%mk_spatialcloudRGB creates spatial cloud stimulus with three color
%channels
if nargin==5; rng(seedn); end

stim = zeros(stimw,stimh,num_frames,3);

% Start with gaussian white noise (on larger square)
L = max([stimw stimh]);
noise_stim = randn(L,L,num_frames,3);

% Make GWN continuously tilable; added by Ethan 05/23/2023
big = cat(2,cat(1,noise_stim,noise_stim),cat(1,noise_stim,noise_stim));
noise_stim = big(stimh/2+1:stimh*1.5,stimw/2+1:stimw*1.5,:,:);

% 2-D Gaussian mask
xs = (0:(L-1))-L/2;
r2s = xs'.^2*ones(1,L) + ones(L,1)*xs.^2;
rad1 = 2*L/pi/spatial_scale;
mask1 = exp(-r2s/(2*rad1^2));
mask1 = mask1/max(mask1(:));

% Fourier transform
for cc = 1:3
    for k =	1:num_frames
	    im1 = squeeze(squeeze(noise_stim(:,:,k,cc)));
	    manip1 = fftshift(fft2(im1)); 
	    manip2 = mask1.*manip1;
	    im2 = (ifft2(ifftshift(manip2),'symmetric'));
	    stim(:,:,k,cc) = im2(1:stimw,1:stimh)./max(abs(im2(:)));
    end
end

%stim = (stim/std(stim(:)))/(max(stim(:)));
%stim = stim./max(stim(:));
end
