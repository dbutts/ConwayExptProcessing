%% loading and saving Geiusler and Burge natural scene statistics dataset
% see https://natural-scenes.cps.utexas.edu/db.shtml
% cite:  Geisler WS & Perry JS (2011). Statistics for optimal point prediction in natural images. Journal of Vision. October 19, 2011 vol. 11 no. 12 article 14.

img_dir = '/home/conwaylab/Documents/Natural images/cps20100428.ppm';
cd(img_dir)

all_images = dir('*.ppm');
num_images=length(all_images);

%% load single m image to test
cur_image = imread(all_images(1).name);

%% load all images and check plots
for ii=1:num_images
    [cur_image, cur_map] = imread(all_images(ii).name, 'ppm');
    figure(1)
    subplot(1,4,1); imagesc(squeeze(cur_image(:,:,1)))
    subplot(1,4,2); imagesc(squeeze(cur_image(:,:,2)))
    subplot(1,4,3); imagesc(squeeze(cur_image(:,:,3)))
    subplot(1,4,4); image(cur_image)
pause
end


%% Loading the UPenn natural image datasets
% from https://web.sas.upenn.edu/upennidb/albums/
% cite: Tkacik G et al, “Natural images from the birthplace of the human eye“, PLoS ONE 6: e20409 (2011).
clear all
%%
cur_set = 'cd05A';
img_dir = ['/home/conwaylab/Documents/Natural images/' cur_set];
cd(img_dir)

all_images = dir('*_RGB.mat');
all_images_jpg = dir('*.JPG');
num_images=length(all_images);

global g_strctParadigm
g_strctParadigm.m_strctConversionMatrices.ldgyb = [    1.0000    1.0000    0.1982
    1.0000   -0.2472   -0.1876
    1.0000   0.0097    1.0000]; %felix added calibration as of jan 10 2024
%% load single image to test

cur_image = load(all_images(1).name);
cur_image_jpg = imread(all_images_jpg(1).name);
figure;

subplot(2,1,1); imagesc(cur_image.RGB_Image./16384)
subplot(2,1,2); image(cur_image_jpg)

% %% load all images from RGb file
% 
% for ii=1:num_images
%     cur_image = load(all_images(ii).name);
%     use_image = cur_image.RGB_Image./2^14;
% end

%% load all images from JPG file

for ii=1:num_images
    cur_image = imread(all_images_jpg(ii).name);
        img_dims = size(cur_image);
    use_image = double(reshape(cur_image./255, [img_dims(1)*img_dims(2), 3]));
    [DKL_image(1,:), DKL_image(2,:), DKL_image(3,:)] = rgb2ldrgyv(use_image');
    images_RGB(ii,:,:,:) = cur_image;
    images_DKL(ii,:,:,:) = reshape(DKL_image',[img_dims(1), img_dims(2),3]);
end

% test image:
%figure; imagesc(squeeze(images_DKL(1,:,:,2)))

%%
save([cur_set '_rgb.mat'], 'images_RGB', '-v7.3');
save([cur_set '_all.mat'], 'images_RGB', 'images_DKL', '-v7.3');
%%
disp('Done saving images!')

