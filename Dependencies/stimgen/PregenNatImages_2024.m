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

img_dir = '/home/conwaylab/Documents/Natural images/cd01A';
cd(img_dir)

all_images = dir('*_RGB.mat');
num_images=length(all_images);

%% load single image to test

cur_image = load(all_images(1).name);
figure;
imagesc(cur_image.RGB_Image./16384)
%% load all images

for ii=1:num_images
    cur_image = load(all_images(ii).name);
    use_image = cur_image.RGB_Image./256;
end