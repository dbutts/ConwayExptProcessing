%% generate clouds
global g_strctParadigm
g_strctParadigm.m_strctConversionMatrices.ldgyb = [    1.0000    1.0000    0.2099
    1.0000   -0.2331   -0.1858
    1.0000   0.0079    1.0000]; %felix added calibration as of jan 28 2022

cloudpix = 100;
cloudscale = 6; 
cloud_n = 10;

init_cloud = mk_spatialcloudRGB(cloudpix, cloudpix, cloud_n, cloudscale);
%%
figure;
subplot(1,4,1); imagesc(squeeze(init_cloud(:,:,1,1))); title('Lum'); pbaspect([1,1,1]); colormap(gray)
subplot(1,4,2); imagesc(squeeze(init_cloud(:,:,1,2))); title('L-M'); pbaspect([1,1,1]); colormap(gray)
subplot(1,4,3); imagesc(squeeze(init_cloud(:,:,1,3))); title('S'); pbaspect([1,1,1]); colormap(gray)
subplot(1,4,4); image(squeeze(init_cloud(:,:,1,:))); pbaspect([1,1,1]); title('Full');
figtitle(['RGB cloud - Spatial Scale' num2str(cloudscale)])
%%
DensenoiseChromcloud_DKlspace=reshape(init_cloud,cloudpix*cloudpix*cloud_n,3);
%DensenoiseChromcloud_sums=sum(abs(DensenoiseChromcloud_DKlspace),2); DensenoiseChromcloud_sums(DensenoiseChromcloud_sums < 1)=1;
%DensenoiseChromcloud_DKlspace=DensenoiseChromcloud_DKlspace./[DensenoiseChromcloud_sums,DensenoiseChromcloud_sums,DensenoiseChromcloud_sums];
DensenoiseChromcloud=reshape(round(255.*ldrgyv2rgb(DensenoiseChromcloud_DKlspace(:,1)',DensenoiseChromcloud_DKlspace(:,2)',DensenoiseChromcloud_DKlspace(:,3)'))',cloudpix,cloudpix,cloud_n,3);
DensenoiseChromcloud_DKlspace=reshape(DensenoiseChromcloud_DKlspace,cloudpix,cloudpix,cloud_n,3);

%%
figure;
subplot(1,4,1); imagesc(squeeze(DensenoiseChromcloud(:,:,1,1))); title('Lum'); pbaspect([1,1,1]); colormap(gray)
subplot(1,4,2); imagesc(squeeze(DensenoiseChromcloud(:,:,1,2))); title('L-M'); pbaspect([1,1,1]); colormap(gray)
subplot(1,4,3); imagesc(squeeze(DensenoiseChromcloud(:,:,1,3))); title('S'); pbaspect([1,1,1]); colormap(gray)
subplot(1,4,4); image(squeeze(DensenoiseChromcloud(:,:,1,:)./255)); pbaspect([1,1,1]); title('Full');
figtitle(['DKL cloud - Spatial Scale' num2str(cloudscale)])

%%
DensenoiseChromcloud_DKlspace=reshape(init_cloud,cloudpix*cloudpix*cloud_n,3);
DensenoiseChromcloud_sums=sum(abs(DensenoiseChromcloud_DKlspace),2); DensenoiseChromcloud_sums(DensenoiseChromcloud_sums < 1)=1;
DensenoiseChromcloud_DKlspace=DensenoiseChromcloud_DKlspace./[DensenoiseChromcloud_sums,DensenoiseChromcloud_sums,DensenoiseChromcloud_sums];
DensenoiseChromcloud=reshape(round(255.*ldrgyv2rgb(DensenoiseChromcloud_DKlspace(:,1)',DensenoiseChromcloud_DKlspace(:,2)',DensenoiseChromcloud_DKlspace(:,3)'))',cloudpix,cloudpix,cloud_n,3);
DensenoiseChromcloud_DKlspace=reshape(DensenoiseChromcloud_DKlspace,cloudpix,cloudpix,cloud_n,3);
%%
figure;
subplot(1,4,1); imagesc(squeeze(DensenoiseChromcloud(:,:,1,1))); title('Lum'); pbaspect([1,1,1]); colormap(gray)
subplot(1,4,2); imagesc(squeeze(DensenoiseChromcloud(:,:,1,2))); title('L-M'); pbaspect([1,1,1]); colormap(gray)
subplot(1,4,3); imagesc(squeeze(DensenoiseChromcloud(:,:,1,3))); title('S'); pbaspect([1,1,1]); colormap(gray)
subplot(1,4,4); image(squeeze(DensenoiseChromcloud(:,:,1,:)./255)); pbaspect([1,1,1]); title('Full');
figtitle(['DKL cloud with norm- Spatial Scale' num2str(cloudscale)])
