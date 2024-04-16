%% pregen Achrom stimuli
cd('Z:\StimulusSet\NTlab_cis\Cloudstims_calib_01_2022')

cloudpix = 60;
for cur_Scale = [1:6]
Dualstim_pregen_achromcloud_n = 4000;
for init=1:10
    DensenoiseAchromcloud = (mk_spatialcloud(cloudpix,cloudpix, Dualstim_pregen_achromcloud_n, cur_Scale)./2 +.5).*255;
    [~,DensenoiseAchromcloud_binned] = histc(DensenoiseAchromcloud,linspace(0,255,256));
    
    save(sprintf('Cloudstims_Achrom_size%d_scale%d_%02d.mat', cloudpix, cur_Scale, init), 'DensenoiseAchromcloud_binned');
    %save(['Cloudstims_Achrom_size' num2str(cloudpix) '_scale' num2str(cur_Scale) '_' num2str( '.mat'])
    clearvars DensenoiseAchromcloud DensenoiseAchromcloud_binned
    
end
end

%% pregen Chromcloud
cd('Z:\StimulusSet\NTlab_cis\Cloudstims_calib_01_2022')

global g_strctParadigm
g_strctParadigm.m_strctConversionMatrices.ldgyb = [    1.0000    1.0000    0.2099
    1.0000   -0.2331   -0.1858
    1.0000   0.0079    1.0000]; %felix added calibration as of jan 28 2022

cloudpix = 60;
for cur_Scale = [1:6]
Dualstim_pregen_chromcloud_n = 4000;
for init=1:10
    %DensenoiseAchromcloud = (mk_spatialcloud(cloudpix,cloudpix, Dualstim_pregen_achromcloud_n, DensenoiseScale)./2 +.5).*255;
    DensenoiseChromcloud_DKlspace=reshape(mk_spatialcloudRGB(cloudpix, cloudpix, Dualstim_pregen_chromcloud_n, cur_Scale),cloudpix*cloudpix*Dualstim_pregen_chromcloud_n,3);
    DensenoiseChromcloud_sums=sum(abs(DensenoiseChromcloud_DKlspace),2); DensenoiseChromcloud_sums(DensenoiseChromcloud_sums < 1)=1;
    DensenoiseChromcloud_DKlspace=DensenoiseChromcloud_DKlspace./[DensenoiseChromcloud_sums,DensenoiseChromcloud_sums,DensenoiseChromcloud_sums];
    DensenoiseChromcloud=reshape(round(255.*ldrgyv2rgb(DensenoiseChromcloud_DKlspace(:,1)',DensenoiseChromcloud_DKlspace(:,2)',DensenoiseChromcloud_DKlspace(:,3)'))',cloudpix,cloudpix,Dualstim_pregen_chromcloud_n,3);
    DensenoiseChromcloud_DKlspace=reshape(DensenoiseChromcloud_DKlspace,cloudpix,cloudpix,Dualstim_pregen_chromcloud_n,3);
    
    save(sprintf('Cloudstims_Chrom_size%d_scale%d_%02d.mat', cloudpix, cur_Scale, init), 'DensenoiseChromcloud','DensenoiseChromcloud_DKlspace');
    clearvars DensenoiseChromcloud DensenoiseChromcloud_DKlspace DensenoiseChromcloud_sums
    
end
disp('Done wtih scale %02d',cur_Scale)
end

%% pregen binary Chromcloud
% CCs = generate_color_cloud_binary( NT, NX, CLRscale, SPscale, THRESH )
% to test: 
% test1=generate_color_cloud_binary(20, 60, 4, 6);
% figure; imagesc([squeeze(test1(:,:,1,:)./2)]+.5)
cd('Z:\StimulusSet\NTlab_cis\Cloudstims_calib_01_2022')

global g_strctParadigm
g_strctParadigm.m_strctConversionMatrices.ldgyb = [    1.0000    1.0000    0.2099
    1.0000   -0.2331   -0.1858
    1.0000   0.0079    1.0000]; %felix added calibration as of jan 28 2022

cloudpix = 60;

cur_SPscale = 6;
for cur_Scale = [1:6]
Dualstim_pregen_chromcloud_n = 4000;
for init=1:10
    %DensenoiseAchromcloud = (mk_spatialcloud(cloudpix,cloudpix, Dualstim_pregen_achromcloud_n, DensenoiseScale)./2 +.5).*255;
    DensenoiseChromcloud_DKlspace=reshape(generate_color_cloud_binary(Dualstim_pregen_chromcloud_n, cloudpix, cur_Scale, cur_SPscale),cloudpix*cloudpix*Dualstim_pregen_chromcloud_n,3);
    DensenoiseChromcloud=reshape(round(255.*ldrgyv2rgb(DensenoiseChromcloud_DKlspace(:,1)',DensenoiseChromcloud_DKlspace(:,2)',DensenoiseChromcloud_DKlspace(:,3)'))',cloudpix,cloudpix,Dualstim_pregen_chromcloud_n,3);
    DensenoiseChromcloud_DKlspace=reshape(DensenoiseChromcloud_DKlspace,cloudpix,cloudpix,Dualstim_pregen_chromcloud_n,3);
    
    save(sprintf('Cloudstims_BinaryChrom_size%d_scale%d_SPscale%d_%02d.mat', cloudpix, cur_Scale, cur_SPscale, init), 'DensenoiseChromcloud','DensenoiseChromcloud_DKlspace');
    clearvars DensenoiseChromcloud DensenoiseChromcloud_DKlspace
    
end
disp('Done wtih scale %02d',cur_Scale)
end

%%
