%% pregen Achrom stimuli
%cd('Z:\StimulusSet\NTlab_cis\Cloudstims_calib_04_2024')
cd('/home/conwaylab/Documents/Cloudstims_calib_04_2024')

global g_strctParadigm
g_strctParadigm.m_strctConversionMatrices.ldgyb = [    1.0000    1.0000    0.1982
    1.0000   -0.2472   -0.1876
    1.0000   0.0097    1.0000]; % Display++ calibration as of jan 10 2024

cloudpix = 60; % resolution of an individual cloud tile

for cur_Scale = [1:10] % loop across scales: note that each scale is multiplied by 2 in the cloud generation function to increase the range of spatial frequencies
    Dualstim_pregen_achromcloud_n = 4000;
    for init=1:10
        DensenoiseAchromcloud = (mk_spatialcloud(cloudpix,cloudpix, Dualstim_pregen_achromcloud_n, cur_Scale*2)./2 +.5).*255;
        [~,DensenoiseAchromcloud_binned] = histc(DensenoiseAchromcloud,linspace(0,255,256)); % bin the cloud to 256 levels for Bits++ compatibility
        
        save(sprintf('Cloudstims_Achrom_size%d_scale%d_%02d.mat', cloudpix, cur_Scale, init), 'DensenoiseAchromcloud_binned');
        clearvars DensenoiseAchromcloud DensenoiseAchromcloud_binned
        
    end
    fprintf('Done wtih lum clouds scale %02d \n',cur_Scale)
end

%% pregen Chromcloud

cloudpix = 60; % resolution of an individual cloud tile

for cur_Scale = [1:10]  % loop across scales: note that each scale is multiplied by 2 in the cloud generation function to increase the range of spatial frequencies

    Dualstim_pregen_chromcloud_n = 4000;
    for init=1:10
        DensenoiseChromcloud_DKlspace=reshape(mk_spatialcloudRGB_wrapping(cloudpix, cloudpix, Dualstim_pregen_chromcloud_n, cur_Scale*2),cloudpix*cloudpix*Dualstim_pregen_chromcloud_n,3);
        
        % correct DKL space projection for the case of all zeros to avoid displaying colors out of bounds
        DensenoiseChromcloud_sums=sum(abs(DensenoiseChromcloud_DKlspace),2); DensenoiseChromcloud_sums(DensenoiseChromcloud_sums < 1)=1; 
        DensenoiseChromcloud_DKlspace=DensenoiseChromcloud_DKlspace./[DensenoiseChromcloud_sums,DensenoiseChromcloud_sums,DensenoiseChromcloud_sums];

        % convert DKL space to RGB space
        DensenoiseChromcloud=reshape(round(255.*ldrgyv2rgb(DensenoiseChromcloud_DKlspace(:,1)',DensenoiseChromcloud_DKlspace(:,2)',DensenoiseChromcloud_DKlspace(:,3)'))',cloudpix,cloudpix,Dualstim_pregen_chromcloud_n,3);
        DensenoiseChromcloud_DKlspace=reshape(DensenoiseChromcloud_DKlspace,cloudpix,cloudpix,Dualstim_pregen_chromcloud_n,3);
        
        save(sprintf('Cloudstims_Chrom_size%d_scale%d_%02d.mat', cloudpix, cur_Scale, init), 'DensenoiseChromcloud','DensenoiseChromcloud_DKlspace');
        clearvars DensenoiseChromcloud DensenoiseChromcloud_DKlspace DensenoiseChromcloud_sums
        
    end
    fprintf('Done wtih color clouds scale %02d \n',cur_Scale)
end

%% pregen binary Chromcloud
% this is the same as above, but with a binary mask that sets certain regions to grey, thus inducing more edges in the cloud to drive cells that are suppressed by high SFs

cloudpix = 60; % resolution of an individual cloud tile

cur_SPscale = 6;
for cur_Scale = [1:10]
Dualstim_pregen_chromcloud_n = 4000;
for init=1:10
    %DensenoiseAchromcloud = (mk_spatialcloud(cloudpix,cloudpix, Dualstim_pregen_achromcloud_n, DensenoiseScale)./2 +.5).*255;
    DensenoiseChromcloud_DKlspace=reshape(mk_spatialcloudRGB_binary(Dualstim_pregen_chromcloud_n, cloudpix, cur_Scale*2, cur_SPscale),cloudpix*cloudpix*Dualstim_pregen_chromcloud_n,3);
    DensenoiseChromcloud=reshape(round(255.*ldrgyv2rgb(DensenoiseChromcloud_DKlspace(:,1)',DensenoiseChromcloud_DKlspace(:,2)',DensenoiseChromcloud_DKlspace(:,3)'))',cloudpix,cloudpix,Dualstim_pregen_chromcloud_n,3);
    DensenoiseChromcloud_DKlspace=reshape(DensenoiseChromcloud_DKlspace,cloudpix,cloudpix,Dualstim_pregen_chromcloud_n,3);
    
    save(sprintf('Cloudstims_BinaryChrom_size%d_scale%d_SPscale%d_%02d.mat', cloudpix, cur_Scale, cur_SPscale, init), 'DensenoiseChromcloud','DensenoiseChromcloud_DKlspace');
    clearvars DensenoiseChromcloud DensenoiseChromcloud_DKlspace
    
end
fprintf('Done wtih binary clouds scale %02d \n',cur_Scale)
end


%% pregen contrast-matched Chromcloud

LumScale = 0.1085; % maximum cone contrast for L-M excursion to scale down Luminance by
% this is based on the latest calibrations DeltaLMMax ./ DeltaLDmax (see calibmethod script for details)

cloudpix = 60;
for cur_Scale = [1:10]
Dualstim_pregen_chromcloud_n = 4000;
for init=1:10
    DensenoiseChromcloud_DKlspace=reshape(mk_spatialcloudRGB_wrapping(cloudpix, cloudpix, Dualstim_pregen_chromcloud_n, cur_Scale*2),cloudpix*cloudpix*Dualstim_pregen_chromcloud_n,3);
    DensenoiseChromcloud_DKlspace(:,1)=DensenoiseChromcloud_DKlspace(:,1)*LumScale;
    DensenoiseChromcloud_sums=sum(abs(DensenoiseChromcloud_DKlspace),2); DensenoiseChromcloud_sums(DensenoiseChromcloud_sums < 1)=1;
    DensenoiseChromcloud_DKlspace=DensenoiseChromcloud_DKlspace./[DensenoiseChromcloud_sums,DensenoiseChromcloud_sums,DensenoiseChromcloud_sums];
    DensenoiseChromcloud=reshape(round(255.*ldrgyv2rgb(DensenoiseChromcloud_DKlspace(:,1)',DensenoiseChromcloud_DKlspace(:,2)',DensenoiseChromcloud_DKlspace(:,3)'))',cloudpix,cloudpix,Dualstim_pregen_chromcloud_n,3);
    DensenoiseChromcloud_DKlspace=reshape(DensenoiseChromcloud_DKlspace,cloudpix,cloudpix,Dualstim_pregen_chromcloud_n,3);
    
    save(sprintf('Cloudstims_ContrastMatched_size%d_scale%d_%02d.mat', cloudpix, cur_Scale, init), 'DensenoiseChromcloud','DensenoiseChromcloud_DKlspace');
    clearvars DensenoiseChromcloud DensenoiseChromcloud_DKlspace DensenoiseChromcloud_sums
    
end
fprintf('Done wtih contrast-matched clouds scale %02d \n',cur_Scale)
end
%%
disp('all done!')
%%