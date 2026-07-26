function [stim_cellArray, cache] = loadCclouds_helper(trial, stimpath, stimseq, imScalar, cache)

persistent ccloudStimFileStrs

if isempty(ccloudStimFileStrs)
    ccloudStimFileStrs = { ...
        'Cloudstims_Chrom_size60_scale%d_%02d.mat',...
        'Cloudstims_BinaryChrom_size60_scale%d_SPscale6_%02d.mat',...
        'Cloudstims_ContrastMatched_size60_scale%d_%02d.mat'};
end

filename = fullfile(stimpath,...
    sprintf(ccloudStimFileStrs{trial.usebinary+1},...
    trial.spatialscale, trial.BlockID));

if isKey(cache, filename)
    stim = cache(filename);
else
    tmp = load(filename);
    stim = tmp.DensenoiseChromcloud_DKlspace;
    cache(filename) = stim;
end

stim = imresize(stim, imScalar, 'bilinear');

switch trial.usebinary
    case {0,1} % full contrast and binary (?) clouds
        stim =int8(127*stim);

    case 2 % matched contrast
        LumScale = 0.1085;
        stim_int8 = zeros(size(stim), 'int8');
        stim_int8(:,:,:,2:3) = int8(127 * stim(:,:,:,2:3));
        stim_int8(:,:,:,1)   = int8((127 / LumScale) * stim(:,:,:,1));
        stim = stim_int8;
end

stim_cellArray = cellfun(@(idx) stim(:,:,idx,:), stimseq, 'UniformOutput', false);
end
