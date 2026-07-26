function [stim_cellArray, cache] = loadHartleys_helper(stimpath, stimseq, imScalar, cache)

persistent hartleysFolder 
if isempty(hartleysFolder )
    hartleysFolder = fullfile(fileparts(stimpath), 'Cloudstims_calib_01_2022');
end

filename = fullfile(hartleysFolder, 'hartleys_60.mat');

if isKey(cache, filename)
    stim = cache(filename);
else
    tmp = load(filename);
    stim = tmp.hartleys60_DKL;
    cache(filename) = stim;
end

stim = permute(stim, [2 3 1 4]);
stim = imresize(stim, imScalar, 'bilinear');
stim =int8(127*stim);

stim_cellArray = cellfun(@(idx) stim(:,:,idx,:), stimseq, 'UniformOutput', false);
end
