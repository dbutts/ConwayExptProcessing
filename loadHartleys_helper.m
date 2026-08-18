function [stim_cellArray, cache] = loadHartleys_helper(stimpath, stimseq, cache)

persistent hartleysFolder 
if isempty(hartleysFolder )
    %hartleysFolder = fullfile(fileparts(stimpath), 'Cloudstims_calib_01_2022');
    hartleysFolder = stimpath;
end

filename = fullfile(hartleysFolder, 'hartleys_60.mat');

if isKey(cache, filename)
    stim = cache(filename);
else
    tmp = load(filename);
    stim = tmp.hartleys60_DKL;
    %stim = tmp.hartleys60_meta;
    cache(filename) = stim;
end

stim = permute(stim, [2 3 1 4]);
%stim = imresize(stim, imScalar, 'bilinear');
stim =int8(127*stim);

stim_cellArray = cellfun(@(idx) stim(:,:,idx,:), stimseq, 'UniformOutput', false);
stim_cellArray = cellfun(@(x) permute(x, [3 1 2 4]), stim_cellArray, 'UniformOutput', false);
stim_cellArray = cellfun(@(x) reshape(x, size(x,1), prod(size(x, 2:4))), stim_cellArray, 'UniformOutput', false);
stim_cellArray = cellfun(@transpose, stim_cellArray, 'UniformOutput', false);

%stim_cellArray = cellfun(@(idx) transpose(stim(idx,:)), stimseq, 'UniformOutput', false);
end
