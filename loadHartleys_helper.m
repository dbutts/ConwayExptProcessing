function [stim_array, cache] = loadHartleys_helper(stimpath, stimseq, cache)

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
    stim = tmp.hartleys60_meta;
    cache(filename) = stim;
end

%stim_cellArray = cellfun(@(idx) stim(idx,:), stimseq, 'UniformOutput', false);
%stim_cellArray = cellfun(@transpose, stim_cellArray, 'UniformOutput', false);
stim_array  = transpose(stim(stimseq,:));

%stim_cellArray = cellfun(@(idx) transpose(stim(idx,:)), stimseq, 'UniformOutput', false);
end
