function Robs = buildRobs(ks_path, stimIntervals, fs, isTrialOfInterest)

if nargin < 4
    isTrialOfInterest = true(numel(stimIntervals)/2,1);
end

spkData = organizeSpikeDataByTrial(ks_path, stimIntervals, fs);
num_arrays = numel(spkData);


stimFrameBinEdges = cellfun(@(startTime, stopTime, nFrames) linspace(startTime, stopTime, nFrames+1), ...
    num2cell(stimStartTimes), ...
    num2cell(stimStopTimes),...
    num2cell(numFrames'), 'UniformOutput', false);

% for each unit, bin spike times by frame bins to get spikes it produced on
% each frame

Robs = cell(numel(spkData),1);

idx = 2*find(isTrialOfInterest);
for a = 1:num_arrays

    nSU = numel(spkData(a).SU_clusters);
    nMU = numel(spkData(a).MU_clusters);

    all_clusters = [spkData(a).SU_clusters; spkData(a).MU_clusters];

    for unit= 1:numel(all_clusters)
        unitID = all_clusters(unit);
        foo = cellfun(@(spk,clust,frameBins) histcounts(spk(clust==unitID), frameBins), spk_times_cellArray(idx), spk_clusters_cellArray(idx), stimFrameBinEdges, 'UniformOutput',false);
        Robs(unit,:) = [foo{:}];
    end

    RobsSU{a} = Robs(1:nSU,:);
    RobsMU{a} = Robs(nSU+1 : nSU + nMU,:);

end

Robs_strct.Robs = [vertcat(RobsSU{:}), vertcat(RobsMU{:})];
Robs_strct.SU_clusters = [spkData.SU_clusters];
Robs_strct.MU_clusters = [spkData.MU_clusters];
Robs_strct.SU_chans = [spkData.SU_chans];
Robs_strct.MU_chans = [spkData.MU_chans];

end
