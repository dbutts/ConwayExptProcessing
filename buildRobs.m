function Robs_strct = buildRobs(spkDataXline, stimTiming, isTrialOfInterest)


% create by spkData = organizeSpikeDataByTrial(stimIntervals,plexon_fname, ks_path)

% input must be either spkData.spkDataOnline or spkData.spkDataffline

if nargin < 3
    isTrialOfInterest = true(numel(stimIntervals)/2,1);
end

stimStartTimes = stimTiming.stimStartTimes;
stimStopTimes = stimTiming.stimStopTimes;
numFrames = stimTiming.numFrames;
num_arrays = numel(spkDataXline);

stimFrameBinEdges = cellfun(@(startTime, stopTime, nFrames) linspace(startTime, stopTime, nFrames+1), ...
    num2cell(stimStartTimes), ...
    num2cell(stimStopTimes),...
    num2cell(numFrames), 'UniformOutput', false);

% for each unit, bin spike times by frame bins to get spikes it produced on
% each frame


idx = find(isTrialOfInterest);
for a = 1:num_arrays

    nSU = numel(spkDataXline(a).SU_clusters);
    nMU = numel(spkDataXline(a).MU_clusters);

    all_clusters = [spkDataXline(a).SU_clusters; spkDataXline(a).MU_clusters];

    for unit= 1:numel(all_clusters)
        unitID = all_clusters(unit);
        foo = cellfun(@(spk,clust,frameBins) histcounts(spk(clust==unitID), frameBins), spkDataXline(a).spk_times_cellArray(2*idx), spkDataXline(a).spk_clusters_cellArray(2*idx), stimFrameBinEdges(idx), 'UniformOutput',false);
        Robs(unit,:) = [foo{:}];
    end

    RobsSU{a} = Robs(1:nSU,:);
    RobsMU{a} = Robs(nSU+1 : nSU + nMU,:);

end

Robs_strct.Robs = [vertcat(RobsSU{:}); vertcat(RobsMU{:})];
Robs_strct.SU_clusters = [spkDataXline.SU_clusters];
Robs_strct.MU_clusters = [spkDataXline.MU_clusters];
Robs_strct.SU_chans = [spkDataXline.SU_chans];
Robs_strct.MU_chans = [spkDataXline.MU_chans];


end
