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

   
    spk_clusters_ontrialsOfInterest = vertcat(spkDataXline(a).spk_clusters_cellArray{2*idx});
    SU_clusters_ontrialsOfInterest = spk_clusters_ontrialsOfInterest(ismember(spk_clusters_ontrialsOfInterest, spkDataXline(a).SU_clusters));
    MU_clusters_ontrialsOfInterest = spk_clusters_ontrialsOfInterest(ismember(spk_clusters_ontrialsOfInterest, spkDataXline(a).MU_clusters));

    SU_clusters{a} = unique(SU_clusters_ontrialsOfInterest);
    MU_clusters{a} = unique(MU_clusters_ontrialsOfInterest);

    SU_chans{a} = spkDataXline(a).SU_chans(ismember(SU_clusters{a}, spkDataXline(a).SU_clusters));
    MU_chans{a} = spkDataXline(a).MU_chans(ismember(MU_clusters{a}, spkDataXline(a).MU_clusters));

    

    nSU = numel(SU_clusters{a});
    nMU = numel(MU_clusters{a});
    SU_arrays{a} = repmat(spkDataXline(a).arrayLabel, nSU,1);
    MU_arrays{a} = repmat(spkDataXline(a).arrayLabel, nMU,1);

    all_clusters = [SU_clusters{a}; MU_clusters{a}];

    for unit= 1:numel(all_clusters)
        unitID = all_clusters(unit);
        foo = cellfun(@(spk,clust,frameBins) histcounts(spk(clust==unitID), frameBins), spkDataXline(a).spk_times_cellArray(2*idx), spkDataXline(a).spk_clusters_cellArray(2*idx), stimFrameBinEdges(idx), 'UniformOutput',false);
        Robs(unit,:) = [foo{:}];
    end

    RobsSU{a} = Robs(1:nSU,:);
    RobsMU{a} = Robs(nSU+1 : nSU + nMU,:);

end

Robs_strct.Robs = [vertcat(RobsSU{:}); vertcat(RobsMU{:})];
Robs_strct.SU_clusters = vertcat(SU_clusters{:});
Robs_strct.MU_clusters = vertcat(MU_clusters{:});
Robs_strct.SU_chans = vertcat(SU_chans{:});
Robs_strct.MU_chans = vertcat(MU_chans{:});
Robs_strct.SU_arrays = vertcat(SU_arrays{:});
Robs_strct.MU_arrays = vertcat(MU_arrays{:});

end
