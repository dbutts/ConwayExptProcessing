function cellArray = binSpikeTimes(spike_ts_by_cluster_by_trial_cellArray, startTimes, stopTimes, dt)

num_clusters =  size(spike_ts_by_cluster_by_trial_cellArray,2);


if isscalar(dt)
    dt = repmat(dt, [size(spike_ts_by_cluster_by_trial_cellArray,1), 1]);
end

if isscalar(startTimes)
     startTimes = repmat(startTimes, [size(spike_ts_by_cluster_by_trial_cellArray,1), 1]);
end


if isscalar(stopTimes)
    stopTimes = repmat(stopTimes, [size(spike_ts_by_cluster_by_trial_cellArray,1), 1]);
end


cellArray = cellfun(@(spk_ts, start, stop, dt) int16(histcounts(spk_ts, 'binedges', linspace(start, stop, (stop-start)/dt + 1))),...
    spike_ts_by_cluster_by_trial_cellArray,...
    num2cell(repmat(startTimes, [1, num_clusters])),...
    num2cell(repmat(stopTimes, [1, num_clusters])),...
     num2cell(repmat(dt, [1, num_clusters])),...
    'UniformOutput', false);

end