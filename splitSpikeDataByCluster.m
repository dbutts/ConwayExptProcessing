function spike_ts_by_cluster_by_trial_cellArray = splitSpikeDataByCluster(spkDataXline, stimTiming, rawOrRelative)

% spikes on each trial for cluster 1
%cellfun(@(spike_ts, cluster) spike_ts(cluster == 1), spkData.spkDataOnline.spk_times_cellArray,  spkData.spkDataOnline.spk_clusters_cellArray, 'UniformOutput', false);

clusters = [spkDataXline.SU_clusters; spkDataXline.MU_clusters];
stim_start_times_rep = repelem(stimTiming.stimStartTimes, 2);
stim_start_times_cell = num2cell([0; stim_start_times_rep(1:end-1)]);
spike_ts_by_cluster_by_trial_cellArray = cell(size(spkDataXline.spk_clusters_cellArray,1), numel(clusters));

for i = 1:numel(clusters)

    if strcmpi(rawOrRelative, 'raw')

        spike_ts_by_cluster_by_trial_cellArray(:,i) =...
            cellfun(@(spike_ts, cluster) spike_ts(cluster == i),...
            spkDataXline.spk_times_cellArray, ...
            spkDataXline.spk_clusters_cellArray, ...
            'UniformOutput', false);

    elseif strcmpi(rawOrRelative, 'relative')

        spike_ts_by_cluster_by_trial_cellArray(:,i) =...
            cellfun(@(spike_ts, cluster, stim_start) spike_ts(cluster == i) - stim_start,...
            spkDataXline.spk_times_cellArray, ...
            spkDataXline.spk_clusters_cellArray, ...
            stim_start_times_cell,...
            'UniformOutput', false);

    end

end