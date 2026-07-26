function spkData = organizeSpikeDataByTrial(ks_path,stimIntervals,fs)

spike_times_dir = dir(fullfile(ks_path, '**/spike_times.npy'));
spike_times_folders = {spike_times_dir(:).folder};
assert(~isempty(spike_times_dir));
[~, ks_folders, ~] =  cellfun(@fileparts, spike_times_folders, 'UniformOutput', false);
tokens = regexp(ks_folders, '^[^_]+_([^_]+)', 'tokens');
array_labels = cellfun(@(t) t{1}{1}, tokens, 'UniformOutput', false);


% for each ks_folder, get spike data and organize it
% Get label of kilosort batch (often corresponding to array name and range of channels processed)

% Read in kilosort outputs

spkData = struct;

chan_offset = 0;
cluster_offset = 0;
for f = 1:numel(ks_folders)

    spike_times_folder = spike_times_folders{f};
    spk_times = readNPY(fullfile(spike_times_folder, 'spike_times.npy'));
    spk_times = double(spk_times)./fs + spk_offset; % convert to seconds
    spk_clusters = readNPY(fullfile(spike_times_folder, 'spike_clusters.npy'))...
        + cluster_offset;
    cluster_KSLabel = tdfread(fullfile(spike_times_folder, 'cluster_KSLabel.tsv'));

    chan_map = readNPY(fullfile(spike_times_folder, 'channel_map.npy'));
    chan_map = chan_map + chan_offset; % make channel numbers unique within array

    if isfile(fullfile(spike_times_folder, 'cluster_info.tsv'))
        cluster_info = tdfread(fullfile(spike_times_folder, 'cluster_info.tsv'));
        cluster_id = cluster_info.cluster_id + cluster_offset;
        group = cluster_info.group;
        n_spikes = cluster_info.n_spikes;
        chan_best = cluster_info.ch + double(chan_offset);
    else
        cluster_group = tdfread(fullfile(spike_times_folder, 'cluster_group.tsv'));
        cluster_group.cluster_id = cluster_group.cluster_id + cluster_offset;

        % account for blank units which may not be in cluster_group
        blank_cluster_id = setdiff(unique(spk_clusters), cluster_group.cluster_id);
        temp_cluster_id = [cluster_group.cluster_id; blank_cluster_id];
        temp_group = cluster_group.group;
        temp_group(end+1:end+length(blank_cluster_id),:) = ' ';
        [temp_cluster_id_sorted, I] = sort(temp_cluster_id);
        temp_group_sorted = temp_group(I,:);
        cluster_id = temp_cluster_id_sorted;
        group = temp_group_sorted;

        % find best channel of each cluster
        templates = readNPY(fullfile(spike_times_folder, 'templates.npy'));
        % n_spikes = accumarray(spk_clusters+1, spk_clusters, [], @numel);
        % n_spikes = n_spikes(n_spikes>0);

        n_spikes = accumarray(spk_clusters+1, 1, [], @sum);
        n_spikes = n_spikes(cluster_id + 1);

        [~,I]= max(sum(templates.^2,2),[],3);
        chan_best = chan_map(I); % best channel for each unique cluster
    end

    % Find indices of units labeled "good", "mua", or ""
    isGood = cellfun(@(x) strcmpi(deblank(x), 'good'), cellstr(group));
    isMua = cellfun(@(x) strcmpi(deblank(x), 'mua'), cellstr(group));
    isBlank = cellfun(@(x) isempty(deblank(x)), cellstr(group));
    hasMinSpikes = n_spikes > minSpikes;


    % Get rid of spike times and clusterIDs that correspond to bad units
    spk_times = spk_times(ismember(spk_clusters, allUnit_clusterIDs));
    spk_clusters = spk_clusters(ismember(spk_clusters, allUnit_clusterIDs));

    % bin spikes by trial
    [~,~,spk_times_bin] = histcounts(spk_times, stimIntervals);

    % each cell gives spike times for stim ON and stim OFF periods

    spk_times_cellArray =accumarray(...
        spk_times_bin(:)+1, ...
        spk_times(:), ...
        [nBins + 1, 1], ...
        @(x){x}, ...
        {[]});

    % corresponding cluster IDs
    spk_clusters_cellArray = accumarray(...
        spk_times_bin(:)+1, ...
        spk_clusters(:), ...
        [nBins + 1, 1], ...
        @(x){x}, ...
        {[]});

    spkData(f).spk_times_cellArray = spk_times_cellArray;
    spkData(f).spk_clusters_cellArray = spk_clusters_cellArray;
    spkData(f).arrayLabel = array_labels{f};
    spkData(f).SU_clusters = cluster_id(isGood & hasMinSpikes);
    spkData(f).MU_clusters = cluster_id((isMua | isBlank) & hasMinSpikes);
    spkData(f).SU_chans = chan_best(isGood & hasMinSpikes);
    spkData(f).MU_chans = chan_best((isMua | isBlank) & hasMinSpikes);

      % update chan_offset
    if f < numel(ks_folders) % if we still got a batch ahead
        next_array_label = array_labels{f+1};
        if strcmpi(next_array_label,  this_array_label)
            chan_offset = max(chan_map)+1; % maximum channel number;
        else
            chan_offset = 0; % new array, reset offset to 0
        end
        %   chan_offsets(ks_batch+1) = chan_offset;
    end
    % update cluster offset, and keep track of them
    cluster_offset = max(cluster_id)+1;
    cluster_offsets(f+1) = cluster_offset;
    chan_offsets(f+1) = chan_offset;


end

end