% spikes on each trial for cluster 1
%cellfun(@(spike_ts, cluster) spike_ts(cluster == 1), spkData.spkDataOnline.spk_times_cellArray,  spkData.spkDataOnline.spk_clusters_cellArray, 'UniformOutput', false);

clusters = [spkData.spkDataOnline.SU_clusters; spkData.spkDataOnline.MU_clusters];
raw_spike_ts_by_cluster_by_trial_cellArray = cell(2*size(trial,2), numel(clusters));
relative_spike_ts_by_cluster_by_trial_cellArray = cell(2*size(trial,2), numel(clusters));
stim_start_times_rep = repelem(stimTiming.stimStartTimes, 2);
stim_start_times_cell = num2cell([0; stim_start_times_rep(1:end-1)]);

for i = 1:numel(clusters)

    raw_spike_ts_by_cluster_by_trial_cellArray(:,i) =...
        cellfun(@(spike_ts, cluster) spike_ts(cluster == i),...
        spkData.spkDataOnline.spk_times_cellArray, ...
        spkData.spkDataOnline.spk_clusters_cellArray, ...
        'UniformOutput', false);
   
    relative_spike_ts_by_cluster_by_trial_cellArray(:,i) =...
        cellfun(@(spike_ts, cluster, stim_start) spike_ts(cluster == i) - stim_start,...
        spkData.spkDataOnline.spk_times_cellArray, ...
        spkData.spkDataOnline.spk_clusters_cellArray, ...
       stim_start_times_cell,...
        'UniformOutput', false);
end

% % combine stim onn  and stim off spike times
relative_spike_ts_by_cluster_by_trial_cellArray = cellfun(@(a,b) [a; b], relative_spike_ts_by_cluster_by_trial_cellArray(1:2:end,:), relative_spike_ts_by_cluster_by_trial_cellArray(2:2:end,:), ...
'UniformOutput', false);

raw_spike_ts_by_cluster_by_trial_cellArray = cellfun(@(a,b) [a; b], ...
    raw_spike_ts_by_cluster_by_trial_cellArray(1:2:end,:),...
    raw_spike_ts_by_cluster_by_trial_cellArray(2:2:end,:), ...
'UniformOutput', false);
%%
figure; hold on
cluster = 25;
for i = 1:numel(find(hartleyTrialIdx))
    plot(relative_spike_ts_by_cluster_by_trial_cellArray{i,cluster},...
        i * ones(size(relative_spike_ts_by_cluster_by_trial_cellArray{i,cluster})), 'ko')

end
xlim([0 5])

temp = cellfun(@(raw_spk_ts, start, stop, nFrames) int8(histcounts(raw_spk_ts, 'binedges', linspace(start, stop, nFrames+1))),...
    raw_spike_ts_by_cluster_by_trial_cellArray,...
    num2cell(repmat(stimTiming.stimStartTimes, [1, numel(clusters)])),...
    num2cell(repmat(stimTiming.stimStopTimes, [1, numel(clusters)])),...
    num2cell(repmat(stimTiming.numFrames, [1, numel(clusters)])),...
    'UniformOutput', false);

% Disc probe
discProbeRGB = data_online(2).stim;

rx   = 0.6280;
ry   = 0.3310;
gx   = 0.3059;
gy   = 0.5826;
bx   = 0.1639;
by   = 0.0617;
R_maxLumCdm2 = 18.6469;
G_maxLumCdm2 = 75.8449;
B_maxLumCdm2 = 10.5313;
xyMat = [rx ry; gx gy; bx by];
lumMat = [R_maxLumCdm2; G_maxLumCdm2; B_maxLumCdm2];
T_DKL2RGB = compute_T_DKL2RGB('lumMat', lumMat, 'xyMat', xyMat);

discProbeDKL = T_DKL2RGB \ ((discProbeRGB - 128)./128);

discProbeHueAngle = mod(rad2deg(atan2(discProbeDKL(3,:), discProbeDKL(2,:))),360);
discProbeSat = sqrt(discProbeDKL(2,:).^2 + discProbeDKL(3,:).^2);

nHues = 16;

hues = linspace(0,360, nHues+1);
hues(hues==360) = [];
lumLevels = [ linspace(-0.8, -0.2, 3), 0,  linspace(0.2,0.8,3)]; % yup, these are the luminance levels for whatever reason
satLevels = [1/5, 1/3, 1/2, 2/3, 4/5, 1]; % Felix, what were you doing??
% figure out the intended hue

D_hue = pdist2(discProbeHueAngle', hues');
[~,I_hue]=min(D_hue, [], 2);
discProbeHueAngle_round = hues(I_hue);
discProbeHueAngle_round(discProbeHueAngle_round == 360) = 0;

D_lum =  pdist2(discProbeDKL(1,:)',lumLevels');
[~,I_lum]=min(D_lum, [], 2);
discProbeLum_round = lumLevels(I_lum);

D_sat =  pdist2(discProbeSat', satLevels');
[~,I_sat]=min(D_sat, [], 2);
discProbeSat_round = satLevels(I_sat);

% get firing rate for each cell and hue at sat = 1, lum = 0

temp_discProbe = cellfun(@(raw_spk_ts, start, stop) int8(histcounts(raw_spk_ts, 'binedges', start:0.02:stop)),...
    raw_spike_ts_by_cluster_by_trial_cellArray(discProbeTrialIdx,:),...
    num2cell(repmat(stimTiming.stimStartTimes(discProbeTrialIdx), [1, numel(clusters)])),...
    num2cell(repmat(stimTiming.stimStopTimes(discProbeTrialIdx), [1, numel(clusters)])),...
    'UniformOutput', false);


spikesPerTrial_discProbe = single(cellfun(@sum, temp_discProbe));


% for each hue get mean firing rate
for hue_idx = 1:numel(hues)
    temp = spikesPerTrial_discProbe(discProbeSat_round == 1 & discProbeLum_round == 0 & discProbeHueAngle_round == hues(hue_idx),:)./0.2;
    means(hue_idx,:) = mean(temp,1);
    stds(hue_idx,:) = std(temp,[],1);

end





Robs_discProbe = zeros(numel(clusters), sum(stimTiming.numFrames(ccloudTrialIdx)), 'int8');
for i = 1:size(temp,2)
    Robs(i,:) = horzcat(temp_discProbe{:,i});
end


