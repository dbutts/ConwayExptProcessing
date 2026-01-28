function [spk_info, spk_times, spk_clusters] = fn_kiloappend(ksFilePath,ch_offset, spk_info, spk_times, spk_clusters)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

if nargin<=2

    spk_info = tdfread([ksFilePath 'cluster_info.tsv']);
        spk_info.ch = spk_info.ch+ch_offset;

    spk_times = readNPY([ksFilePath 'spike_times_seconds.npy']);
    spk_clusters = readNPY([ksFilePath 'spike_clusters.npy']);

else
    
    spk_times = [spk_times; readNPY([ksFilePath 'spike_times_seconds.npy'])];
    
    nclusts=max(spk_info.cluster_id)+1;
    spk_clusters = [spk_clusters; readNPY([ksFilePath 'spike_clusters.npy'])+nclusts];
    
    spk_info_temp = tdfread([ksFilePath 'cluster_info.tsv']);
    spk_info.cluster_id = [spk_info.cluster_id; spk_info_temp.cluster_id+nclusts];
    spk_info.Amplitude = [spk_info.Amplitude; spk_info_temp.Amplitude];
    spk_info.ContamPct = [spk_info.ContamPct; spk_info_temp.ContamPct];
    spk_info.KSLabel = [spk_info.KSLabel; spk_info_temp.KSLabel];
    %spk_info.Rating = [spk_info.Rating; spk_info_temp.r];
    spk_info.amp = [spk_info.amp; spk_info_temp.amp];
    spk_info.ch = [spk_info.ch; spk_info_temp.ch+ch_offset];
    spk_info.depth = [spk_info.depth; spk_info_temp.depth];
    spk_info.fr = [spk_info.fr; spk_info_temp.fr];
    try    
        spk_info.group = [spk_info.group; spk_info_temp.group];
    catch
        spk_info_temp.group=spk_info_temp.group(:,1:5)
        spk_info.group = [spk_info.group; spk_info_temp.group];
    end
    spk_info.n_spikes = [spk_info.n_spikes; spk_info_temp.n_spikes];
    spk_info.sh = [spk_info.sh; spk_info_temp.sh];

end

return