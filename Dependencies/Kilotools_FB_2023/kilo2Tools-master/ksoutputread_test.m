%%
cd('/media/felix/Internal_1/Data/BevilColor/220209_163451_Jacomo/kilosorting_laminar')
%spk_templates = readNPY('spike_templates.npy');
spk_times = readNPY('spike_times_seconds.npy');
spk_clusters = readNPY('spike_clusters.npy');
spk_labels = tdfread('cluster_KSLabel.tsv');
spk_clustIDs = unique(spk_clusters); nclusts=length(spk_clustIDs);

spk_labels_SU=[]; spk_labels_MU=[];
for cc=1:nclusts
    if strcmp(spk_labels.KSLabel(cc,:), 'good')
        spk_labels_SU = [spk_labels_SU,cc];
    elseif strcmp(spk_labels.KSLabel(cc,:), 'mua ')
        spk_labels_MU = [spk_labels_MU,cc];
    end
end
%%
unique(spk_clusters)