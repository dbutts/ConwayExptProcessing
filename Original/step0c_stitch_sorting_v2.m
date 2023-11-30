%%
clear all

%%
filenameP = '220722_144441_Jacomo';
ksFilePath= ['/media/felix/Internal_1/Data/BevilColor/' filenameP '/kilosorting_laminar/'];
strStitchPath = ['/media/felix/Internal_1/Data/BevilColor/' filenameP '/kilosorting_stitch/'];

spk_times = readNPY([ksFilePath 'spike_times_seconds.npy']);
spk_clusters = readNPY([ksFilePath 'spike_clusters.npy']);
spk_info = tdfread([ksFilePath 'cluster_info.tsv']);
spk_info.Rating = spk_info.r;
    spk_info.group=spk_info.group(:,1:5)

%% temp tools for fixing info arrays
% for tempclust_i=1:length(spk_info_temp.KSLabel)
%     spk_info_temp.KSLabel2(tempclust_i,:)=[spk_info_temp.KSLabel(tempclust_i,:) ' '];
% end
%%
ksFilePath= ['/media/felix/Internal_1/Data/BevilColor/' filenameP '/kilosorting_NForm/'];
NForm_chns = [33,40,46,47,52,53,54,59,65,67,71,81,83,89,90,95,98,102,103,109,112,131,138,139,145,146,152,158]';

spk_times = [spk_times; readNPY([ksFilePath 'spike_times_seconds.npy'])];

nclusts=max(spk_info.cluster_id)+1;
spk_info_temp = tdfread([ksFilePath 'cluster_info.tsv']);

for temp_c=1:length(spk_info_temp.KSLabel)
    spk_info_temp.KSLabel2(temp_c,:)=pad(spk_info_temp.KSLabel(temp_c,:),4,'right');
end
spk_info_temp.KSLabel=spk_info_temp.KSLabel2;

spk_info.cluster_id = [spk_info.cluster_id; spk_info_temp.cluster_id+nclusts];
spk_info.Amplitude = [spk_info.Amplitude; spk_info_temp.Amplitude];
spk_info.ContamPct = [spk_info.ContamPct; spk_info_temp.ContamPct];
spk_info.KSLabel = [spk_info.KSLabel; spk_info_temp.KSLabel];
%spk_info.Rating = [spk_info.Rating; spk_info_temp.Rating];
spk_info.Rating = [spk_info.Rating; spk_info_temp.r];
spk_info.amp = [spk_info.amp; spk_info_temp.amp];
spk_info.ch = [spk_info.ch; NForm_chns(spk_info_temp.ch+1)-1];
spk_info.depth = [spk_info.depth; spk_info_temp.depth];
spk_info.fr = [spk_info.fr; spk_info_temp.fr];
    spk_info_temp.group=spk_info_temp.group(:,1:5)
spk_info.group = [spk_info.group; spk_info_temp.group];
spk_info.n_spikes = [spk_info.n_spikes; spk_info_temp.n_spikes];
spk_info.sh = [spk_info.sh; spk_info_temp.sh];
spk_clusters = [spk_clusters; readNPY([ksFilePath 'spike_clusters.npy'])+nclusts];
%%
ksFilePath= ['/media/felix/Internal_1/Data/BevilColor/' filenameP '/kilosorting_Utah_1to24/'];

spk_times = [spk_times; readNPY([ksFilePath 'spike_times_seconds.npy'])];

nclusts=max(spk_info.cluster_id)+1;
spk_info_temp = tdfread([ksFilePath 'cluster_info.tsv']);

spk_info.cluster_id = [spk_info.cluster_id; spk_info_temp.cluster_id+nclusts];
spk_info.Amplitude = [spk_info.Amplitude; spk_info_temp.Amplitude];
spk_info.ContamPct = [spk_info.ContamPct; spk_info_temp.ContamPct];
spk_info.KSLabel = [spk_info.KSLabel; spk_info_temp.KSLabel];
%spk_info.Rating = [spk_info.Rating; spk_info_temp.Rating];
spk_info.Rating = [spk_info.Rating; spk_info_temp.r];
spk_info.amp = [spk_info.amp; spk_info_temp.amp];
spk_info.ch = [spk_info.ch; spk_info_temp.ch+160];
spk_info.depth = [spk_info.depth; spk_info_temp.depth];
spk_info.fr = [spk_info.fr; spk_info_temp.fr];
    spk_info_temp.group=spk_info_temp.group(:,1:5)
spk_info.group = [spk_info.group; spk_info_temp.group];
spk_info.n_spikes = [spk_info.n_spikes; spk_info_temp.n_spikes];
spk_info.sh = [spk_info.sh; spk_info_temp.sh];
spk_clusters = [spk_clusters; readNPY([ksFilePath 'spike_clusters.npy'])+nclusts];
%%
ksFilePath= ['/media/felix/Internal_1/Data/BevilColor/' filenameP '/kilosorting_Utah_25to48/'];
%ksFilePath= ['/media/felix/Internal_1/Data/BevilColor/' filenameP '/kilosorting_Utah_Ethan2/kilosorting2_Utah_25to48/'];

spk_times = [spk_times; readNPY([ksFilePath 'spike_times_seconds.npy'])];

nclusts=max(spk_info.cluster_id)+1;
spk_info_temp = tdfread([ksFilePath 'cluster_info.tsv']);

for temp_c=1:length(spk_info_temp.KSLabel)
    spk_info_temp.KSLabel2(temp_c,:)=pad(spk_info_temp.KSLabel(temp_c,:),4,'right');
end
spk_info_temp.KSLabel=spk_info_temp.KSLabel2;

spk_info.cluster_id = [spk_info.cluster_id; spk_info_temp.cluster_id+nclusts];
spk_info.Amplitude = [spk_info.Amplitude; spk_info_temp.Amplitude];
spk_info.ContamPct = [spk_info.ContamPct; spk_info_temp.ContamPct];
spk_info.KSLabel = [spk_info.KSLabel; spk_info_temp.KSLabel];
%spk_info.Rating = [spk_info.Rating; spk_info_temp.Rating];
spk_info.Rating = [spk_info.Rating; spk_info_temp.r];
spk_info.amp = [spk_info.amp; spk_info_temp.amp];
spk_info.ch = [spk_info.ch; spk_info_temp.ch+24+160];
spk_info.depth = [spk_info.depth; spk_info_temp.depth];
spk_info.fr = [spk_info.fr; spk_info_temp.fr];
    spk_info_temp.group=spk_info_temp.group(:,1:5)
spk_info.group = [spk_info.group; spk_info_temp.group];
spk_info.n_spikes = [spk_info.n_spikes; spk_info_temp.n_spikes];
spk_info.sh = [spk_info.sh; spk_info_temp.sh];
spk_clusters = [spk_clusters; readNPY([ksFilePath 'spike_clusters.npy'])+nclusts];
%%
ksFilePath= ['/media/felix/Internal_1/Data/BevilColor/' filenameP '/kilosorting_Utah_49to72/'];
%ksFilePath= ['/media/felix/Internal_1/Data/BevilColor/' filenameP '/kilosorting_Utah_Ethan2/kilosorting2_Utah_49to72/'];

spk_times = [spk_times; readNPY([ksFilePath 'spike_times_seconds.npy'])];

nclusts=max(spk_info.cluster_id)+1;
spk_info_temp = tdfread([ksFilePath 'cluster_info.tsv']);

for temp_c=1:length(spk_info_temp.KSLabel)
    spk_info_temp.KSLabel2(temp_c,:)=pad(spk_info_temp.KSLabel(temp_c,:),4,'right');
end
spk_info_temp.KSLabel=spk_info_temp.KSLabel2;

spk_info.cluster_id = [spk_info.cluster_id; spk_info_temp.cluster_id+nclusts];
spk_info.Amplitude = [spk_info.Amplitude; spk_info_temp.Amplitude];
spk_info.ContamPct = [spk_info.ContamPct; spk_info_temp.ContamPct];
spk_info.KSLabel = [spk_info.KSLabel; spk_info_temp.KSLabel];
%spk_info.Rating = [spk_info.Rating; spk_info_temp.Rating];
spk_info.Rating = [spk_info.Rating; spk_info_temp.r];
spk_info.amp = [spk_info.amp; spk_info_temp.amp];
spk_info.ch = [spk_info.ch; spk_info_temp.ch+48+160];
spk_info.depth = [spk_info.depth; spk_info_temp.depth];
spk_info.fr = [spk_info.fr; spk_info_temp.fr];
    spk_info_temp.group=spk_info_temp.group(:,1:5)
spk_info.group = [spk_info.group; spk_info_temp.group];
spk_info.n_spikes = [spk_info.n_spikes; spk_info_temp.n_spikes];
spk_info.sh = [spk_info.sh; spk_info_temp.sh];
spk_clusters = [spk_clusters; readNPY([ksFilePath 'spike_clusters.npy'])+nclusts];
%%
ksFilePath= ['/media/felix/Internal_1/Data/BevilColor/' filenameP '/kilosorting_Utah_73to96/'];
%ksFilePath= ['/media/felix/Internal_1/Data/BevilColor/' filenameP '/kilosorting_Utah_Ethan2/kilosorting2_Utah_73to96/'];

spk_times = [spk_times; readNPY([ksFilePath 'spike_times_seconds.npy'])];

nclusts=max(spk_info.cluster_id)+1;
spk_info_temp = tdfread([ksFilePath 'cluster_info.tsv']);

for temp_c=1:length(spk_info_temp.KSLabel)
    spk_info_temp.KSLabel2(temp_c,:)=pad(spk_info_temp.KSLabel(temp_c,:),4,'right');
end
spk_info_temp.KSLabel=spk_info_temp.KSLabel2;

spk_info.cluster_id = [spk_info.cluster_id; spk_info_temp.cluster_id+nclusts];
spk_info.Amplitude = [spk_info.Amplitude; spk_info_temp.Amplitude];
spk_info.ContamPct = [spk_info.ContamPct; spk_info_temp.ContamPct];
spk_info.KSLabel = [spk_info.KSLabel; spk_info_temp.KSLabel];
%spk_info.Rating = [spk_info.Rating; spk_info_temp.Rating];
spk_info.Rating = [spk_info.Rating; spk_info_temp.r];
%spk_info.Rating = [spk_info.Rating; repmat(NaN, [length(spk_info_temp.KSLabel)],1)];

spk_info.amp = [spk_info.amp; spk_info_temp.amp];
spk_info.ch = [spk_info.ch; spk_info_temp.ch+72+160];
spk_info.depth = [spk_info.depth; spk_info_temp.depth];
spk_info.fr = [spk_info.fr; spk_info_temp.fr];
    spk_info_temp.group=spk_info_temp.group(:,1:5)
spk_info.group = [spk_info.group; spk_info_temp.group];
spk_info.n_spikes = [spk_info.n_spikes; spk_info_temp.n_spikes];
spk_info.sh = [spk_info.sh; spk_info_temp.sh];
spk_clusters = [spk_clusters; readNPY([ksFilePath 'spike_clusters.npy'])+nclusts];

%%
if ~exist(strStitchPath,'dir');
    mkdir(strStitchPath);
end
save([strStitchPath 'KS_stitched.mat'], "spk_clusters", "spk_times", "spk_info")
%%
disp('saving ks_stich - done')
%%