%%
disp('Loading samples-to-seconds map.');

ops.root = ['/media/felix/Internal_1/Data/BevilColor/220914_160534_Jacomo/kilosorting_NForm/'];

load(fullfile([ops.root  'sampsToSecsMap.mat']));
disp('Map loaded - transforming spike times from samples to seconds.')
%%
spikeTimesSamples = readNPY(fullfile(ops.root, ...
    'spike_times.npy'));
spikeTimesSeconds = sampsToSecsMap(spikeTimesSamples);
disp('Saving spike times in seconds to .npy file.')
writeNPY(spikeTimesSeconds, fullfile(ops.root, ...
    'spike_times_seconds.npy'));
disp('done')
%%