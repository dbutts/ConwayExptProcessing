function [] = step0b_masterMegaFile_ks3_FB(datPath, fs, nCh, probeGeometry, ops, Nfilt)
% The kilosort repo comes with 3 files that need to be editted:
%       masterFile.m, configFile.m, & createChanMap.m
% 3 files are unnecessarily messy and I'm "1 file to rule them all" kinda 
% guy so I've combined them all into this masterMegaFile.m. 
%
% !!! this m file should copied into each directory you wish to sort. The 
% copying can be done manually or automatically, as I have implemented, for
% example, in batch_convertAndSort.m.
%
% !!! this file has set many of the less-frequently-used parameters of
% kiloSort to my(!) liking, in the 'ops' struct. In all likelihood, you'll
% want the same settings. But it is your responsibility to make sure. 
%
% INPUT:
%   datPath         - full path to binary dat file which you wish to sort
%   fs              - sampling rate (Hz)
%   nCh             - number of channels in your dat file
%   probeGeometry   - string. eg 'linear50, 'dualLinear'... 
%                     see probeGeometry2coords.m for further details
%   Nfilt           - number of template filters you wish to use for sort. 
%                     kiloSort documentation recommends 2-4 times the
%                     number of channels in your file, and multiples of 32.
% OUTPUT:
%   no matlab output. 
%   Function saves a bunch of kiloSort output (m-files and npy files) into
%   the same folder as the .dat file (likely 'kiloSorted' folder, created
%   during conversion)

%% SELECT THY DATASET:
% you can either input the full path to your dat file (datPath) to the
% function or hardcode it here
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dbstop if error

if exist('datPath', 'var')
    [datFolder, datFile, ext] = fileparts(datPath);
    assert(strcmp(ext, '.dat'), ['ERROR: must provide a .dat file. You provided: ' ext])
else
     [datFile, datFolder] = uigetfile('*.dat', 'Select dat file for kiloSorting', pwd);
     datPath = fullfile(datFolder, datFile);
     datFile = datFile(1:end-4);
end
assert(exist(datPath, 'file')~=0, 'ERROR: can''t find your dat file');
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%% Define paths to your dataset:
addPathsForSpikeSorting_ks3;

%% 'createChanMap.m' section:

connected       = true(nCh, 1);
chanMap         = 1:nCh;
chanMap0ind     = chanMap - 1;

[xcoords, ycoords, kcoords] = probeGeometry2coords(probeGeometry, nCh);

% save:
save(fullfile(datFolder, 'chanMap.mat'), ...
    'chanMap','connected', 'xcoords', 'ycoords', 'kcoords', 'chanMap0ind', 'fs', 'probeGeometry')

%% 'standardConfig.m' section:

% lnk: I added this here:
ops.root                = datFolder; 
% this section is taken from the \configFiles\StandardConfig_MOVEME.m file.
ops.chanMap             = fullfile(datFolder, 'chanMap.mat');
% sample rate
ops.fs = fs;  



%% danger, changing these settings can lead to fatal errors
% options for determining PCs
ops.spkTh           = -6;      % spike threshold in standard deviations (-6)
ops.reorder         = 1;       % whether to reorder batches for drift correction. 
ops.nskip           = 25;  % how many batches to skip for determining spike PCs

ops.GPU                 = 1; % has to be 1, no CPU version yet, sorry

% ops.Nfilt               = 1024; % max number of clusters
ops.nfilt_factor        = 4; % max number of clusters per good channel (even temporary ones)
ops.ntbuff              = 64;    % samples of symmetrical buffer for whitening and spike detection
ops.NT                  = 2*64*1024+ ops.ntbuff; % must be multiple of 32 + ntbuff. This is the batch size (try decreasing if out of memory). 

ops.whiteningRange      = 24; % number of channels to use for whitening each channel
ops.nSkipCov            = 25; % compute whitening matrix from every N-th batch
ops.scaleproc           = 200;   % int16 scaling of whitened data
ops.nPCs                = 3; % how many PCs to project the spikes into
ops.useRAM              = 0; % not yet available

save(fullfile(ops.root,'ops.mat'), 'ops')

%% master_kilosort section:
rootZ = datFolder; % the raw data binary file is in this folder

[~, hostName] = system('hostname');
if contains(hostName, 'IT', 'IgnoreCase', 1) 
    rootH = '/home/felix/kilodata'; % path to temporary binary file (same size as data, should be on fast SSD)
elseif contains(hostName, 'mt', 'IgnoreCase', 1) 
    rootH = '/home/fellixbartsch/kilodata'; % path to temporary binary file (same size as data, should be on fast SSD)
end

% pathToYourConfigFile = 'D:\GitHub\KiloSort2\configFiles'; % take from Github folder and put it somewhere else (together with the master_file)
% chanMapFile = 'neuropixPhase3A_kilosortChanMap.mat';

ops.trange      = [0 Inf]; % time range to sort
ops.NchanTOT    = nCh; % total number of channels in your recording

% run(fullfile(pathToYourConfigFile, 'configFile384.m'))
ops.fproc       = fullfile(rootH, 'temp_wh.dat'); % proc file on a fast SSD
% ops.chanMap = fullfile(pathToYourConfigFile, chanMapFile);

%% this block runs all the steps of the algorithm
fprintf('Looking for data inside %s \n', rootZ)

% is there a channel map file in this folder?
fs = dir(fullfile(rootZ, 'chan*.mat'));
if ~isempty(fs)
    ops.chanMap = fullfile(rootZ, fs(1).name);
end

% find the binary file
fs          = [dir(fullfile(rootZ, '*.bin')) dir(fullfile(rootZ, '*.dat'))];
ops.fbinary = fullfile(rootZ, fs(1).name);
%% KS2.5 approach: 
% % preprocess data to create temp_wh.dat
% rez = preprocessDataSub(ops);
% % time-reordering as a function of drift
% rez = clusterSingleBatches(rez);
% % saving here is a good idea, because the rest can be resumed after loading rez
% save(fullfile(rootZ, 'rez.mat'), 'rez', '-v7.3');
% % main tracking and template matching algorithm
% rez = learnAndSolve8b(rez);
% % final merges
% rez = find_merges(rez, 1);
% % final splits by SVD
% rez = splitAllClusters(rez, 1);
% % final splits by amplitudes
% rez = splitAllClusters(rez, 0);
% % decide on cutoff
% rez = set_cutoff(rez);
% fprintf('found %d good units \n', sum(rez.good>0))
% 
% % write to Phy
% fprintf('Saving results to Phy  \n')
% rezToPhy(rez, rootZ);
%% KS3 appproach:
rez                = preprocessDataSub(ops);
rez                = datashift2(rez, 1);

[rez, st3, tF]     = extract_spikes(rez);
rez                = template_learning(rez, tF, st3);
[rez, st3, tF]     = trackAndSort(rez);
rez                = final_clustering(rez, tF, st3);
rez                = find_merges(rez, 1);

fprintf('found %d good units \n', sum(rez.good>0))
fprintf('Saving results to Phy  \n')
rezToPhy2(rez, rootZ);
%%
% if you want to save the results to a Matlab file...

% discard features in final rez file (too slow to save)
rez.cProj = [];
rez.cProjPC = [];

% save final results as rez2
fprintf('Saving final results in rez2  \n')
fname = fullfile(rootZ, 'rez2.mat');
save(fname, 'rez', '-v7.3');

%% below is from version 1 kiloSrot
% % % %% sort that sweet sweet data:
% % % % taken from the masterFile.m
% % % 
% % % if ops.GPU     
% % %     gpuDevice(1); % initialize GPU (will erase any existing GPU arrays)
% % % end
% % % 
% % % tic;
% % % [rez, DATA, uproj] = preprocessData(ops); % preprocess data and extract spikes for initialization
% % % rez                = fitTemplates(rez, DATA, uproj);  % fit templates iteratively
% % % rez                = fullMPMU(rez, DATA);% extract final spike times (overlapping extraction)
% % % 
% % % % save sort results
% % % save(fullfile(ops.root,  'rez.mat'), 'rez', '-v7.3');
% % % rezToPhy(rez, ops.root);
% % % 
% % % % AutoMerge:
% % % % rez2Phy will use clusters from the new 5th column of st3 if you run this
% % % if isfield(ops, 'performAutoMerge') && ops.performAutoMerge
% % %     rez = merge_posthoc2(rez);
% % %     % save post autoMerge results in root folder:
% % %     autoMergeFolder = fullfile(datFolder, 'autoMerged');
% % %     mkdir(autoMergeFolder)
% % %     save(fullfile(autoMergeFolder,  'rez.mat'), 'rez', '-v7.3');
% % %     rezToPhy(rez, autoMergeFolder);
% % % end
% % % 
% % % % remove temporary file
% % % delete(ops.fproc);

%% sampsToSecsMap:
% load "samples to seconds map" (created during 'convertRawToDat.m'), 
% transform spike times in samples to spike times in seconds and store.


disp('Loading samples-to-seconds map.');
load(fullfile(ops.root, 'sampsToSecsMap.mat'));
disp('Map loaded - transforming spike times from samples to seconds.')
spikeTimesSamples = readNPY(fullfile(ops.root, ...
    'spike_times.npy'));
spikeTimesSeconds = sampsToSecsMap(spikeTimesSamples);
disp('Saving spike times in seconds to .npy file.')
writeNPY(spikeTimesSeconds, fullfile(ops.root, ...
    'spike_times_seconds.npy'));

%% ding!
try
    load gong
    soundsc(y, Fs)
%     evalc('soundsc(double(aiffread(''/System/Library/Sounds/Glass.aiff''))'',50000)');
catch
end

