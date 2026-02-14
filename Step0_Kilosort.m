function [datPath, droptestcheck] = Step0_Kilosort(dirpath,filenameP,pl2path, opts)
%Step0_KilosortLaminar Packages and then kilosorts laminar probe data
%   Based on custom code by Leor Katz, adapted by Felix Bartsch, 2024

plexon_fname = fullfile(pl2path, [filenameP '.pl2']);
if nargin <3
    preconverted=0;
end

fs = 40000; %default
try
    nCh = opts.nChans; % = size(samples,1);
catch
    opts.nChans = length(opts.curchannels);
    nCh = opts.nChans;
end
%% ops fields

% frequency for high pass filtering (150)
ops.fshigh = 300;

% minimum firing rate on a "good" channel (0 to skip)
ops.minfr_goodchannels = 0;

% how important is the amplitude penalty (like in Kilosort1, 0 means not used, 10 is average, 50 is a lot)
ops.lam = 10;

% splitting a cluster at the end requires at least this much isolation for each sub-cluster (max = 1)
ops.AUCsplit = 0.9;

% minimum spike rate (Hz), if a cluster falls below this for too long it gets removed
ops.minFR = 0;

% number of samples to average over (annealed from first to second value)
ops.momentum = [20 400];

% spatial constant in um for computing residual variance of spike
ops.sigmaMask  = 30;

% Only used for KS3
% spatial smoothness constant for registration
ops.sig        = 20;

% blocks for registration. 0 turns it off, 1 does rigid registration. Replaces "datashift" option.
ops.nblocks    = 10;  % default = 5

% threshold crossings for pre-clustering (in PCA projection space)
ops.ThPre = 8; %defauilt=[8]

% are we dealing with laminar or utah array?
if contains(opts.ArrayLabel, 'lam', 'IgnoreCase', true) % laminar

    outputFolder = fullfile(dirpath, filenameP, 'kilosorting_laminar/');
    % threshold on projections (like in Kilosort1, can be different for last pass like [10 4])

    ops.Th = [9 5];
    probeGeometry = ['linear' num2str(opts.arraySpacing(2))];
    probeGeometryFields = {'xcoords', 'ycoords', 'kcoords'};

    if isfield(opts, 'chInfo') && all(isfield(opts.chInfo, probeGeometryFields))
            ops = loadChInfoIntoOps(opts, ops);
            ops.connected = true(nCh, 1);
    else
            % make coordinates

            opts.chInfo.connected = true(nCh, 1);
            opts.chInfo.chanMap = 1:nCh;
            opts.chInfo.chanMap0ind  = opts.chInfo.chanMap - 1;
            opts.chInfo.xcoords = zeros(nCh,1);
            opts.chInfo.ycoords = opts.arraySpacing(2).*linspace(nCh-1, 0, nCh);
            opts.chInfo.kcoords = ones(nCh,1);
            opts.specificChannels = opts.chInfo.chanMap(opts.curchannels)+opts.ChnOffset;
            opts.chInfo.fs = fs;

            ops = loadChInfoIntoOps(opts, ops);
            


    end
  
elseif contains(opts.ArrayLabel, 'UT', 'IgnoreCase', true) % array

    probeGeometry = 'array';
    outputFolder = fullfile(dirpath, filenameP, ['kilosorting_' opts.ArrayLabel '_' num2str(min(opts.curchannels)) 'to' num2str(max(opts.curchannels))]);

    % threshold on projections (like in Kilosort1, can be different for last pass like [10 4])
    ops.Th = [9 4]; %MJG: no clue why this is different from lam case

    ops = loadChInfoIntoOps(opts, ops);
    opts.specificChannels = opts.chInfo.chanMap(opts.curchannels)+opts.ChnOffset;

end

%% opts fields
opts.commonAverageReferencing = false; % we do not want to do digital rereferencing here, better to do in Plexon
opts.removeArtifacts = false; % enter at own risk!
opts.plotProbeVoltage = false;
opts.extractLfp = false;
opts.outputFolder = outputFolder;

%%
%cd(dirpath) % may need to uncomment?
if dirpath(end) ~= filesep
    dirpath = [dirpath filesep];
end

%opts.specificChannels = [];  % user can select which plexon channels to use for
%                               conversion. remember, this must be in plexon-numbering.
%%
pl2 = PL2ReadFileIndex(plexon_fname);
numDigitsInLastSpkChan = ceil(log10(length(pl2.SpikeChannels)));

[fs, n, ts, fn, ~] = plx_ad_v(plexon_fname,...
    ['SPKC' num2str(1, ['%0' num2str(numDigitsInLastSpkChan) '.f'])]);

[fs_aux, n_aux, ts_aux, fn_aux, ~] = plx_ad_v(plexon_fname, ['AI01'] );

if n<2
    error('error reading laminar probe')
end
%%
droptestcheck = n_aux/fs_aux - n/fs;
%%
disp('Starting conversion')
if ~opts.preconverted
    %convert pl2 file into dat format
    %[samples, datPath] = convertRawToDat_Laminar([pl2path filenameP filesep filenameP '.pl2'],opts);
    if contains(opts.ArrayLabel, 'UT', 'IgnoreCase', true)
        [samples, datPath] = convertRawToDat_Array(plexon_fname,opts);
    elseif contains(opts.ArrayLabel, 'lam', 'IgnoreCase', true)
        [samples, datPath] = convertRawToDat_Array(plexon_fname, opts);
    end
else
    datPath = fullfile(outputFolder, [filenameP '.dat']);
end

disp('Done with conversion')

%% edit the Megafile, save it into the target directory
%

if exist('datPath', 'var')
    [datFolder, datFile, ext] = fileparts(datPath);
    assert(strcmp(ext, '.dat'), ['ERROR: must provide a .dat file. You provided: ' ext])
else
    [datFile, datFolder] = uigetfile('*.dat', 'Select dat file for kiloSorting', pwd);
    datPath = fullfile(datFolder, datFile);
    datFile = datFile(1:end-4);
end
assert(exist(datPath, 'file')~=0, 'ERROR: can''t find your dat file');


%% Define paths to your dataset:
%addPathsForSpikeSorting; % MJG: dont think I need this




%% danger, changing these settings can lead to fatal errors
% options for determining PCs
ops.spkTh           = -6;      % spike threshold in standard deviations (-6)
ops.reorder         = 1;       % whether to reorder batches for drift correction.
ops.nskip           = 25;  % how many batches to skip for determining spike PCs

ops.GPU                 = 1; % has to be 1, no CPU version yet, sorry

% ops.Nfilt               = 1024; % max number of clusters
ops.nfilt_factor        = 4; % max number of clusters per good channel (even temporary ones)
ops.ntbuff              = 64;    % samples of symmetrical buffer for whitening and spike detection
ops.NT                  = 1*64*1024+ ops.ntbuff; % must be multiple of 32 + ntbuff. This is the batch size (try decreasing if out of memory). Originally 2*64*1024
%%

ops.whiteningRange      = 24; % number of channels to use for whitening each channel
ops.nSkipCov            = 25; % compute whitening matrix from every N-th batch
ops.scaleproc           = 200;   % int16 scaling of whitened data
ops.nPCs                = 3; % how many PCs to project the spikes into
ops.useRAM              = 0; % not yet available

connected       = ops.connected;
chanMap         = ops.chanMap-min(ops.chanMap)+1;
chanMap0ind     = chanMap - 1;

xcoords         = ops.xcoords;
ycoords         = ops.ycoords;
kcoords         = ops.kcoords;

save(fullfile(datFolder, 'chanMap.mat'), ...
    'chanMap','connected', 'xcoords', 'ycoords', 'kcoords', 'chanMap0ind', 'fs', 'probeGeometry');

%% 'standardConfig.m' section:

% lnk: I added this here:
ops.root                = datFolder;
% this section is taken from the \configFiles\StandardConfig_MOVEME.m file.
ops.chanMap             = fullfile(datFolder, 'chanMap.mat');
% sample rate
ops.fs = fs;

save(fullfile(ops.root,'ops.mat'), 'ops')

%% master_kilosort section:
rootZ = datFolder; % the raw data binary file is in this folder

[~, hostName] = system('hostname');
if contains(hostName, 'IT', 'IgnoreCase', 1)
    rootH = '/home/felix/kilodata'; % path to temporary binary file (same size as data, should be on fast SSD)
elseif contains(hostName, 'mt', 'IgnoreCase', 1)
    rootH = '/home/fellixbartsch/kilodata'; % path to temporary binary file (same size as data, should be on fast SSD)
elseif contains(hostName, 'conwaylab', 'IgnoreCase', 1)
    rootH = '/home/conwaylab/kilodata';
elseif contains(hostName, 'dl', 'IgnoreCase', 1)
    rootH = '/home/bizon/kilodata';
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

% preprocess data to create temp_wh.dat
rez = preprocessDataSub(ops);

% time-reordering as a function of drift
rez = clusterSingleBatches(rez);

%
% saving here is a good idea, because the rest can be resumed after loading rez
save(fullfile(rootZ, 'rez.mat'), 'rez', '-v7.3');

% main tracking and template matching algorithm
rez = learnAndSolve8b(rez);

% final merges
rez = find_merges(rez, 1);

% final splits by SVD
rez = splitAllClusters(rez, 1); %commented out to test memory errors

% final splits by amplitudes
rez = splitAllClusters(rez, 0);

% decide on cutoff
rez = set_cutoff(rez);

fprintf('found %d good units \n', sum(rez.good>0))

% write to Phy
fprintf('Saving results to Phy  \n')
rezToPhy(rez, rootZ);

% if you want to save the results to a Matlab file...

% discard features in final rez file (too slow to save)
rez.cProj = [];
rez.cProjPC = [];

% save final results as rez2
fprintf('Saving final results in rez2  \n')
fname = fullfile(rootZ, 'rez2.mat');
save(fname, 'rez', '-v7.3');

disp('Loading samples-to-seconds map.');
load(fullfile(ops.root, 'sampsToSecsMap.mat'));
disp('Map loaded - transforming spike times from samples to seconds.')
spikeTimesSamples = readNPY(fullfile(ops.root, ...
    'spike_times.npy'));
spikeTimesSeconds = sampsToSecsMap(spikeTimesSamples);
disp('Saving spike times in seconds to .npy file.')
writeNPY(spikeTimesSeconds, fullfile(ops.root, ...
    'spike_times_seconds.npy'));
%%
disp('DONE with the first sorting pass WOOOOOOOO0')

function ops = loadChInfoIntoOps(opts, ops)

if nargin == 1
    ops = struct;
end

if ~isfield(opts, 'curchannels')
    opts.curchannels = 1:opts.nChans;
end
ops.connected = opts.chInfo.connected(opts.curchannels);
ops.chanMap = 1:opts.nChans;
ops.xcoords = opts.chInfo.xcoords(opts.curchannels);
ops.ycoords = opts.chInfo.ycoords(opts.curchannels);
ops.kcoords =  opts.chInfo.kcoords(opts.curchannels);


