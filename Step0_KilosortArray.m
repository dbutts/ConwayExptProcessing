function [datPath, droptestcheck] = Step0_KilosortArray(dirpath,filenameP,pl2path, outputFolder, opts)
%Step0_KilosortLaminar Packages and then kilosorts laminar probe data
%   Based on custom code by Leor Katz, adapted by Felix Bartsch, 2024

if nargin <3
    preconverted=0;
end
%%
cd(dirpath)
if dirpath(end) ~= filesep
	dirpath = [dirpath filesep];
end

disp('Starting conversion')

opts.commonAverageReferencing = false; % we do not want to do digital rereferencing here, better to do in Plexon
opts.removeArtifacts = false; % enter at own risk!
opts.plotProbeVoltage = false;
opts.extractLfp = false;
opts.outputFolder = outputFolder;

% load in channel map info to load batches of array channels
opts.nChans = length(opts.curchannels);
opts.specificChannels = opts.chInfo.chanMap(opts.curchannels)+opts.ChnOffset;
opts.connected = opts.chInfo.connected(opts.curchannels);
opts.chanMap = opts.chInfo.chanMap(opts.curchannels);
opts.xcoords = opts.chInfo.xcoords(opts.curchannels);
opts.ycoords = opts.chInfo.ycoords(opts.curchannels);
opts.kcoords = opts.chInfo.kcoords(opts.curchannels);
        
%%
% try
%     [fs, n, ts, fn, ~] = plx_ad_v([pl2path filenameP '.pl2'], ['SPKC001'] );
% catch
%     [fs, n, ts, fn, ~] = plx_ad_v([pl2path filenameP '.pl2'], ['SPKC01'] );
% end

try
    [fs, n, ts, fn, ~] = plx_ad_v([pl2path filenameP '.pl2'], ['SPKC001'] );
end
if n<2
    try
        [fs, n, ts, fn, ~] = plx_ad_v([pl2path filenameP '.pl2'], ['SPKC01'] );
    end
end

[~, n_aux, ts_aux, fn_aux, ~] = plx_ad_v([pl2path filenameP '.pl2'], ['AI01'] );
if n<2; error('error reading laminar probe'); end
%%
droptestcheck = [n/40000- n_aux/1000];
%%
if ~opts.preconverted
    %convert pl2 file into dat format
    [samples, datPath] = convertRawToDat_Array([pl2path filenameP '.pl2'],opts);
else
    datPath = [dirpath filenameP '/kilosorting_' opts.ArrayLabel '_' num2str(opts.curchannels(1)) 'to' num2str(opts.curchannels(end)) '/' filenameP '.dat'];
end

disp('Done with conversion')
%% Probe geometry
%[xcoords, ycoords, kcoords] = probeGeometry2coords('linear50', 24); % now
%lives in mastermegafile

%% edit the Megafile, save it into the target directory
%

fs = 40000; %default
nCh = opts.nChans; % = size(samples,1);

ops.connected   = opts.connected;
ops.chanMap     = 1:nCh; %because the indexing has reset following conversion to .dat format
ops.xcoords     = opts.xcoords;
ops.ycoords     = opts.ycoords;
ops.kcoords     = opts.kcoords;

% set Parameters: 
% frequency for high pass filtering (150)
%ops.fshigh = 150; % KS2.5  
ops.fshigh = 300; % KS3: high-pass more aggresively

% minimum firing rate on a "good" channel (0 to skip)
ops.minfr_goodchannels = 0; 

% threshold on projections (like in Kilosort1, can be different for last pass like [10 4])
ops.Th = [9 4]; % KS2.5
%ops.Th = [9 9];  % KS3

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

%% Only used for KS3:
% spatial smoothness constant for registration
ops.sig        = 20;  % default = 20

% blocks for registration. 0 turns it off, 1 does rigid registration. Replaces "datashift" option. 
ops.nblocks    = 10;  % default = 5

% threshold crossings for pre-clustering (in PCA projection space)
ops.ThPre = 8; %defauilt=[8]

%%
step0b_ks_array(datPath, fs, nCh, ops)
%%
disp('DONE with the first sorting pass WOOOOOOOO0')

end

