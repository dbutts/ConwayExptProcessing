%% Kilo pipeline
clear all
%%
filenameP = '230510_141725_Jacomo'; iso_SUs=[];

%%
cd('C:\SpkSort2023\AAActiveData')
outputFolder = ['C:\SpkSort2023\AAActiveData\' filenameP '\kilosorting_laminar/'];

disp('Starting conversion')

opts.commonAverageReferencing = false; % TODO: set to true as experiment later
opts.removeArtifacts = false; % enter at own risk!
opts.plotProbeVoltage = false;
opts.extractLfp = false;
opts.outputFolder = outputFolder;

opts.nChans = 24;
opts.ChnOffset=0;
opts.batch_size = 300000000;
%opts.specificChannels = [];  % user can select which plexon channels to use for
%                               conversion. remember, this must be in
%                               plexon-numbering, eg SPKC1 is usually ch num 65.
%
%%
[fs, n, ts, fn, ~] = plx_ad_v([filenameP '.pl2'], ['SPKC001'] );
[~, n_aux, ts_aux, fn_aux, ~] = plx_ad_v([filenameP '.pl2'], ['AI01'] );
if n<2; error('error reading laminar probe'); end
%%
droptestcheck = [n/40000- n_aux/1000]
%%
[samples, datPath] = convertRawToDat_Laminar([filenameP '.pl2'],opts);%uncommented for spike sorting a session for the first time

%datPath = ['C:\SpkSort2023\AAActiveData\' filenameP
%'\kilosorting_laminar\' filenameP '.dat'];%uncommented for spike sorting a session for subsequent times

disp('Done with conversion')
%% Probe geometry
%[xcoords, ycoords, kcoords] = probeGeometry2coords('linear50', 24); % now
%lives in mastermegafile

%% edit the Megafile, save it into the target directory
%edit masterMegaFile_base_FB

%
% if ~exist(['/media/felix/Internal_1/Data/BevilColor/' filenameP '/kilosorting/kiloSorted2/'],'dir')
%     mkdir(['/media/felix/Internal_1/Data/BevilColor/' filenameP '/kilosorting/kiloSorted2/']);
% end
fs = 40000; %from above
nCh = opts.nChans; % = size(samples,1);

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
ops.AUCsplit = 0.5; 

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
step0b_masterMegaFile_base_FBv2(datPath, fs, nCh, ['linear50'], ops, 128)
%%
disp('DONE with the first sorting pass WOOOOO')
%%
